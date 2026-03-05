#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }
    out[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

options <- parse_args(args)

required <- c("input", "output", "groupby", "initial-color")
missing <- required[!required %in% names(options)]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo_root, "R", "helpers.R"))
source(file.path(repo_root, "R", "source.R"))
source(file.path(repo_root, "R", "payload.R"))
source(file.path(repo_root, "R", "export.R"))

split_csv <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return(NULL)
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

format_percentage <- function(x) {
  sprintf("%.1f%%", 100 * x)
}

format_column_coverage <- function(obs, columns, max_columns = 8L) {
  columns <- intersect(as.character(columns %||% character()), names(obs))
  if (length(columns) == 0L) {
    return("<none>")
  }

  coverage <- vapply(columns, function(column_name) {
    fraction_non_missing_values(obs[[column_name]])
  }, numeric(1))
  order_idx <- order(coverage, columns, decreasing = TRUE)
  columns <- columns[order_idx]
  coverage <- coverage[order_idx]

  limit <- min(length(columns), max_columns)
  preview <- sprintf("%s %s", columns[seq_len(limit)], format_percentage(coverage[seq_len(limit)]))
  paste(preview, collapse = ", ")
}

report_metadata_merge <- function(obj) {
  report <- extract_metadata_merge_report(obj)
  if (is.null(report)) {
    return(invisible(NULL))
  }

  obs <- extract_metadata_table(obj)
  cat(
    "Metadata overlap: ",
    report$overlap_rows,
    "/",
    report$primary_rows,
    " rows (",
    format_percentage(report$overlap_fraction),
    ")\n",
    sep = ""
  )

  coverage_table <- report$column_coverage
  partial <- coverage_table[coverage_table$coverage < 0.999999, , drop = FALSE]
  if (nrow(partial) > 0L) {
    cat(
      "Partially annotated merged columns: ",
      format_column_coverage(obs, partial$column),
      "\n",
      sep = ""
    )
  }

  invisible(report)
}

warn_low_coverage_selected_colors <- function(obj, selected_columns, min_coverage = 0.9) {
  report <- extract_metadata_merge_report(obj)
  if (is.null(report)) {
    return(invisible(NULL))
  }

  low_coverage_columns <- report$column_coverage$column[
    report$column_coverage$coverage < min_coverage
  ]
  selected_low_coverage <- intersect(as.character(selected_columns %||% character()), low_coverage_columns)
  if (length(selected_low_coverage) > 0L) {
    obs <- extract_metadata_table(obj)
    cat(
      "Selected low-coverage color columns: ",
      format_column_coverage(obs, selected_low_coverage),
      "\n",
      sep = ""
    )
  }

  invisible(report)
}

prepared_input <- prepare_karospace_input(
  input = options[["input"]],
  metadata_input = options[["metadata-input"]],
  metadata_input_columns = split_csv(options[["metadata-input-columns"]]),
  metadata_prefix = options[["metadata-prefix"]]
)
report_metadata_merge(prepared_input)
warn_low_coverage_selected_colors(
  obj = prepared_input,
  selected_columns = c(options[["initial-color"]], split_csv(options[["additional-colors"]]))
)

export_karospace_viewer(
  input = prepared_input,
  output_path = options[["output"]],
  groupby = options[["groupby"]],
  initial_color = options[["initial-color"]],
  additional_colors = split_csv(options[["additional-colors"]]),
  genes = split_csv(options[["genes"]]),
  assay = options[["assay"]],
  metadata_input = NULL,
  metadata_input_columns = NULL,
  metadata_prefix = NULL,
  neighbor_mode = options[["neighbor-mode"]] %||% "spatial",
  neighbor_graph = options[["neighbor-graph"]],
  neighbor_k = as.integer(options[["neighbor-k"]] %||% 6L),
  metadata_columns = split_csv(options[["metadata-columns"]]),
  outline_by = options[["outline-by"]],
  title = options[["title"]] %||% "KaroSpace",
  theme = options[["theme"]] %||% "light",
  min_panel_size = as.numeric(options[["min-panel-size"]] %||% 150),
  spot_size = as.numeric(options[["spot-size"]] %||% 2)
)

cat("Wrote viewer to ", normalizePath(options[["output"]], mustWork = FALSE), "\n", sep = "")
