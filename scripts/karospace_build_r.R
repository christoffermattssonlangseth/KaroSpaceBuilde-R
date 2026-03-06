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

    next_is_value <- i < length(args) && !startsWith(args[[i + 1L]], "--")
    if (next_is_value) {
      out[[substring(key, 3L)]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      out[[substring(key, 3L)]] <- TRUE
      i <- i + 1L
    }
  }
  out
}

options <- parse_args(args)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo_root, "R", "helpers.R"))
source(file.path(repo_root, "R", "source.R"))
source(file.path(repo_root, "R", "payload.R"))
source(file.path(repo_root, "R", "export.R"))

if (isTRUE(options$help) || length(options) == 0L) {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/karospace_build_r.R --config path/to/build.json",
      "Rscript scripts/karospace_build_r.R --input input.rds --output viewer.html --groupby sample_id --initial-color cell_type",
      "",
      "Optional overrides:",
      "[--config build.json] [--write-config merged-build.json] [--additional-colors course,condition]",
      "[--metadata-input metadata.rds] [--metadata-input-columns col1,col2] [--metadata-prefix ext_]",
      "[--assay SCT] [--genes GENE1,GENE2] [--top-genes 200] [--lightweight]",
      "[--neighbor-mode spatial] [--neighbor-graph integrated_snn] [--neighbor-k 6]",
      "[--metadata-columns course,condition] [--outline-by sample_name]",
      "[--marker-genes-groupby auto] [--marker-genes-top-n 20] [--marker-test mean_diff|wilcoxon]",
      "[--interaction-markers-groupby cell_type] [--interaction-markers-top-targets 8] [--interaction-markers-top-genes 12]",
      "[--interaction-markers-min-cells 30] [--interaction-markers-min-neighbors 1]",
      "[--neighbor-stats-permutations 0] [--neighbor-stats-seed 42]",
      "[--title MyViewer] [--theme light] [--min-panel-size 150] [--spot-size 2]",
      sep = "\n"
    ),
    "\n"
  )
  quit(save = "no", status = 0L)
}

split_csv <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return(NULL)
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

parse_bool_option <- function(value, default = FALSE) {
  if (is.null(value)) {
    return(default)
  }
  normalized <- tolower(trimws(as.character(value)[[1]]))
  if (normalized %in% c("true", "t", "1", "yes", "y", "on")) {
    return(TRUE)
  }
  if (normalized %in% c("false", "f", "0", "no", "n", "off")) {
    return(FALSE)
  }
  stop("Could not parse boolean option value: ", value)
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

config_source_path <- if (!is.null(options$config)) {
  normalizePath(options$config, mustWork = TRUE)
} else {
  NULL
}
config <- if (!is.null(config_source_path)) {
  read_karospace_build_config(config_source_path)
} else {
  list(version = 1L)
}

if (!is.null(options$input)) {
  config$input <- options$input
}
if (!is.null(options$output)) {
  config$output <- options$output
}
if (!is.null(options$groupby)) {
  config$groupby <- options$groupby
}
if (!is.null(options[["initial-color"]])) {
  config$initial_color <- options[["initial-color"]]
}
if (!is.null(options[["additional-colors"]])) {
  config$additional_colors <- split_csv(options[["additional-colors"]])
}
if (!is.null(options[["metadata-input"]])) {
  config$metadata_input <- options[["metadata-input"]]
}
if (!is.null(options[["metadata-input-columns"]])) {
  config$metadata_input_columns <- split_csv(options[["metadata-input-columns"]])
}
if (!is.null(options[["metadata-prefix"]])) {
  config$metadata_prefix <- options[["metadata-prefix"]]
}
if (!is.null(options[["assay"]])) {
  config$assay <- options[["assay"]]
}
if (!is.null(options[["genes"]])) {
  config$genes <- split_csv(options[["genes"]])
}
if (!is.null(options[["top-genes"]])) {
  config$top_genes_n <- suppressWarnings(as.integer(options[["top-genes"]]))
}
if (!is.null(options[["lightweight"]])) {
  config$lightweight <- parse_bool_option(options[["lightweight"]], default = FALSE)
}
if (!is.null(options[["neighbor-mode"]])) {
  config$neighbor_mode <- options[["neighbor-mode"]]
}
if (!is.null(options[["neighbor-graph"]])) {
  config$neighbor_graph <- options[["neighbor-graph"]]
}
if (!is.null(options[["neighbor-k"]])) {
  config$neighbor_k <- suppressWarnings(as.integer(options[["neighbor-k"]]))
}
if (!is.null(options[["metadata-columns"]])) {
  config$metadata_columns <- split_csv(options[["metadata-columns"]])
}
if (!is.null(options[["outline-by"]])) {
  config$outline_by <- options[["outline-by"]]
}
if (!is.null(options[["marker-genes-groupby"]])) {
  config$marker_genes_groupby <- split_csv(options[["marker-genes-groupby"]])
}
if (!is.null(options[["marker-genes-top-n"]])) {
  config$marker_genes_top_n <- suppressWarnings(as.integer(options[["marker-genes-top-n"]]))
}
if (!is.null(options[["interaction-markers-groupby"]])) {
  config$interaction_markers_groupby <- split_csv(options[["interaction-markers-groupby"]])
}
if (!is.null(options[["interaction-markers-top-targets"]])) {
  config$interaction_markers_top_targets <- suppressWarnings(as.integer(options[["interaction-markers-top-targets"]]))
}
if (!is.null(options[["interaction-markers-top-genes"]])) {
  config$interaction_markers_top_genes <- suppressWarnings(as.integer(options[["interaction-markers-top-genes"]]))
}
if (!is.null(options[["interaction-markers-min-cells"]])) {
  config$interaction_markers_min_cells <- suppressWarnings(as.integer(options[["interaction-markers-min-cells"]]))
}
if (!is.null(options[["interaction-markers-min-neighbors"]])) {
  config$interaction_markers_min_neighbors <- suppressWarnings(as.integer(options[["interaction-markers-min-neighbors"]]))
}
if (!is.null(options[["marker-test"]])) {
  config$marker_test <- options[["marker-test"]]
}
if (!is.null(options[["neighbor-stats-permutations"]])) {
  config$neighbor_stats_permutations <- suppressWarnings(as.integer(options[["neighbor-stats-permutations"]]))
}
if (!is.null(options[["neighbor-stats-seed"]])) {
  config$neighbor_stats_seed <- suppressWarnings(as.integer(options[["neighbor-stats-seed"]]))
}
if (!is.null(options[["title"]])) {
  config$title <- options[["title"]]
}
if (!is.null(options[["theme"]])) {
  config$theme <- options[["theme"]]
}
if (!is.null(options[["min-panel-size"]])) {
  config$min_panel_size <- suppressWarnings(as.numeric(options[["min-panel-size"]]))
}
if (!is.null(options[["spot-size"]])) {
  config$spot_size <- suppressWarnings(as.numeric(options[["spot-size"]]))
}

resolved_config <- normalize_karospace_build_config(
  config = config,
  config_path = config_source_path,
  require_required = TRUE,
  resolve_paths = TRUE
)

prepared_input <- prepare_karospace_input(
  input = resolved_config$input,
  metadata_input = resolved_config$metadata_input,
  metadata_input_columns = resolved_config$metadata_input_columns,
  metadata_prefix = resolved_config$metadata_prefix
)
report_metadata_merge(prepared_input)
warn_low_coverage_selected_colors(
  obj = prepared_input,
  selected_columns = c(resolved_config$initial_color, resolved_config$additional_colors)
)

if (!is.null(options[["write-config"]])) {
  write_karospace_build_config(
    config = resolved_config,
    path = options[["write-config"]]
  )
  cat(
    "Wrote config to ",
    normalizePath(options[["write-config"]], mustWork = FALSE),
    "\n",
    sep = ""
  )
}

export_karospace_viewer(
  input = prepared_input,
  output_path = resolved_config$output,
  groupby = resolved_config$groupby,
  initial_color = resolved_config$initial_color,
  additional_colors = resolved_config$additional_colors,
  genes = resolved_config$genes,
  top_genes_n = resolved_config$top_genes_n,
  assay = resolved_config$assay,
  metadata_input = NULL,
  metadata_input_columns = NULL,
  metadata_prefix = NULL,
  neighbor_mode = resolved_config$neighbor_mode,
  neighbor_graph = resolved_config$neighbor_graph,
  neighbor_k = resolved_config$neighbor_k,
  metadata_columns = resolved_config$metadata_columns,
  outline_by = resolved_config$outline_by,
  lightweight = resolved_config$lightweight,
  marker_genes_groupby = resolved_config$marker_genes_groupby,
  marker_genes_top_n = resolved_config$marker_genes_top_n,
  interaction_markers_groupby = resolved_config$interaction_markers_groupby,
  interaction_markers_top_targets = resolved_config$interaction_markers_top_targets,
  interaction_markers_top_genes = resolved_config$interaction_markers_top_genes,
  interaction_markers_min_cells = resolved_config$interaction_markers_min_cells,
  interaction_markers_min_neighbors = resolved_config$interaction_markers_min_neighbors,
  marker_test = resolved_config$marker_test,
  neighbor_stats_permutations = resolved_config$neighbor_stats_permutations,
  neighbor_stats_seed = resolved_config$neighbor_stats_seed,
  title = resolved_config$title %||% "KaroSpace",
  theme = resolved_config$theme,
  min_panel_size = resolved_config$min_panel_size,
  spot_size = resolved_config$spot_size
)

cat("Wrote viewer to ", normalizePath(resolved_config$output, mustWork = FALSE), "\n", sep = "")
