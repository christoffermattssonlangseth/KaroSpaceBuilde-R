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

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(repo_root, "R", "helpers.R"))
source(file.path(repo_root, "R", "source.R"))
source(file.path(repo_root, "R", "payload.R"))
source(file.path(repo_root, "R", "export.R"))

options <- parse_args(args)

if (isTRUE(options$help) || is.null(options$input)) {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/example_export.R --input path/to/object.rds [--output viewer.html]",
      "[--groupby sample_id] [--initial-color cell_type] [--additional-colors course,condition]",
      "[--title MyViewer] [--theme light] [--inspect]",
      sep = "\n"
    ),
    "\n"
  )
  quit(save = "no", status = if (is.null(options$input)) 1L else 0L)
}

split_csv <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return(NULL)
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

extract_obs <- function(x) {
  if (inherits(x, "Seurat")) {
    return(x@meta.data)
  }

  if (inherits(x, "SpatialExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SpatialExperiment input requires SummarizedExperiment to inspect colData.")
    }
    return(as.data.frame(SummarizedExperiment::colData(x)))
  }

  if (inherits(x, "SingleCellExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SingleCellExperiment input requires SummarizedExperiment to inspect colData.")
    }
    return(as.data.frame(SummarizedExperiment::colData(x)))
  }

  if (is.list(x) && is.data.frame(x$obs)) {
    return(x$obs)
  }

  stop(
    "Could not inspect the input. Supported inputs are list with obs, SingleCellExperiment, and SpatialExperiment."
  )
}

is_color_candidate <- function(column) {
  is.factor(column) || is.character(column) || is.logical(column) || is.numeric(column)
}

is_groupby_candidate <- function(column) {
  values <- unique_non_missing(column)
  n_unique <- length(values)
  if (is.factor(column) || is.character(column) || is.logical(column)) {
    return(n_unique > 1L && n_unique <= max(500L, floor(length(column) * 0.95)))
  }

  if (is.numeric(column)) {
    is_whole <- all(is.na(column) | abs(column - round(column)) < 1e-8)
    return(is_whole && n_unique > 1L && n_unique <= min(100L, floor(length(column) * 0.25)))
  }

  FALSE
}

pick_preferred_column <- function(names_vec, preferred) {
  lower <- tolower(names_vec)
  hits <- match(preferred, lower, nomatch = 0L)
  hits <- hits[hits > 0L]
  if (length(hits) == 0L) {
    return(NULL)
  }
  names_vec[[hits[[1L]]]]
}

detect_groupby <- function(obs) {
  preferred <- c(
    "sample_id", "section_id", "section", "sample", "library_id", "orig.ident",
    "imageid", "fov", "field_of_view", "slice"
  )
  direct_preferred <- pick_preferred_column(names(obs), preferred)
  if (!is.null(direct_preferred)) {
    return(direct_preferred)
  }

  candidates <- names(obs)[vapply(obs, is_groupby_candidate, logical(1))]
  if (length(candidates) > 0L) {
    return(pick_preferred_column(candidates, preferred) %||% candidates[[1L]])
  }

  categorical_cols <- names(obs)[vapply(obs, function(column) {
    is.factor(column) || is.character(column) || is.logical(column)
  }, logical(1))]
  if (length(categorical_cols) == 0L) {
    stop("Could not auto-detect a groupby column. Pass --groupby explicitly.")
  }

  pick_preferred_column(categorical_cols, preferred) %||% categorical_cols[[1L]]
}

detect_initial_color <- function(obs, groupby) {
  all_candidates <- setdiff(names(obs)[vapply(obs, is_color_candidate, logical(1))], groupby)
  lower_names <- tolower(all_candidates)
  is_coord_like <- lower_names %in% c(
    "x", "y", "z", "coord_x", "coord_y", "spatial_x", "spatial_y",
    "centroid_x", "centroid_y", "row", "col"
  )
  candidates <- all_candidates[!is_coord_like]
  if (length(candidates) == 0L) {
    candidates <- all_candidates
  }
  if (length(candidates) == 0L) {
    stop("Could not auto-detect an initial color column. Pass --initial-color explicitly.")
  }

  categorical_candidates <- candidates[vapply(obs[candidates], function(column) {
    is.factor(column) || is.character(column) || is.logical(column)
  }, logical(1))]

  preferred <- c(
    "cell_type", "celltype", "celltypes", "annotation", "annotations",
    "predicted_cell_type", "predicted.celltype", "cluster", "clusters",
    "leiden", "seurat_clusters", "subclass", "class"
  )

  pick_preferred_column(categorical_candidates, preferred) %||%
    pick_preferred_column(candidates, preferred) %||%
    if (length(categorical_candidates) > 0L) categorical_candidates[[1L]] else candidates[[1L]]
}

detect_additional_colors <- function(obs, groupby, initial_color) {
  candidates <- setdiff(names(obs)[vapply(obs, is_color_candidate, logical(1))], c(groupby, initial_color))
  if (length(candidates) == 0L) {
    return(NULL)
  }

  categorical <- candidates[vapply(obs[candidates], function(column) {
    is.factor(column) || is.character(column) || is.logical(column)
  }, logical(1))]
  numeric <- setdiff(candidates, categorical)

  utils::head(c(categorical, numeric), 3L)
}

default_output_path <- function(input_path) {
  input_abs <- normalizePath(input_path, mustWork = TRUE)
  outdir <- dirname(input_abs)
  stem <- tools::file_path_sans_ext(basename(input_abs))
  file.path(outdir, paste0(stem, "_karospace_buildr.html"))
}

input_path <- normalizePath(options$input, mustWork = TRUE)
obj <- read_karospace_source(input_path)
obs <- extract_obs(obj)

groupby <- options$groupby %||% detect_groupby(obs)
initial_color <- options[["initial-color"]] %||% detect_initial_color(obs, groupby)
additional_colors <- split_csv(options[["additional-colors"]]) %||%
  detect_additional_colors(obs, groupby, initial_color)
missing_additional_colors <- setdiff(additional_colors %||% character(), names(obs))
if (length(missing_additional_colors) > 0L) {
  warning(
    "Dropping missing additional colors: ",
    paste(missing_additional_colors, collapse = ", "),
    call. = FALSE
  )
  additional_colors <- setdiff(additional_colors, missing_additional_colors)
}
output_path <- options$output %||% default_output_path(input_path)
title <- options$title %||% tools::file_path_sans_ext(basename(input_path))
theme <- options$theme %||% "light"

cat("Input: ", input_path, "\n", sep = "")
cat("Detected groupby: ", groupby, "\n", sep = "")
cat("Detected initial color: ", initial_color, "\n", sep = "")
cat(
  "Additional colors: ",
  if (length(additional_colors %||% character()) > 0L) paste(additional_colors, collapse = ", ") else "<none>",
  "\n",
  sep = ""
)
cat("Output: ", normalizePath(output_path, mustWork = FALSE), "\n", sep = "")

categorical_cols <- names(obs)[vapply(obs, function(column) {
  is.factor(column) || is.character(column) || is.logical(column)
}, logical(1))]
numeric_cols <- names(obs)[vapply(obs, is.numeric, logical(1))]

cat(
  "Categorical columns: ",
  if (length(categorical_cols) > 0L) paste(categorical_cols, collapse = ", ") else "<none>",
  "\n",
  sep = ""
)
cat(
  "Numeric columns: ",
  if (length(numeric_cols) > 0L) paste(utils::head(numeric_cols, 12L), collapse = ", ") else "<none>",
  "\n",
  sep = ""
)

if (isTRUE(options$inspect)) {
  quit(save = "no", status = 0L)
}

export_karospace_viewer(
  input = obj,
  output_path = output_path,
  groupby = groupby,
  initial_color = initial_color,
  additional_colors = additional_colors,
  genes = split_csv(options$genes),
  metadata_columns = split_csv(options[["metadata-columns"]]),
  outline_by = options[["outline-by"]],
  title = title,
  theme = theme,
  min_panel_size = as.numeric(options[["min-panel-size"]] %||% 150),
  spot_size = as.numeric(options[["spot-size"]] %||% 2)
)

cat("Viewer written to ", normalizePath(output_path, mustWork = FALSE), "\n", sep = "")
