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
      "[--metadata-input other_object.rds] [--metadata-input-columns col1,col2] [--metadata-prefix ext_]",
      "[--assay SCT] [--genes GENE1,GENE2] [--inspect] [--inspect-genes]",
      "[--gene-query COL] [--gene-limit 50] [--title MyViewer] [--theme light]",
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
    return(augment_seurat_obs(x))
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
    "Could not inspect the input. Supported inputs are list with obs, Seurat, SingleCellExperiment, and SpatialExperiment."
  )
}

extract_assay_names <- function(x) {
  if (inherits(x, "Seurat")) {
    return(names(x@assays))
  }

  if (inherits(x, "SpatialExperiment") || inherits(x, "SingleCellExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      return(character())
    }
    return(SummarizedExperiment::assayNames(x))
  }

  if (is.list(x) && is.list(x$assays)) {
    return(names(x$assays))
  }

  character()
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
  hits <- match(tolower(preferred), lower, nomatch = 0L)
  hits <- hits[hits > 0L]
  if (length(hits) == 0L) {
    return(NULL)
  }
  names_vec[[hits[[1L]]]]
}

detect_groupby <- function(obs) {
  preferred <- c(
    "sample_name", "sample_id", "sampleid", "section_id", "section", "sample",
    "library_id", "orig.ident",
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

  candidates <- candidates[vapply(obs[candidates], function(column) {
    length(unique_non_missing(column)) > 1L
  }, logical(1))]
  if (length(candidates) == 0L) {
    return(groupby)
  }

  categorical_candidates <- candidates[vapply(obs[candidates], function(column) {
    is.factor(column) || is.character(column) || is.logical(column)
  }, logical(1))]

  preferred <- c(
    "cell_type", "celltype", "celltypes", "annotation", "annotations",
    "predicted_cell_type", "predicted.celltype", "cluster", "clusters",
    "leiden", "seurat_clusters", "subclass", "class"
  )

  if (length(categorical_candidates) == 0L &&
      (is.factor(obs[[groupby]]) || is.character(obs[[groupby]]) || is.logical(obs[[groupby]]))) {
    return(groupby)
  }

  pick_preferred_column(categorical_candidates, preferred) %||%
    pick_preferred_column(candidates, preferred) %||%
    if (length(categorical_candidates) > 0L) categorical_candidates[[1L]] else candidates[[1L]]
}

detect_additional_colors <- function(obs, groupby, initial_color) {
  candidates <- setdiff(names(obs)[vapply(obs, is_color_candidate, logical(1))], c(groupby, initial_color))
  if (length(candidates) == 0L) {
    return(NULL)
  }

  candidates <- candidates[vapply(obs[candidates], function(column) {
    length(unique_non_missing(column)) > 1L
  }, logical(1))]
  if (length(candidates) == 0L) {
    return(NULL)
  }

  categorical <- candidates[vapply(obs[candidates], function(column) {
    is.factor(column) || is.character(column) || is.logical(column)
  }, logical(1))]
  numeric <- setdiff(candidates, categorical)

  utils::head(c(categorical, numeric), 3L)
}

detect_assay <- function(x) {
  assay_names <- extract_assay_names(x)
  if (length(assay_names) == 0L) {
    return(NULL)
  }

  preferred <- c("SCT", "logcounts", "Spatial", "RNA", "integrated", "counts")
  pick_preferred_column(assay_names, preferred) %||% assay_names[[1L]]
}

default_output_path <- function(input_path) {
  input_abs <- normalizePath(input_path, mustWork = TRUE)
  outdir <- dirname(input_abs)
  stem <- tools::file_path_sans_ext(basename(input_abs))
  file.path(outdir, paste0(stem, "_karospace_buildr.html"))
}

filter_gene_names <- function(gene_names, query = NULL) {
  gene_names <- as.character(gene_names %||% character())
  if (is.null(query) || !nzchar(query)) {
    return(gene_names)
  }

  keep <- grepl(query, gene_names, ignore.case = TRUE, fixed = TRUE)
  gene_names[keep]
}

format_gene_preview <- function(gene_names, limit = 50L) {
  gene_names <- as.character(gene_names %||% character())
  if (length(gene_names) == 0L) {
    return("<none>")
  }

  limit <- max(1L, as.integer(limit))
  preview <- utils::head(gene_names, limit)
  paste(preview, collapse = ", ")
}

input_path <- normalizePath(options$input, mustWork = TRUE)
metadata_input_path <- options[["metadata-input"]]
obj <- prepare_karospace_input(
  input = input_path,
  metadata_input = metadata_input_path,
  metadata_input_columns = split_csv(options[["metadata-input-columns"]]),
  metadata_prefix = options[["metadata-prefix"]]
)
obs <- extract_obs(obj)
available_assays <- extract_assay_names(obj)

groupby <- options$groupby %||% detect_groupby(obs)
initial_color <- options[["initial-color"]] %||% detect_initial_color(obs, groupby)
assay_name <- options[["assay"]] %||% detect_assay(obj)
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
if (!is.null(metadata_input_path) && nzchar(metadata_input_path)) {
  cat(
    "Metadata input: ",
    normalizePath(metadata_input_path, mustWork = TRUE),
    "\n",
    sep = ""
  )
}
cat("Detected groupby: ", groupby, "\n", sep = "")
cat("Detected initial color: ", initial_color, "\n", sep = "")
cat("Assay: ", assay_name %||% "<none>", "\n", sep = "")
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
cat(
  "Available assays: ",
  if (length(available_assays) > 0L) paste(available_assays, collapse = ", ") else "<none>",
  "\n",
  sep = ""
)

if (isTRUE(options$inspect) || isTRUE(options[["inspect-genes"]])) {
  gene_limit <- suppressWarnings(as.integer(options[["gene-limit"]] %||% 50L))
  if (is.na(gene_limit) || gene_limit < 1L) {
    gene_limit <- 50L
  }

  gene_info <- resolve_input_gene_names(
    x = obj,
    requested_assay = assay_name
  )
  filtered_genes <- filter_gene_names(
    gene_names = gene_info$gene_names,
    query = options[["gene-query"]]
  )

  cat("Gene assay: ", gene_info$assay_name %||% "<none>", "\n", sep = "")
  cat("Gene count: ", length(gene_info$gene_names), "\n", sep = "")
  if (!is.null(options[["gene-query"]]) && nzchar(options[["gene-query"]])) {
    cat("Gene query: ", options[["gene-query"]], "\n", sep = "")
    cat("Matched genes: ", length(filtered_genes), "\n", sep = "")
  }
  cat("Gene preview: ", format_gene_preview(filtered_genes, gene_limit), "\n", sep = "")
}

if (isTRUE(options$inspect) || isTRUE(options[["inspect-genes"]])) {
  quit(save = "no", status = 0L)
}

export_karospace_viewer(
  input = obj,
  output_path = output_path,
  groupby = groupby,
  initial_color = initial_color,
  additional_colors = additional_colors,
  genes = split_csv(options$genes),
  assay = assay_name,
  metadata_input = NULL,
  metadata_input_columns = NULL,
  metadata_prefix = NULL,
  metadata_columns = split_csv(options[["metadata-columns"]]),
  outline_by = options[["outline-by"]],
  title = title,
  theme = theme,
  min_panel_size = as.numeric(options[["min-panel-size"]] %||% 150),
  spot_size = as.numeric(options[["spot-size"]] %||% 2)
)

cat("Viewer written to ", normalizePath(output_path, mustWork = FALSE), "\n", sep = "")
