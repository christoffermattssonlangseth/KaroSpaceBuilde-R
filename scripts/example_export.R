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

if (isTRUE(options$help) || (is.null(options$input) && is.null(options$config))) {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/example_export.R --input path/to/object.rds [--output viewer.html]",
      "[--config build.json] [--write-config build.json]",
      "[--groupby sample_id] [--initial-color cell_type] [--additional-colors course,condition]",
      "[--metadata-input other_object.rds] [--metadata-input-columns col1,col2] [--metadata-prefix ext_]",
      "[--assay SCT] [--genes GENE1,GENE2] [--top-genes 200] [--lightweight] [--neighbor-mode spatial] [--neighbor-graph SCT_snn] [--neighbor-k 6] [--inspect] [--inspect-genes]",
      "[--marker-genes-groupby auto] [--marker-genes-top-n 20] [--marker-test mean_diff|wilcoxon]",
      "[--interaction-markers-groupby cell_type] [--interaction-markers-top-targets 8] [--interaction-markers-top-genes 12]",
      "[--interaction-markers-min-cells 30] [--interaction-markers-min-neighbors 1]",
      "[--neighbor-stats-permutations 0] [--neighbor-stats-seed 42]",
      "[--gene-query COL] [--gene-limit 50] [--title MyViewer] [--theme light]",
      sep = "\n"
    ),
    "\n"
  )
  quit(save = "no", status = if (isTRUE(options$help)) 0L else 1L)
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
  if (isTRUE(value)) {
    return(TRUE)
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

extract_graph_names <- function(x) {
  if (inherits(x, "Seurat")) {
    return(names(x@graphs))
  }
  character()
}

is_color_candidate <- function(column) {
  is.factor(column) || is.character(column) || is.logical(column) || is.numeric(column)
}

column_coverage <- function(column) {
  fraction_non_missing_values(column)
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
    column_coverage(obs[[column_name]])
  }, numeric(1))
  order_idx <- order(coverage, columns, decreasing = TRUE)
  columns <- columns[order_idx]
  coverage <- coverage[order_idx]

  limit <- min(length(columns), max_columns)
  preview <- sprintf("%s %s", columns[seq_len(limit)], format_percentage(coverage[seq_len(limit)]))
  paste(preview, collapse = ", ")
}

report_metadata_merge <- function(obj, obs) {
  report <- extract_metadata_merge_report(obj)
  if (is.null(report)) {
    return(invisible(NULL))
  }

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
  if (!is.data.frame(coverage_table) || nrow(coverage_table) == 0L) {
    return(invisible(report))
  }

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
  min_coverage <- 0.9
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

  coverage <- vapply(obs[candidates], column_coverage, numeric(1))
  high_coverage_candidates <- candidates[coverage[candidates] >= min_coverage]
  if (length(high_coverage_candidates) == 0L &&
      (is.factor(obs[[groupby]]) || is.character(obs[[groupby]]) || is.logical(obs[[groupby]]))) {
    return(groupby)
  }
  if (length(high_coverage_candidates) > 0L) {
    candidates <- high_coverage_candidates
  } else {
    candidates <- candidates[order(coverage[candidates], candidates, decreasing = TRUE)]
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
  min_coverage <- 0.9
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

  coverage <- vapply(obs[candidates], column_coverage, numeric(1))
  high_coverage_candidates <- candidates[coverage[candidates] >= min_coverage]
  if (length(high_coverage_candidates) > 0L) {
    candidates <- high_coverage_candidates
  } else {
    candidates <- candidates[order(coverage[candidates], candidates, decreasing = TRUE)]
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
  file.path(outdir, paste0(stem, "_karospacer.html"))
}

filter_gene_names <- function(gene_names, query = NULL) {
  gene_names <- as.character(gene_names %||% character())
  if (is.null(query) || !nzchar(query)) {
    return(gene_names)
  }

  keep <- grepl(tolower(query), tolower(gene_names), fixed = TRUE)
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

config_source_path <- if (!is.null(options$config)) {
  normalizePath(options$config, mustWork = TRUE)
} else {
  NULL
}
base_config <- if (!is.null(config_source_path)) {
  read_karospace_build_config(config_source_path)
} else {
  list(version = 1L)
}

input_path <- if (!is.null(options$input)) {
  normalizePath(options$input, mustWork = TRUE)
} else {
  base_config$input
}
if (is.null(input_path)) {
  stop("Could not resolve an input path. Pass --input or --config.")
}

metadata_input_path <- if (!is.null(options[["metadata-input"]])) {
  normalizePath(options[["metadata-input"]], mustWork = TRUE)
} else {
  base_config$metadata_input
}
metadata_input_columns <- if (!is.null(options[["metadata-input-columns"]])) {
  split_csv(options[["metadata-input-columns"]])
} else {
  base_config$metadata_input_columns
}
metadata_prefix <- options[["metadata-prefix"]] %||% base_config$metadata_prefix
obj <- prepare_karospace_input(
  input = input_path,
  metadata_input = metadata_input_path,
  metadata_input_columns = metadata_input_columns,
  metadata_prefix = metadata_prefix
)
obs <- extract_obs(obj)
available_assays <- extract_assay_names(obj)
available_graphs <- extract_graph_names(obj)
neighbor_mode <- options[["neighbor-mode"]] %||% base_config$neighbor_mode %||% "spatial"
neighbor_graph <- options[["neighbor-graph"]] %||% base_config$neighbor_graph
neighbor_k <- suppressWarnings(as.integer(
  if (!is.null(options[["neighbor-k"]])) {
    options[["neighbor-k"]]
  } else {
    base_config$neighbor_k %||% 6L
  }
))
if (is.na(neighbor_k) || neighbor_k < 1L) {
  neighbor_k <- 6L
}
merge_report <- report_metadata_merge(obj, obs)

groupby <- options$groupby %||% base_config$groupby %||% detect_groupby(obs)
initial_color <- options[["initial-color"]] %||% base_config$initial_color %||% detect_initial_color(obs, groupby)
assay_name <- options[["assay"]] %||% base_config$assay %||% detect_assay(obj)
genes <- if (!is.null(options$genes)) {
  split_csv(options$genes)
} else {
  base_config$genes
}
top_genes_n <- if (!is.null(options[["top-genes"]])) {
  suppressWarnings(as.integer(options[["top-genes"]]))
} else {
  base_config$top_genes_n
}
if (!is.null(options[["top-genes"]]) && (is.na(top_genes_n) || top_genes_n < 1L)) {
  stop("--top-genes must be a positive integer.")
}
if (!is.null(genes) && length(top_genes_n) > 0L) {
  warning("--top-genes is ignored because --genes was provided explicitly.", call. = FALSE)
  top_genes_n <- NULL
}
additional_colors <- split_csv(options[["additional-colors"]]) %||%
  base_config$additional_colors %||%
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
output_path <- if (!is.null(options$output)) {
  resolve_output_path(options$output)
} else {
  base_config$output %||% default_output_path(input_path)
}
title <- options$title %||% base_config$title %||% tools::file_path_sans_ext(basename(input_path))
theme <- options$theme %||% base_config$theme %||% "light"
lightweight <- if (!is.null(options[["lightweight"]])) {
  parse_bool_option(options[["lightweight"]], default = FALSE)
} else {
  isTRUE(base_config$lightweight %||% FALSE)
}
marker_genes_groupby <- split_csv(options[["marker-genes-groupby"]]) %||%
  base_config$marker_genes_groupby %||%
  if (isTRUE(lightweight)) "none" else "auto"
marker_genes_top_n <- suppressWarnings(as.integer(
  if (!is.null(options[["marker-genes-top-n"]])) {
    options[["marker-genes-top-n"]]
  } else {
    base_config$marker_genes_top_n %||% 20L
  }
))
if (is.na(marker_genes_top_n) || marker_genes_top_n < 1L) {
  marker_genes_top_n <- 20L
}
interaction_markers_groupby <- if (!is.null(options[["interaction-markers-groupby"]])) {
  split_csv(options[["interaction-markers-groupby"]])
} else {
  base_config$interaction_markers_groupby
}
if (isTRUE(lightweight) && is.null(options[["interaction-markers-groupby"]])) {
  interaction_markers_groupby <- "none"
}
interaction_markers_top_targets <- suppressWarnings(as.integer(
  if (!is.null(options[["interaction-markers-top-targets"]])) {
    options[["interaction-markers-top-targets"]]
  } else {
    base_config$interaction_markers_top_targets %||% 8L
  }
))
if (is.na(interaction_markers_top_targets) || interaction_markers_top_targets < 1L) {
  interaction_markers_top_targets <- 8L
}
interaction_markers_top_genes <- suppressWarnings(as.integer(
  if (!is.null(options[["interaction-markers-top-genes"]])) {
    options[["interaction-markers-top-genes"]]
  } else {
    base_config$interaction_markers_top_genes %||% 12L
  }
))
if (is.na(interaction_markers_top_genes) || interaction_markers_top_genes < 1L) {
  interaction_markers_top_genes <- 12L
}
interaction_markers_min_cells <- suppressWarnings(as.integer(
  if (!is.null(options[["interaction-markers-min-cells"]])) {
    options[["interaction-markers-min-cells"]]
  } else {
    base_config$interaction_markers_min_cells %||% 30L
  }
))
if (is.na(interaction_markers_min_cells) || interaction_markers_min_cells < 2L) {
  interaction_markers_min_cells <- 30L
}
interaction_markers_min_neighbors <- suppressWarnings(as.integer(
  if (!is.null(options[["interaction-markers-min-neighbors"]])) {
    options[["interaction-markers-min-neighbors"]]
  } else {
    base_config$interaction_markers_min_neighbors %||% 1L
  }
))
if (is.na(interaction_markers_min_neighbors) || interaction_markers_min_neighbors < 1L) {
  interaction_markers_min_neighbors <- 1L
}
marker_test <- tolower(trimws(options[["marker-test"]] %||% base_config$marker_test %||% "mean_diff"))
if (!marker_test %in% c("mean_diff", "wilcoxon")) {
  warning("Unknown --marker-test value '", marker_test, "'. Falling back to 'mean_diff'.", call. = FALSE)
  marker_test <- "mean_diff"
}
neighbor_stats_permutations <- suppressWarnings(as.integer(
  if (!is.null(options[["neighbor-stats-permutations"]])) {
    options[["neighbor-stats-permutations"]]
  } else {
    base_config$neighbor_stats_permutations %||% 0L
  }
))
if (is.na(neighbor_stats_permutations) || neighbor_stats_permutations < 0L) {
  neighbor_stats_permutations <- 0L
}
neighbor_stats_seed <- suppressWarnings(as.integer(
  if (!is.null(options[["neighbor-stats-seed"]])) {
    options[["neighbor-stats-seed"]]
  } else {
    base_config$neighbor_stats_seed %||% 42L
  }
))
if (is.na(neighbor_stats_seed)) {
  neighbor_stats_seed <- 42L
}
metadata_columns <- if (!is.null(options[["metadata-columns"]])) {
  split_csv(options[["metadata-columns"]])
} else {
  base_config$metadata_columns
}
outline_by <- options[["outline-by"]] %||% base_config$outline_by
min_panel_size <- suppressWarnings(as.numeric(
  if (!is.null(options[["min-panel-size"]])) {
    options[["min-panel-size"]]
  } else {
    base_config$min_panel_size %||% 150
  }
))
if (is.na(min_panel_size) || min_panel_size < 1) {
  min_panel_size <- 150
}
spot_size <- suppressWarnings(as.numeric(
  if (!is.null(options[["spot-size"]])) {
    options[["spot-size"]]
  } else {
    base_config$spot_size %||% 2
  }
))
if (is.na(spot_size) || spot_size < 0) {
  spot_size <- 2
}

resolved_config <- normalize_karospace_build_config(
  config = list(
    version = 1L,
    input = input_path,
    output = output_path,
    groupby = groupby,
    initial_color = initial_color,
    additional_colors = additional_colors,
    genes = genes,
    top_genes_n = top_genes_n,
    assay = assay_name,
    metadata_input = metadata_input_path,
    metadata_input_columns = metadata_input_columns,
    metadata_prefix = metadata_prefix,
    neighbor_mode = neighbor_mode,
    neighbor_graph = neighbor_graph,
    neighbor_k = neighbor_k,
    metadata_columns = metadata_columns,
    outline_by = outline_by,
    lightweight = lightweight,
    marker_genes_groupby = marker_genes_groupby,
    marker_genes_top_n = marker_genes_top_n,
    interaction_markers_groupby = interaction_markers_groupby,
    interaction_markers_top_targets = interaction_markers_top_targets,
    interaction_markers_top_genes = interaction_markers_top_genes,
    interaction_markers_min_cells = interaction_markers_min_cells,
    interaction_markers_min_neighbors = interaction_markers_min_neighbors,
    marker_test = marker_test,
    neighbor_stats_permutations = neighbor_stats_permutations,
    neighbor_stats_seed = neighbor_stats_seed,
    title = title,
    theme = theme,
    min_panel_size = min_panel_size,
    spot_size = spot_size
  ),
  require_required = TRUE,
  resolve_paths = FALSE
)
phase4_color_columns <- unique(c(initial_color, additional_colors %||% character()))
phase4_color_columns <- intersect(phase4_color_columns, names(obs))
phase4_color_data <- lapply(phase4_color_columns, function(column_name) build_color_column(obs[[column_name]]))
names(phase4_color_data) <- phase4_color_columns
resolved_marker_groupby <- resolve_phase4_groupby_columns(
  requested_groupby = unique(c(marker_genes_groupby, interaction_markers_groupby)),
  color_data = phase4_color_data
)
resolved_interaction_groupby <- resolve_phase4_groupby_columns(
  requested_groupby = interaction_markers_groupby,
  color_data = phase4_color_data
)

cat("Input: ", input_path, "\n", sep = "")
if (!is.null(metadata_input_path) && nzchar(metadata_input_path)) {
  cat(
    "Metadata input: ",
    normalizePath(metadata_input_path, mustWork = TRUE),
    "\n",
    sep = ""
  )
}
if (!is.null(config_source_path)) {
  cat("Config: ", config_source_path, "\n", sep = "")
}
cat("Detected groupby: ", groupby, "\n", sep = "")
cat("Detected initial color: ", initial_color, "\n", sep = "")
cat("Assay: ", assay_name %||% "<none>", "\n", sep = "")
cat("Top genes: ", top_genes_n %||% "<default>", "\n", sep = "")
cat("Lightweight: ", if (isTRUE(lightweight)) "true" else "false", "\n", sep = "")
cat("Neighbor mode: ", neighbor_mode, "\n", sep = "")
cat("Neighbor graph: ", neighbor_graph %||% "<auto>", "\n", sep = "")
cat("Neighbor k: ", neighbor_k, "\n", sep = "")
cat(
  "Marker colors: ",
  if (length(resolved_marker_groupby) > 0L) paste(resolved_marker_groupby, collapse = ", ") else "<none>",
  "\n",
  sep = ""
)
cat("Marker top N: ", marker_genes_top_n, "\n", sep = "")
cat(
  "Interaction colors: ",
  if (length(resolved_interaction_groupby) > 0L) paste(resolved_interaction_groupby, collapse = ", ") else "<none>",
  "\n",
  sep = ""
)
cat("Interaction top targets: ", interaction_markers_top_targets, "\n", sep = "")
cat("Interaction top genes: ", interaction_markers_top_genes, "\n", sep = "")
cat("Interaction min cells: ", interaction_markers_min_cells, "\n", sep = "")
cat("Interaction min neighbors: ", interaction_markers_min_neighbors, "\n", sep = "")
cat("Marker test: ", marker_test, "\n", sep = "")
cat("Neighbor stats permutations: ", neighbor_stats_permutations, "\n", sep = "")
cat("Neighbor stats seed: ", neighbor_stats_seed, "\n", sep = "")
cat(
  "Additional colors: ",
  if (length(additional_colors %||% character()) > 0L) paste(additional_colors, collapse = ", ") else "<none>",
  "\n",
  sep = ""
)
cat("Output: ", normalizePath(output_path, mustWork = FALSE), "\n", sep = "")

if (!is.null(options[["write-config"]])) {
  write_karospace_build_config(
    config = resolved_config,
    path = options[["write-config"]]
  )
  cat(
    "Config written to ",
    normalizePath(options[["write-config"]], mustWork = FALSE),
    "\n",
    sep = ""
  )
}

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
cat(
  "Available graphs: ",
  if (length(available_graphs) > 0L) paste(available_graphs, collapse = ", ") else "<none>",
  "\n",
  sep = ""
)

if (!is.null(merge_report) && is.null(options[["initial-color"]])) {
  low_coverage_merged <- merge_report$column_coverage$column[
    merge_report$column_coverage$coverage < 0.9
  ]
  low_coverage_merged <- intersect(low_coverage_merged, setdiff(names(obs)[vapply(obs, is_color_candidate, logical(1))], groupby))
  if (length(low_coverage_merged) > 0L) {
    cat(
      "Low-coverage merged color columns skipped from defaults: ",
      format_column_coverage(obs, low_coverage_merged),
      "\n",
      sep = ""
    )
  }
}

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
  title = resolved_config$title,
  theme = resolved_config$theme,
  min_panel_size = resolved_config$min_panel_size,
  spot_size = resolved_config$spot_size
)

cat("Viewer written to ", normalizePath(resolved_config$output, mustWork = FALSE), "\n", sep = "")
