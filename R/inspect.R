extract_inspect_obs <- function(x) {
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

extract_inspect_assay_names <- function(x) {
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

extract_inspect_graph_names <- function(x) {
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

is_groupby_candidate <- function(column) {
  values <- unique_non_missing(column)
  n_unique <- length(values)
  if (is_categorical(column)) {
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

  categorical_cols <- names(obs)[vapply(obs, is_categorical, logical(1))]
  if (length(categorical_cols) == 0L) {
    stop("Could not auto-detect a groupby column. Pass groupby explicitly.")
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
    stop("Could not auto-detect an initial color column. Pass initial_color explicitly.")
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
      is_categorical(obs[[groupby]])) {
    return(groupby)
  }
  if (length(high_coverage_candidates) > 0L) {
    candidates <- high_coverage_candidates
  } else {
    candidates <- candidates[order(coverage[candidates], candidates, decreasing = TRUE)]
  }

  categorical_candidates <- candidates[vapply(obs[candidates], is_categorical, logical(1))]

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

  categorical <- candidates[vapply(obs[candidates], is_categorical, logical(1))]
  numeric <- setdiff(candidates, categorical)

  utils::head(c(categorical, numeric), 3L)
}

detect_assay <- function(x) {
  assay_names <- extract_inspect_assay_names(x)
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

build_column_summary <- function(obs, columns) {
  columns <- as.character(columns %||% character())
  lapply(columns, function(column_name) {
    column <- obs[[column_name]]
    non_missing <- unique_non_missing(column)
    list(
      name = column_name,
      type = if (is_categorical(column)) {
        "categorical"
      } else if (is.numeric(column)) {
        "numeric"
      } else {
        "other"
      },
      coverage = column_coverage(column),
      n_unique = length(non_missing),
      preview = as_json_array(utils::head(non_missing, 10L))
    )
  })
}

serialize_metadata_merge_report <- function(report) {
  if (is.null(report)) {
    return(NULL)
  }

  coverage_entries <- list()
  coverage_table <- report$column_coverage
  if (is.data.frame(coverage_table) && nrow(coverage_table) > 0L) {
    coverage_entries <- lapply(seq_len(nrow(coverage_table)), function(i) {
      list(
        column = as.character(coverage_table$column[[i]]),
        coverage = as.numeric(coverage_table$coverage[[i]])
      )
    })
  }

  list(
    primary_rows = report$primary_rows,
    secondary_rows = report$secondary_rows,
    overlap_rows = report$overlap_rows,
    overlap_fraction = report$overlap_fraction,
    added_columns = as_json_array(report$added_columns),
    skipped_columns = as_json_array(report$skipped_columns),
    column_coverage = coverage_entries
  )
}

inspect_karospace_input <- function(
  input,
  metadata_input = NULL,
  metadata_input_columns = NULL,
  metadata_prefix = NULL,
  assay = NULL,
  gene_query = NULL,
  gene_limit = 50L
) {
  input_path <- if (is.character(input) && length(input) == 1L && file.exists(input)) {
    normalizePath(input, mustWork = TRUE)
  } else {
    NULL
  }
  metadata_input_path <- if (is.character(metadata_input) && length(metadata_input) == 1L && file.exists(metadata_input)) {
    normalizePath(metadata_input, mustWork = TRUE)
  } else {
    NULL
  }

  obj <- prepare_karospace_input(
    input = input,
    metadata_input = metadata_input,
    metadata_input_columns = metadata_input_columns,
    metadata_prefix = metadata_prefix
  )
  obs <- extract_inspect_obs(obj)
  categorical_cols <- names(obs)[vapply(obs, is_categorical, logical(1))]
  numeric_cols <- names(obs)[vapply(obs, is.numeric, logical(1))]
  available_assays <- extract_inspect_assay_names(obj)
  available_graphs <- extract_inspect_graph_names(obj)
  merge_report <- extract_metadata_merge_report(obj)

  groupby <- detect_groupby(obs)
  initial_color <- detect_initial_color(obs, groupby)
  additional_colors <- detect_additional_colors(obs, groupby, initial_color)
  assay_name <- assay %||% detect_assay(obj)
  gene_limit <- suppressWarnings(as.integer(gene_limit %||% 50L))
  if (is.na(gene_limit) || gene_limit < 1L) {
    gene_limit <- 50L
  }

  gene_info <- resolve_input_gene_names(
    x = obj,
    requested_assay = assay_name
  )
  filtered_genes <- filter_gene_names(
    gene_names = gene_info$gene_names,
    query = gene_query
  )

  low_coverage_merged <- character()
  if (!is.null(merge_report)) {
    low_coverage_merged <- merge_report$column_coverage$column[
      merge_report$column_coverage$coverage < 0.9
    ]
    low_coverage_merged <- intersect(
      low_coverage_merged,
      setdiff(names(obs)[vapply(obs, is_color_candidate, logical(1))], groupby)
    )
  }

  default_config <- normalize_karospace_build_config(
    config = list(
      version = 1L,
      input = input_path %||% input,
      output = if (!is.null(input_path)) default_output_path(input_path) else NULL,
      groupby = groupby,
      initial_color = initial_color,
      additional_colors = additional_colors,
      assay = assay_name,
      metadata_input = metadata_input_path %||% metadata_input,
      metadata_input_columns = metadata_input_columns,
      metadata_prefix = metadata_prefix,
      neighbor_mode = "spatial",
      marker_genes_groupby = "auto",
      title = if (!is.null(input_path)) tools::file_path_sans_ext(basename(input_path)) else "KaroSpace",
      theme = "light"
    ),
    require_required = FALSE,
    resolve_paths = FALSE
  )

  list(
    version = 1L,
    input = input_path,
    metadata_input = metadata_input_path,
    source_class = as_json_array(class(obj)),
    default_config = karospace_build_config_for_json(default_config),
    defaults = list(
      groupby = groupby,
      initial_color = initial_color,
      additional_colors = as_json_array(additional_colors),
      assay = assay_name,
      output = default_config$output,
      neighbor_mode = "spatial"
    ),
    columns = list(
      categorical = as_json_array(categorical_cols),
      numeric = as_json_array(numeric_cols),
      summary = build_column_summary(obs, names(obs))
    ),
    assays = list(
      available = as_json_array(available_assays),
      selected = assay_name
    ),
    graphs = list(
      available = as_json_array(available_graphs)
    ),
    genes = list(
      assay = gene_info$assay_name,
      count = length(gene_info$gene_names),
      query = gene_query %||% NULL,
      matched_count = length(filtered_genes),
      preview = as_json_array(utils::head(filtered_genes, gene_limit))
    ),
    metadata_merge = serialize_metadata_merge_report(merge_report),
    warnings = list(
      low_coverage_merged_color_columns = as_json_array(low_coverage_merged)
    )
  )
}
