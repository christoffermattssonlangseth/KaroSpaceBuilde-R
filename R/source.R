read_karospace_source <- function(input) {
  if (is.character(input) && length(input) == 1 && file.exists(input)) {
    rds_result <- tryCatch(
      readRDS(input),
      error = function(err) err
    )
    if (!inherits(rds_result, "error")) {
      return(rds_result)
    }

    load_env <- new.env(parent = emptyenv())
    load_result <- tryCatch(
      load(input, envir = load_env),
      error = function(err) err
    )
    if (!inherits(load_result, "error")) {
      if (length(load_result) != 1) {
        stop(
          "Input file is an R workspace/archive, not a single-object RDS, and contains ",
          length(load_result),
          " objects: ",
          paste(load_result, collapse = ", "),
          "."
        )
      }
      return(get(load_result[[1]], envir = load_env, inherits = FALSE))
    }

    stop(
      "Could not read input as either an RDS file or an R workspace. ",
      "readRDS error: ", conditionMessage(rds_result), ". ",
      "load error: ", conditionMessage(load_result), "."
    )
  }
  input
}

prepare_karospace_input <- function(
  input,
  metadata_input = NULL,
  metadata_input_columns = NULL,
  metadata_prefix = NULL
) {
  obj <- read_karospace_source(input)
  if (is.null(metadata_input)) {
    return(obj)
  }

  metadata_obj <- read_karospace_source(metadata_input)
  merge_external_metadata(
    primary = obj,
    secondary = metadata_obj,
    columns = metadata_input_columns,
    prefix = metadata_prefix
  )
}

count_non_missing_values <- function(x) {
  if (is_categorical(x)) {
    values <- as.character(x)
    return(sum(!is.na(values) & nzchar(values)))
  }

  if (is.numeric(x)) {
    return(sum(is.finite(x)))
  }

  sum(!is.na(x))
}

fraction_non_missing_values <- function(x) {
  total <- length(x)
  if (total == 0L) {
    return(0)
  }
  count_non_missing_values(x) / total
}

extract_metadata_table <- function(x) {
  if (inherits(x, "Seurat")) {
    obs <- augment_seurat_obs(x)
    ids <- tryCatch(SeuratObject::Idents(x), error = function(err) NULL)
    if (!is.null(ids) && length(unique_non_missing(ids)) > 1L && !("active_ident" %in% names(obs))) {
      obs$active_ident <- as.character(ids)
    }
    return(obs)
  }

  if (inherits(x, "SpatialExperiment") || inherits(x, "SingleCellExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment is required to extract metadata from this object.")
    }
    return(as.data.frame(SummarizedExperiment::colData(x)))
  }

  if (is.list(x) && is.data.frame(x$obs)) {
    return(x$obs)
  }

  if (is.data.frame(x)) {
    return(x)
  }

  stop(
    "Could not extract metadata from input class: ",
    paste(class(x), collapse = ", ")
  )
}

assign_metadata_table <- function(x, obs) {
  if (inherits(x, "Seurat")) {
    x@meta.data <- obs
    return(x)
  }

  if (inherits(x, "SpatialExperiment") || inherits(x, "SingleCellExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
        !requireNamespace("S4Vectors", quietly = TRUE)) {
      stop("SummarizedExperiment and S4Vectors are required to assign metadata to this object.")
    }
    SummarizedExperiment::colData(x) <- S4Vectors::DataFrame(obs)
    return(x)
  }

  if (is.list(x)) {
    x$obs <- obs
    return(x)
  }

  if (is.data.frame(x)) {
    return(obs)
  }

  stop(
    "Could not assign metadata back to input class: ",
    paste(class(x), collapse = ", ")
  )
}

attach_metadata_merge_report <- function(x, report) {
  if (inherits(x, "Seurat") && "misc" %in% slotNames(x)) {
    x@misc$karospace_metadata_merge_report <- report
    return(x)
  }

  if ((inherits(x, "SpatialExperiment") || inherits(x, "SingleCellExperiment")) &&
      requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    meta <- SummarizedExperiment::metadata(x)
    meta$karospace_metadata_merge_report <- report
    SummarizedExperiment::metadata(x) <- meta
    return(x)
  }

  attr(x, "karospace_metadata_merge_report") <- report
  x
}

extract_metadata_merge_report <- function(x) {
  if (inherits(x, "Seurat") && "misc" %in% slotNames(x)) {
    return(x@misc$karospace_metadata_merge_report %||% NULL)
  }

  if ((inherits(x, "SpatialExperiment") || inherits(x, "SingleCellExperiment")) &&
      requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    return(SummarizedExperiment::metadata(x)$karospace_metadata_merge_report %||% NULL)
  }

  attr(x, "karospace_metadata_merge_report") %||% NULL
}

merge_external_metadata <- function(primary, secondary, columns = NULL, prefix = NULL) {
  primary_obs <- extract_metadata_table(primary)
  secondary_obs <- extract_metadata_table(secondary)

  primary_ids <- rownames(primary_obs)
  secondary_ids <- rownames(secondary_obs)
  if (is.null(primary_ids) || length(primary_ids) == 0) {
    stop("Primary input metadata must have row names for metadata merging.")
  }
  if (is.null(secondary_ids) || length(secondary_ids) == 0) {
    stop("Secondary metadata input must have row names for metadata merging.")
  }

  overlap_ids <- intersect(primary_ids, secondary_ids)
  if (length(overlap_ids) == 0) {
    stop("No overlapping row names found between primary input and metadata input.")
  }

  merge_columns <- colnames(secondary_obs)
  if (!is.null(columns) && length(columns) > 0) {
    columns <- as.character(columns)
    missing_columns <- setdiff(columns, merge_columns)
    if (length(missing_columns) > 0) {
      warning(
        "Dropping missing metadata_input columns: ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
    merge_columns <- intersect(columns, merge_columns)
  }

  if (length(merge_columns) == 0) {
    stop("No metadata columns available to merge from metadata input.")
  }

  idx <- match(primary_ids, secondary_ids)
  aligned <- secondary_obs[idx, merge_columns, drop = FALSE]
  rownames(aligned) <- primary_ids

  merged <- primary_obs
  skipped_columns <- character()
  added_columns <- character()
  for (column_name in merge_columns) {
    target_name <- column_name
    if (target_name %in% names(merged)) {
      if (is.null(prefix) || !nzchar(prefix)) {
        skipped_columns <- c(skipped_columns, column_name)
        next
      }
      target_name <- make.unique(c(names(merged), paste0(prefix, column_name)))[length(names(merged)) + 1L]
    }
    merged[[target_name]] <- aligned[[column_name]]
    added_columns <- c(added_columns, target_name)
  }

  if (length(skipped_columns) > 0) {
    warning(
      "Skipping metadata_input columns already present in the primary input: ",
      paste(skipped_columns, collapse = ", "),
      ". Pass metadata_prefix to keep both versions.",
      call. = FALSE
    )
  }

  report <- list(
    primary_rows = length(primary_ids),
    secondary_rows = length(secondary_ids),
    overlap_rows = length(overlap_ids),
    overlap_fraction = length(overlap_ids) / length(primary_ids),
    requested_columns = as.character(columns %||% character()),
    added_columns = added_columns,
    skipped_columns = skipped_columns,
    column_coverage = if (length(added_columns) == 0) {
      data.frame(
        column = character(),
        non_missing = integer(),
        total = integer(),
        coverage = numeric(),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        column = added_columns,
        non_missing = vapply(added_columns, function(column_name) {
          count_non_missing_values(merged[[column_name]])
        }, integer(1)),
        total = rep.int(length(primary_ids), length(added_columns)),
        coverage = vapply(added_columns, function(column_name) {
          fraction_non_missing_values(merged[[column_name]])
        }, numeric(1)),
        stringsAsFactors = FALSE
      )
    }
  )

  attach_metadata_merge_report(
    x = assign_metadata_table(primary, merged),
    report = report
  )
}

normalize_input_source <- function(
  x,
  groupby,
  initial_color,
  additional_colors = NULL,
  genes = NULL,
  top_genes_n = NULL,
  assay = NULL,
  neighbor_mode = "spatial",
  neighbor_graph = NULL,
  neighbor_k = 6L,
  metadata_columns = NULL,
  outline_by = NULL
) {
  if (inherits(x, "Seurat")) {
    return(normalize_seurat_object(
      x = x,
      groupby = groupby,
      initial_color = initial_color,
      additional_colors = additional_colors,
      genes = genes,
      top_genes_n = top_genes_n,
      assay = assay,
      neighbor_mode = neighbor_mode,
      neighbor_graph = neighbor_graph,
      neighbor_k = neighbor_k,
      metadata_columns = metadata_columns,
      outline_by = outline_by
    ))
  }

  if (inherits(x, "SpatialExperiment")) {
    return(normalize_spatial_experiment(
      x = x,
      groupby = groupby,
      initial_color = initial_color,
      additional_colors = additional_colors,
      genes = genes,
      top_genes_n = top_genes_n,
      assay = assay,
      neighbor_mode = neighbor_mode,
      neighbor_graph = neighbor_graph,
      neighbor_k = neighbor_k,
      metadata_columns = metadata_columns,
      outline_by = outline_by
    ))
  }

  if (inherits(x, "SingleCellExperiment")) {
    return(normalize_single_cell_experiment(
      x = x,
      groupby = groupby,
      initial_color = initial_color,
      additional_colors = additional_colors,
      genes = genes,
      top_genes_n = top_genes_n,
      assay = assay,
      neighbor_mode = neighbor_mode,
      neighbor_graph = neighbor_graph,
      neighbor_k = neighbor_k,
      metadata_columns = metadata_columns,
      outline_by = outline_by
    ))
  }

  if (is.list(x)) {
    return(normalize_list_source(
      x = x,
      groupby = groupby,
      initial_color = initial_color,
      additional_colors = additional_colors,
      genes = genes,
      top_genes_n = top_genes_n,
      assay = assay,
      neighbor_mode = neighbor_mode,
      neighbor_graph = neighbor_graph,
      neighbor_k = neighbor_k,
      metadata_columns = metadata_columns,
      outline_by = outline_by
    ))
  }

  stop(
    "Unsupported input class: ",
    paste(class(x), collapse = ", "),
    ". Supported inputs are list, Seurat, SingleCellExperiment, and SpatialExperiment."
  )
}

normalize_seurat_object <- function(
  x,
  groupby,
  initial_color,
  additional_colors = NULL,
  genes = NULL,
  top_genes_n = NULL,
  assay = NULL,
  neighbor_mode = "spatial",
  neighbor_graph = NULL,
  neighbor_k = 6L,
  metadata_columns = NULL,
  outline_by = NULL
) {
  obs <- augment_seurat_obs(x)
  if (!is.data.frame(obs)) {
    stop("Seurat object does not have a usable meta.data frame.")
  }

  cell_names <- rownames(obs)
  if (is.null(cell_names) || length(cell_names) == 0) {
    stop("Seurat meta.data must have row names matching cell names.")
  }

  coord_df <- resolve_seurat_coordinate_df(
    x = x,
    cell_names = cell_names
  )

  coords <- as.matrix(coord_df[, c("x", "y"), drop = FALSE])
  mode(coords) <- "numeric"

  umap <- NULL
  if ("umap" %in% names(x@reductions)) {
    embeddings <- x@reductions$umap@cell.embeddings
    embeddings <- embeddings[cell_names, , drop = FALSE]
    if (ncol(embeddings) >= 2) {
      umap <- as.matrix(embeddings[, seq_len(2), drop = FALSE])
      mode(umap) <- "numeric"
    }
  }

  expression_info <- resolve_seurat_expression(
    x = x,
    requested_assay = assay
  )

  neighbor_info <- resolve_seurat_neighbor_info(
    x = x,
    obs = obs,
    cell_names = cell_names,
    coords = coords,
    groupby = groupby,
    neighbor_mode = neighbor_mode,
    neighbor_graph = neighbor_graph,
    neighbor_k = neighbor_k
  )

  normalize_list_source(
    x = list(
      obs = obs,
      coordinates = coords,
      umap = umap,
      expression = expression_info$expression,
      gene_names = expression_info$gene_names,
      expression_assay = expression_info$assay_name,
      neighbor_edges_by_section = neighbor_info$section_edges,
      neighbors_key = neighbor_info$key
    ),
    groupby = groupby,
    initial_color = initial_color,
    additional_colors = additional_colors,
    genes = genes,
    top_genes_n = top_genes_n,
    assay = expression_info$assay_name,
    neighbor_mode = neighbor_mode,
    neighbor_graph = neighbor_graph,
    neighbor_k = neighbor_k,
    metadata_columns = metadata_columns,
    outline_by = outline_by
  )
}

augment_seurat_obs <- function(x) {
  obs <- x@meta.data
  if (!is.data.frame(obs)) {
    return(obs)
  }

  staffli_obj <- NULL
  if (length(x@tools) > 0 && "Staffli" %in% names(x@tools)) {
    staffli_obj <- x@tools[["Staffli"]]
  }
  staffli_obs <- extract_staffli_obs_df(staffli_obj)
  if (is.null(staffli_obs) || nrow(staffli_obs) == 0) {
    return(obs)
  }

  cell_names <- rownames(obs)
  if (is.null(cell_names) || length(cell_names) == 0) {
    return(obs)
  }

  keep <- intersect(cell_names, rownames(staffli_obs))
  if (length(keep) == 0) {
    return(obs)
  }

  aligned <- staffli_obs[cell_names, , drop = FALSE]
  merged <- obs
  for (column_name in colnames(aligned)) {
    if (!(column_name %in% names(merged))) {
      merged[[column_name]] <- aligned[[column_name]]
      next
    }

    current <- merged[[column_name]]
    fill <- is.na(current) | (!nzchar(as.character(current)))
    if (any(fill)) {
      current[fill] <- aligned[[column_name]][fill]
      merged[[column_name]] <- current
    }
  }

  merged
}

resolve_seurat_coordinate_df <- function(x, cell_names) {
  image_coords <- extract_seurat_image_coordinate_df(x)
  staffli_obj <- NULL
  if (length(x@tools) > 0 && "Staffli" %in% names(x@tools)) {
    staffli_obj <- x@tools[["Staffli"]]
  }
  staffli_coords <- extract_staffli_coordinate_df(staffli_obj)

  coord_df <- image_coords
  if (!is.null(staffli_coords)) {
    if (is.null(coord_df)) {
      coord_df <- staffli_coords
    } else {
      missing_from_images <- setdiff(rownames(staffli_coords), rownames(coord_df))
      if (length(missing_from_images) > 0) {
        coord_df <- rbind(coord_df, staffli_coords[missing_from_images, , drop = FALSE])
      }
    }
  }

  if (is.null(coord_df) || nrow(coord_df) == 0) {
    stop(
      "Seurat object has no spatial coordinates in @images or tools$Staffli."
    )
  }

  missing_cells <- setdiff(cell_names, rownames(coord_df))
  if (length(missing_cells) > 0) {
    checked_sources <- c(
      if (length(x@images) > 0) "@images" else NULL,
      if (!is.null(staffli_coords)) "tools$Staffli" else NULL
    )
    stop(
      "Seurat spatial coordinates do not cover all cells in meta.data. Missing ",
      length(missing_cells),
      " cells after checking ",
      paste(checked_sources, collapse = " and "),
      "."
    )
  }

  coord_df[cell_names, , drop = FALSE]
}

extract_seurat_image_coordinate_df <- function(x) {
  if (length(x@images) == 0) {
    return(NULL)
  }

  coord_frames <- list()
  for (image_name in names(x@images)) {
    image_obj <- x@images[[image_name]]
    if (!("coordinates" %in% slotNames(image_obj))) {
      next
    }

    image_coords <- image_obj@coordinates
    if (!is.data.frame(image_coords) || nrow(image_coords) == 0) {
      next
    }

    image_x_col <- pick_first_existing(colnames(image_coords), c("imagecol", "col", "x"))
    image_y_col <- pick_first_existing(colnames(image_coords), c("imagerow", "row", "y"))
    if (is.null(image_x_col) || is.null(image_y_col)) {
      numeric_cols <- colnames(image_coords)[vapply(image_coords, is.numeric, logical(1))]
      if (length(numeric_cols) < 2) {
        next
      }
      image_x_col <- image_x_col %||% numeric_cols[[1]]
      image_y_col <- image_y_col %||% numeric_cols[[2]]
    }

    coord_frames[[image_name]] <- data.frame(
      x = as.numeric(image_coords[[image_x_col]]),
      y = as.numeric(image_coords[[image_y_col]]),
      row.names = rownames(image_coords),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (length(coord_frames) == 0) {
    return(NULL)
  }

  coord_df <- do.call(rbind, unname(coord_frames))
  coord_df[!duplicated(rownames(coord_df)), , drop = FALSE]
}

extract_staffli_coordinate_df <- function(staffli) {
  if (is.null(staffli)) {
    return(NULL)
  }

  staffli_attrs <- attributes(staffli)
  meta_data <- staffli_attrs$meta_data %||% NULL
  if (!is.data.frame(meta_data) || nrow(meta_data) == 0) {
    return(NULL)
  }

  barcode_col <- pick_first_existing(
    colnames(meta_data),
    c("barcode", "cell", "cell_id", "spot", "spot_id")
  )
  x_col <- pick_first_existing(
    colnames(meta_data),
    c("pxl_col_in_fullres", "imagecol", "col", "x")
  )
  y_col <- pick_first_existing(
    colnames(meta_data),
    c("pxl_row_in_fullres", "imagerow", "row", "y")
  )

  if (is.null(barcode_col) || is.null(x_col) || is.null(y_col)) {
    return(NULL)
  }

  coord_df <- data.frame(
    barcode = as.character(meta_data[[barcode_col]]),
    x = as.numeric(meta_data[[x_col]]),
    y = as.numeric(meta_data[[y_col]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  coord_df <- coord_df[!is.na(coord_df$barcode) & nzchar(coord_df$barcode), , drop = FALSE]
  coord_df <- coord_df[!duplicated(coord_df$barcode), , drop = FALSE]
  rownames(coord_df) <- coord_df$barcode
  coord_df$barcode <- NULL
  coord_df
}

extract_staffli_obs_df <- function(staffli) {
  if (is.null(staffli)) {
    return(NULL)
  }

  staffli_attrs <- attributes(staffli)
  meta_data <- staffli_attrs$meta_data %||% NULL
  if (!is.data.frame(meta_data) || nrow(meta_data) == 0) {
    return(NULL)
  }

  barcode_col <- pick_first_existing(
    colnames(meta_data),
    c("barcode", "cell", "cell_id", "spot", "spot_id")
  )
  if (is.null(barcode_col)) {
    return(NULL)
  }

  out <- data.frame(
    row.names = as.character(meta_data[[barcode_col]]),
    stringsAsFactors = FALSE
  )

  sample_id_col <- pick_first_existing(colnames(meta_data), c("sampleID", "sample_id", "sample"))
  if (!is.null(sample_id_col)) {
    out$sampleID <- as.character(meta_data[[sample_id_col]])
  }

  sample_map <- extract_staffli_sample_map(staffli_attrs)
  if (!is.null(sample_map) && "sampleID" %in% names(out)) {
    matched <- match(out$sampleID, sample_map$sampleID)
    out$sample_name <- sample_map$sample_name[matched]
  }

  out
}

extract_staffli_sample_map <- function(staffli_attrs) {
  imgs <- as.character(staffli_attrs$imgs %||% character())
  if (length(imgs) == 0) {
    return(NULL)
  }

  sample_ids <- NULL
  image_info <- staffli_attrs$image_info %||% NULL
  if (is.data.frame(image_info) && "sampleID" %in% colnames(image_info)) {
    sample_ids <- as.character(image_info$sampleID)
  }
  if (is.null(sample_ids) || length(sample_ids) != length(imgs)) {
    sample_ids <- as.character(seq_along(imgs))
  }

  sample_names <- vapply(
    imgs,
    function(path) basename(dirname(dirname(path))),
    character(1)
  )

  data.frame(
    sampleID = sample_ids,
    sample_name = sample_names,
    stringsAsFactors = FALSE
  )
}

normalize_single_cell_experiment <- function(
  x,
  groupby,
  initial_color,
  additional_colors = NULL,
  genes = NULL,
  top_genes_n = NULL,
  assay = NULL,
  neighbor_mode = "spatial",
  neighbor_graph = NULL,
  neighbor_k = 6L,
  metadata_columns = NULL,
  outline_by = NULL
) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment support requires the SingleCellExperiment package.")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment support requires the SummarizedExperiment package.")
  }

  obs <- as.data.frame(SummarizedExperiment::colData(x))
  spatial <- NULL
  if ("spatial" %in% SingleCellExperiment::reducedDimNames(x)) {
    spatial <- SingleCellExperiment::reducedDim(x, "spatial")
  }
  if (is.null(spatial) && "Spatial" %in% SingleCellExperiment::reducedDimNames(x)) {
    spatial <- SingleCellExperiment::reducedDim(x, "Spatial")
  }
  if (is.null(spatial)) {
    stop("No spatial coordinates found. Expected reducedDim named 'spatial' or 'Spatial'.")
  }

  umap <- NULL
  if ("UMAP" %in% SingleCellExperiment::reducedDimNames(x)) {
    umap <- SingleCellExperiment::reducedDim(x, "UMAP")
  } else if ("X_umap" %in% SingleCellExperiment::reducedDimNames(x)) {
    umap <- SingleCellExperiment::reducedDim(x, "X_umap")
  }

  expression_info <- resolve_summarized_experiment_expression(
    x = x,
    requested_assay = assay
  )

  normalize_list_source(
    x = list(
      obs = obs,
      coordinates = spatial,
      umap = umap,
      expression = expression_info$expression,
      gene_names = expression_info$gene_names,
      expression_assay = expression_info$assay_name
    ),
    groupby = groupby,
    initial_color = initial_color,
    additional_colors = additional_colors,
    genes = genes,
    top_genes_n = top_genes_n,
    assay = expression_info$assay_name,
    neighbor_mode = neighbor_mode,
    neighbor_graph = neighbor_graph,
    neighbor_k = neighbor_k,
    metadata_columns = metadata_columns,
    outline_by = outline_by
  )
}

normalize_spatial_experiment <- function(
  x,
  groupby,
  initial_color,
  additional_colors = NULL,
  genes = NULL,
  top_genes_n = NULL,
  assay = NULL,
  neighbor_mode = "spatial",
  neighbor_graph = NULL,
  neighbor_k = 6L,
  metadata_columns = NULL,
  outline_by = NULL
) {
  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("SpatialExperiment support requires the SpatialExperiment package.")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SpatialExperiment support requires the SummarizedExperiment package.")
  }

  spatial <- SpatialExperiment::spatialCoords(x)
  obs <- as.data.frame(SummarizedExperiment::colData(x))

  umap <- NULL
  if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    rd_names <- SingleCellExperiment::reducedDimNames(x)
    if ("UMAP" %in% rd_names) {
      umap <- SingleCellExperiment::reducedDim(x, "UMAP")
    } else if ("X_umap" %in% rd_names) {
      umap <- SingleCellExperiment::reducedDim(x, "X_umap")
    }
  }

  expression_info <- resolve_summarized_experiment_expression(
    x = x,
    requested_assay = assay
  )

  normalize_list_source(
    x = list(
      obs = obs,
      coordinates = spatial,
      umap = umap,
      expression = expression_info$expression,
      gene_names = expression_info$gene_names,
      expression_assay = expression_info$assay_name
    ),
    groupby = groupby,
    initial_color = initial_color,
    additional_colors = additional_colors,
    genes = genes,
    top_genes_n = top_genes_n,
    assay = expression_info$assay_name,
    neighbor_mode = neighbor_mode,
    neighbor_graph = neighbor_graph,
    neighbor_k = neighbor_k,
    metadata_columns = metadata_columns,
    outline_by = outline_by
  )
}

normalize_list_source <- function(
  x,
  groupby,
  initial_color,
  additional_colors = NULL,
  genes = NULL,
  top_genes_n = NULL,
  assay = NULL,
  neighbor_mode = "spatial",
  neighbor_graph = NULL,
  neighbor_k = 6L,
  metadata_columns = NULL,
  outline_by = NULL
) {
  obs <- x$obs
  if (!is.data.frame(obs)) {
    stop("List input must include an 'obs' data.frame.")
  }

  coords <- as_plain_matrix(x$coordinates %||% x$coords)
  if (is.null(coords) || ncol(coords) < 2) {
    stop("List input must include 'coordinates' with at least two columns.")
  }
  coords <- coords[, seq_len(2), drop = FALSE]
  mode(coords) <- "numeric"

  if (nrow(obs) != nrow(coords)) {
    stop("The number of rows in obs must match the number of coordinate rows.")
  }

  if (!groupby %in% names(obs)) {
    stop("groupby column not found in obs: ", groupby)
  }
  if (!initial_color %in% names(obs)) {
    stop("initial_color column not found in obs: ", initial_color)
  }

  additional_colors <- unique(as.character(additional_colors %||% character()))
  missing_colors <- setdiff(additional_colors, names(obs))
  if (length(missing_colors) > 0) {
    warning(
      "Dropping missing additional_colors: ",
      paste(missing_colors, collapse = ", "),
      call. = FALSE
    )
    additional_colors <- setdiff(additional_colors, missing_colors)
  }

  umap <- as_plain_matrix(x$umap)
  if (!is.null(umap)) {
    if (nrow(umap) != nrow(obs) || ncol(umap) < 2) {
      stop("umap must have the same number of rows as obs and at least two columns.")
    }
    umap <- umap[, seq_len(2), drop = FALSE]
    mode(umap) <- "numeric"
  }

  expression_info <- normalize_expression(
    expression = x$expression %||% x$expr,
    gene_names = x$gene_names %||% x$genes_available %||% rownames(x$expression %||% x$expr),
    n_cells = nrow(obs),
    cell_names = rownames(obs)
  )

  neighbor_info <- resolve_list_neighbor_info(
    x = x,
    obs = obs,
    coords = coords,
    groupby = groupby,
    neighbor_mode = neighbor_mode,
    neighbor_graph = neighbor_graph,
    neighbor_k = neighbor_k
  )

  list(
    obs = obs,
    coords = coords,
    umap = umap,
    expression = expression_info$expression,
    gene_names = expression_info$gene_names,
    selected_genes = resolve_selected_genes(
      genes = genes,
      gene_names = expression_info$gene_names,
      expression = expression_info$expression,
      top_genes_n = top_genes_n
    ),
    expression_assay = x$expression_assay %||% assay,
    groupby = groupby,
    initial_color = initial_color,
    additional_colors = unique(c(initial_color, additional_colors)),
    metadata_columns = resolve_metadata_columns(
      obs = obs,
      groupby = groupby,
      metadata_columns = metadata_columns
    ),
    section_edges = neighbor_info$section_edges,
    neighbors_key = neighbor_info$key,
    outline_by = outline_by
  )
}

resolve_seurat_neighbor_info <- function(
  x,
  obs,
  cell_names,
  coords,
  groupby,
  neighbor_mode = "spatial",
  neighbor_graph = NULL,
  neighbor_k = 6L
) {
  mode <- normalize_neighbor_mode(neighbor_mode)
  if (identical(mode, "none")) {
    return(empty_neighbor_info())
  }

  if (!groupby %in% names(obs)) {
    stop("groupby column not found in obs: ", groupby)
  }

  if (mode %in% c("auto", "existing")) {
    graph_info <- resolve_seurat_graph_neighbor_info(
      x = x,
      cell_names = cell_names,
      group_values = obs[[groupby]],
      requested_graph = neighbor_graph
    )
    if (has_neighbor_edges(graph_info$section_edges)) {
      return(graph_info)
    }
    if (identical(mode, "existing")) {
      return(graph_info)
    }
  }

  derive_spatial_neighbor_info(
    coords = coords,
    group_values = obs[[groupby]],
    neighbor_k = neighbor_k
  )
}

resolve_list_neighbor_info <- function(
  x,
  obs,
  coords,
  groupby,
  neighbor_mode = "spatial",
  neighbor_graph = NULL,
  neighbor_k = 6L
) {
  mode <- normalize_neighbor_mode(neighbor_mode)
  if (identical(mode, "none")) {
    return(empty_neighbor_info())
  }

  precomputed <- normalize_precomputed_section_edges(
    section_edges = x$neighbor_edges_by_section %||% x$section_edges %||% NULL
  )
  if (has_neighbor_edges(precomputed)) {
    return(list(
      section_edges = precomputed,
      key = x$neighbors_key %||% neighbor_graph %||% "provided"
    ))
  }
  if (identical(mode, "existing")) {
    return(empty_neighbor_info())
  }

  derive_spatial_neighbor_info(
    coords = coords,
    group_values = obs[[groupby]],
    neighbor_k = neighbor_k
  )
}

empty_neighbor_info <- function() {
  list(
    section_edges = empty_named_list(),
    key = NULL
  )
}

normalize_neighbor_mode <- function(neighbor_mode = "spatial") {
  mode <- tolower(as.character(neighbor_mode %||% "spatial")[[1]])
  if (!(mode %in% c("auto", "existing", "spatial", "none"))) {
    stop(
      "Unsupported neighbor_mode: ",
      neighbor_mode,
      ". Expected one of auto, existing, spatial, none."
    )
  }
  mode
}

has_neighbor_edges <- function(section_edges) {
  if (is.null(section_edges) || length(section_edges) == 0) {
    return(FALSE)
  }
  any(vapply(section_edges, length, integer(1)) > 0L)
}

normalize_precomputed_section_edges <- function(section_edges) {
  if (is.null(section_edges) || length(section_edges) == 0) {
    return(empty_named_list())
  }

  out <- vector("list", length(section_edges))
  names(out) <- names(section_edges)
  for (i in seq_along(section_edges)) {
    edges <- section_edges[[i]]
    if (is.null(edges) || length(edges) == 0) {
      out[[i]] <- integer()
      next
    }
    if (is.matrix(edges) || is.data.frame(edges)) {
      edge_matrix <- as.matrix(edges)
      if (ncol(edge_matrix) < 2) {
        stop("Precomputed section edge matrices must have at least two columns.")
      }
      out[[i]] <- as.integer(as.vector(t(edge_matrix[, seq_len(2), drop = FALSE])))
      next
    }
    out[[i]] <- as.integer(edges)
  }
  out
}

resolve_seurat_graph_neighbor_info <- function(
  x,
  cell_names,
  group_values,
  requested_graph = NULL
) {
  graph_names <- names(x@graphs)
  if (length(graph_names) == 0) {
    if (!is.null(requested_graph) && nzchar(requested_graph)) {
      stop("Requested Seurat graph not found because the object has no graphs: ", requested_graph)
    }
    return(empty_neighbor_info())
  }

  graph_name <- requested_graph
  if (is.null(graph_name) || !nzchar(graph_name)) {
    preferred <- unique(c(
      paste0(x@active.assay, "_snn"),
      paste0(x@active.assay, "_nn"),
      "SCT_snn",
      "SCT_nn",
      "Spatial_snn",
      "Spatial_nn",
      "RNA_snn",
      "RNA_nn"
    ))
    graph_name <- pick_first_existing(graph_names, preferred) %||% graph_names[[1]]
  } else if (!(graph_name %in% graph_names)) {
    stop(
      "Requested Seurat graph not found: ",
      graph_name,
      ". Available graphs: ",
      paste(graph_names, collapse = ", ")
    )
  }

  graph <- x@graphs[[graph_name]]
  list(
    section_edges = extract_sparse_graph_edges_by_section(
      graph = graph,
      cell_names = cell_names,
      group_values = group_values
    ),
    key = graph_name
  )
}

# Extract neighbor graph edges and split them by section (sample/FOV).
#
# Handles two neighbor graph representations:
#   - Seurat-style sparse Matrix stored in x@graphs[[name]]: a cell×cell
#     adjacency matrix where non-zero entries indicate edges.
#   - Other matrix-like objects would raise an error in extract_sparse_graph_edge_pairs.
#
# The graph is first reordered to match cell_names (important when the graph
# was computed on a superset of cells), then converted to an undirected edge
# list as a 2-column integer matrix of global (1-based) cell indices.
#
# Output: named list (one entry per section) of flat integer vectors in the
# format expected by the payload builder — pairs interleaved as
# [src0, dst0, src1, dst1, ...] with 0-based local (within-section) indices.
extract_sparse_graph_edges_by_section <- function(graph, cell_names, group_values) {
  edge_pairs <- extract_sparse_graph_edge_pairs(
    graph = graph,
    cell_names = cell_names
  )
  split_global_edges_by_section(
    edge_pairs = edge_pairs,
    group_values = group_values
  )
}

extract_sparse_graph_edge_pairs <- function(graph, cell_names) {
  if (!methods::is(graph, "Matrix")) {
    stop(
      "Unsupported graph object for neighbor export: ",
      paste(class(graph), collapse = ", ")
    )
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Neighbor graph export from sparse matrices requires the Matrix package.")
  }

  graph_aligned <- graph
  graph_rows <- rownames(graph_aligned)
  graph_cols <- colnames(graph_aligned)
  if (!is.null(graph_rows) && !is.null(graph_cols) &&
      all(cell_names %in% graph_rows) && all(cell_names %in% graph_cols)) {
    graph_aligned <- graph_aligned[cell_names, cell_names, drop = FALSE]
  }

  summary_df <- Matrix::summary(graph_aligned)
  if (nrow(summary_df) == 0L) {
    return(matrix(integer(), ncol = 2L))
  }

  from <- as.integer(summary_df$i)
  to <- as.integer(summary_df$j)
  keep <- from != to
  if (!any(keep)) {
    return(matrix(integer(), ncol = 2L))
  }

  edge_pairs <- cbind(
    from = pmin(from[keep], to[keep]),
    to = pmax(from[keep], to[keep])
  )
  edge_pairs <- edge_pairs[!duplicated(edge_pairs), , drop = FALSE]
  edge_pairs[order(edge_pairs[, 1], edge_pairs[, 2]), , drop = FALSE]
}

# Partition a global edge list into per-section (sample/FOV) local edge lists.
#
# Why split by section: the neighbor graph may span multiple tissue sections
# that are stored together in one object. Cross-section edges are biologically
# meaningless (cells in different slides cannot be spatial neighbors), so they
# are discarded. The payload builder processes each section independently.
#
# Input:
#   edge_pairs   — [n_edges × 2] integer matrix of global 1-based cell indices
#   group_values — character/factor vector (length = n_cells) giving each cell's
#                  section label (e.g. sample_id or orig.ident)
#
# Output: named list (section_id → integer vector). Each vector is a flat,
# interleaved sequence of 0-based local index pairs [src0, dst0, src1, dst1, ...]
# where indices are local to that section (i.e. rank within section, not global).
#
# Variables:
#   same_edges  — subset of edge_pairs where both endpoints share the same section
#   same_labels — section label for each row in same_edges (used to route edges)
#   local_lookup — per-section integer vector: global index → 0-based local index
split_global_edges_by_section <- function(edge_pairs, group_values) {
  section_ids <- unique(as.character(group_values))
  section_edges <- stats::setNames(vector("list", length(section_ids)), section_ids)
  if (length(section_ids) == 0L) {
    return(empty_named_list())
  }

  group_values <- as.character(group_values)
  for (section_id in section_ids) {
    section_edges[[section_id]] <- integer()
  }
  if (is.null(edge_pairs) || length(edge_pairs) == 0L || nrow(edge_pairs) == 0L) {
    return(section_edges)
  }

  # Keep only intra-section edges (discard cross-section edges).
  edge_section <- group_values[edge_pairs[, 1]]
  same_section <- edge_section == group_values[edge_pairs[, 2]]
  if (!any(same_section)) {
    return(section_edges)
  }

  same_edges <- edge_pairs[same_section, , drop = FALSE]
  same_labels <- edge_section[same_section]

  # Build a global→local index map for each section.
  # local_lookup[[s]][global_i] gives the 0-based local position of cell global_i
  # within section s (or 0 / undefined for cells outside section s).
  local_lookup <- lapply(
    section_ids,
    function(section_id) {
      idx <- which(group_values == section_id)
      lookup <- integer(length(group_values))
      lookup[idx] <- seq_along(idx) - 1L
      lookup
    }
  )
  names(local_lookup) <- section_ids

  for (section_id in section_ids) {
    section_pair_idx <- which(same_labels == section_id)
    if (length(section_pair_idx) == 0L) {
      next
    }
    section_pairs <- same_edges[section_pair_idx, , drop = FALSE]
    lookup <- local_lookup[[section_id]]
    # Convert global indices to 0-based local indices, then flatten row-wise
    # into the interleaved format [src0, dst0, src1, dst1, ...].
    local_pairs <- cbind(
      lookup[section_pairs[, 1]] ,
      lookup[section_pairs[, 2]]
    )
    section_edges[[section_id]] <- as.integer(as.vector(t(local_pairs)))
  }

  section_edges
}

derive_spatial_neighbor_info <- function(coords, group_values, neighbor_k = 6L) {
  list(
    section_edges = derive_spatial_neighbor_edges_by_section(
      coords = coords,
      group_values = group_values,
      neighbor_k = neighbor_k
    ),
    key = paste0("spatial_knn_k", as.integer(neighbor_k))
  )
}

derive_spatial_neighbor_edges_by_section <- function(coords, group_values, neighbor_k = 6L) {
  section_ids <- unique(as.character(group_values))
  section_edges <- stats::setNames(vector("list", length(section_ids)), section_ids)
  group_values <- as.character(group_values)

  for (section_id in section_ids) {
    idx <- which(group_values == section_id)
    section_edges[[section_id]] <- derive_section_knn_edges(
      coords = coords[idx, , drop = FALSE],
      neighbor_k = neighbor_k
    )
  }

  section_edges
}

derive_section_knn_edges <- function(coords, neighbor_k = 6L, chunk_size = 512L) {
  coords <- as.matrix(coords)
  n <- nrow(coords)
  if (n < 2L) {
    return(integer())
  }

  k <- max(1L, min(as.integer(neighbor_k), n - 1L))
  if (!is.finite(k) || k < 1L) {
    return(integer())
  }

  norms <- rowSums(coords * coords)
  from <- integer(n * k)
  to <- integer(n * k)
  pos <- 1L

  for (start in seq.int(1L, n, by = chunk_size)) {
    end <- min(start + chunk_size - 1L, n)
    block_idx <- seq.int(start, end)
    block <- coords[block_idx, , drop = FALSE]
    dist_sq <- outer(norms[block_idx], norms, "+") - 2 * (block %*% t(coords))
    dist_sq[cbind(seq_along(block_idx), block_idx)] <- Inf

    for (row_pos in seq_along(block_idx)) {
      nearest <- head(order(dist_sq[row_pos, ]), k)
      row_end <- pos + k - 1L
      from[pos:row_end] <- block_idx[[row_pos]]
      to[pos:row_end] <- nearest
      pos <- row_end + 1L
    }
  }

  if (pos == 1L) {
    return(integer())
  }

  from <- from[seq_len(pos - 1L)] - 1L
  to <- to[seq_len(pos - 1L)] - 1L
  undirected <- cbind(
    from = pmin(from, to),
    to = pmax(from, to)
  )
  undirected <- undirected[undirected[, 1] != undirected[, 2], , drop = FALSE]
  if (nrow(undirected) == 0L) {
    return(integer())
  }

  undirected <- undirected[!duplicated(undirected), , drop = FALSE]
  undirected <- undirected[order(undirected[, 1], undirected[, 2]), , drop = FALSE]
  as.integer(as.vector(t(undirected)))
}

normalize_expression <- function(expression, gene_names = NULL, n_cells, cell_names = NULL) {
  if (is.null(expression)) {
    return(list(expression = NULL, gene_names = character()))
  }

  expr <- as_plain_matrix(expression)
  dims <- dim(expr)
  if (length(dims) != 2) {
    stop("expression must be two-dimensional.")
  }

  if (dims[[2]] == n_cells) {
    if (!is.null(cell_names) && !is.null(colnames(expr)) && all(cell_names %in% colnames(expr))) {
      expr <- expr[, cell_names, drop = FALSE]
    }
    inferred_genes <- coalesce_character(gene_names, rownames(expr))
    if (length(inferred_genes) == 0) {
      inferred_genes <- sprintf("Gene%05d", seq_len(dims[[1]]))
    }
    return(list(
      expression = expr,
      gene_names = as.character(inferred_genes)
    ))
  }

  if (dims[[1]] == n_cells) {
    if (!is.null(cell_names) && !is.null(rownames(expr)) && all(cell_names %in% rownames(expr))) {
      expr <- expr[cell_names, , drop = FALSE]
    }
    inferred_genes <- coalesce_character(gene_names, colnames(expr))
    if (length(inferred_genes) == 0) {
      inferred_genes <- sprintf("Gene%05d", seq_len(dims[[2]]))
    }
    return(list(
      expression = t(expr),
      gene_names = as.character(inferred_genes)
    ))
  }

  stop("expression dimensions must align with the number of cells.")
}

resolve_seurat_expression <- function(x, requested_assay = NULL) {
  assay_name <- resolve_seurat_assay_name(
    x = x,
    requested_assay = requested_assay
  )
  assay_obj <- x@assays[[assay_name]]
  expression_info <- extract_expression_from_assay_object(assay_obj)
  if (is.null(expression_info$expression)) {
    stop(
      "Could not resolve expression data from Seurat assay '",
      assay_name,
      "'."
    )
  }
  expression_info$assay_name <- assay_name
  expression_info
}

resolve_input_gene_names <- function(x, requested_assay = NULL) {
  if (inherits(x, "Seurat")) {
    assay_name <- resolve_seurat_assay_name(x = x, requested_assay = requested_assay)
    assay_obj <- x@assays[[assay_name]]
    gene_names <- extract_assay_feature_names(assay_obj)
    return(list(
      assay_name = assay_name,
      gene_names = unique(as.character(gene_names))
    ))
  }

  if (inherits(x, "SpatialExperiment") || inherits(x, "SingleCellExperiment")) {
    expression_info <- resolve_summarized_experiment_expression(
      x = x,
      requested_assay = requested_assay
    )
    return(list(
      assay_name = expression_info$assay_name,
      gene_names = unique(as.character(expression_info$gene_names))
    ))
  }

  if (is.list(x)) {
    expression <- x$expression %||% x$expr
    expression_info <- if (is.null(expression)) {
      list(expression = NULL, gene_names = character())
    } else {
      normalize_expression(
        expression = expression,
        gene_names = x$gene_names %||% x$genes_available %||% rownames(expression),
        n_cells = nrow(x$obs %||% data.frame())
      )
    }
    return(list(
      assay_name = requested_assay,
      gene_names = unique(as.character(expression_info$gene_names))
    ))
  }

  stop(
    "Could not resolve gene names for input class: ",
    paste(class(x), collapse = ", ")
  )
}

resolve_seurat_assay_name <- function(x, requested_assay = NULL) {
  assay_names <- names(x@assays)
  if (length(assay_names) == 0) {
    return(NULL)
  }

  if (!is.null(requested_assay) && nzchar(requested_assay)) {
    if (!(requested_assay %in% assay_names)) {
      stop(
        "Requested assay not found in Seurat object: ",
        requested_assay,
        ". Available assays: ",
        paste(assay_names, collapse = ", ")
      )
    }
    return(requested_assay)
  }

  preferred <- unique(c(
    "SCT",
    "Spatial",
    x@active.assay,
    "RNA",
    "integrated"
  ))
  chosen <- pick_first_existing(assay_names, preferred)
  chosen %||% assay_names[[1]]
}

resolve_summarized_experiment_expression <- function(x, requested_assay = NULL) {
  assay_names <- SummarizedExperiment::assayNames(x)
  if (length(assay_names) == 0) {
    return(list(
      expression = NULL,
      gene_names = character(),
      assay_name = NULL
    ))
  }

  assay_name <- requested_assay
  if (is.null(assay_name) || !nzchar(assay_name)) {
    assay_name <- pick_first_existing(
      assay_names,
      c("logcounts", "data", "normalized", "counts")
    ) %||% assay_names[[1]]
  }

  if (!(assay_name %in% assay_names)) {
    stop(
      "Requested assay not found: ",
      assay_name,
      ". Available assays: ",
      paste(assay_names, collapse = ", ")
    )
  }

  expression <- SummarizedExperiment::assay(x, assay_name)
  list(
    expression = expression,
    gene_names = coalesce_character(rownames(x), rownames(expression)),
    assay_name = assay_name
  )
}

extract_expression_from_assay_object <- function(assay_obj) {
  feature_names <- extract_assay_feature_names(assay_obj)

  if ("layers" %in% slotNames(assay_obj)) {
    layer_names <- names(assay_obj@layers)
    layer_name <- pick_first_existing(layer_names, c("data", "counts")) %||%
      if (length(layer_names) > 0) layer_names[[1]] else NULL
    if (!is.null(layer_name)) {
      expression <- assay_obj@layers[[layer_name]]
      dims <- dim(expression)
      if (!is.null(dims) && length(dims) == 2 && all(dims > 0)) {
        return(list(
          expression = expression,
          gene_names = coalesce_character(rownames(expression), feature_names),
          source = paste0("layer:", layer_name)
        ))
      }
    }
  }

  for (slot_name in c("data", "counts")) {
    if (!(slot_name %in% slotNames(assay_obj))) {
      next
    }
    expression <- slot(assay_obj, slot_name)
    dims <- dim(expression)
    if (!is.null(dims) && length(dims) == 2 && all(dims > 0)) {
      return(list(
        expression = expression,
        gene_names = coalesce_character(rownames(expression), feature_names),
        source = paste0("slot:", slot_name)
      ))
    }
  }

  list(
    expression = NULL,
    gene_names = character(),
    source = NULL
  )
}

extract_assay_feature_names <- function(assay_obj) {
  if ("features" %in% slotNames(assay_obj)) {
    feature_names <- rownames(assay_obj@features)
    if (length(feature_names) > 0) {
      return(as.character(feature_names))
    }
  }

  if (!is.null(rownames(assay_obj)) && length(rownames(assay_obj)) > 0) {
    return(as.character(rownames(assay_obj)))
  }

  character()
}

rank_top_expressed_gene_names <- function(expression, gene_names, top_genes_n = 20L) {
  gene_names <- as.character(gene_names %||% character())
  if (length(gene_names) == 0L || is.null(expression)) {
    return(character())
  }

  top_genes_n <- suppressWarnings(as.integer(top_genes_n %||% 20L))
  if (is.na(top_genes_n) || top_genes_n < 1L) {
    return(character())
  }

  mean_expression <- matrix_rowmeans(expression)
  if (length(mean_expression) != length(gene_names)) {
    stop("Expression matrix rows do not align with gene names.")
  }

  valid <- is.finite(mean_expression) & !is.na(gene_names) & nzchar(gene_names)
  if (!any(valid)) {
    return(character())
  }

  valid_idx <- which(valid)
  ord <- order(
    -mean_expression[valid_idx],
    gene_names[valid_idx]
  )
  gene_names[valid_idx][utils::head(ord, top_genes_n)]
}

resolve_selected_genes <- function(genes, gene_names, expression = NULL, top_genes_n = NULL) {
  gene_names <- as.character(gene_names %||% character())
  if (length(gene_names) == 0) {
    return(character())
  }

  if (!is.null(genes) && length(genes) > 0) {
    selected <- intersect(as.character(genes), gene_names)
    if (length(selected) == 0) {
      stop("None of the requested genes were found in the expression matrix.")
    }
    return(selected)
  }

  if (!is.null(top_genes_n) && length(top_genes_n) > 0) {
    top_selected <- rank_top_expressed_gene_names(
      expression = expression,
      gene_names = gene_names,
      top_genes_n = top_genes_n
    )
    if (length(top_selected) == 0L) {
      stop("Could not resolve top expressed genes from the expression matrix.")
    }
    return(top_selected)
  }

  gene_names[seq_len(min(20, length(gene_names)))]
}

resolve_metadata_columns <- function(obs, groupby, metadata_columns = NULL) {
  if (!is.null(metadata_columns) && length(metadata_columns) > 0) {
    missing <- setdiff(metadata_columns, names(obs))
    if (length(missing) > 0) {
      stop("metadata_columns not found in obs: ", paste(missing, collapse = ", "))
    }
    return(as.character(metadata_columns))
  }

  is_metadata <- vapply(obs, is_categorical, logical(1))
  metadata <- setdiff(names(obs)[is_metadata], groupby)
  if (length(metadata) == 0) {
    return(character())
  }

  group_values <- as.character(obs[[groupby]])
  section_ids <- unique(group_values)
  keep <- vapply(
    metadata,
    function(column_name) {
      all(vapply(
        section_ids,
        function(section_id) {
          idx <- which(group_values == section_id)
          values <- unique_non_missing(obs[[column_name]][idx])
          length(values) <= 1
        },
        logical(1)
      ))
    },
    logical(1)
  )

  metadata[keep]
}

pick_first_existing <- function(choices, preferred) {
  matches <- preferred[preferred %in% choices]
  if (length(matches) == 0) {
    return(NULL)
  }
  matches[[1]]
}
