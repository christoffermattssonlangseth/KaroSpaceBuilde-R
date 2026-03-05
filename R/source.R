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

normalize_input_source <- function(
  x,
  groupby,
  initial_color,
  additional_colors = NULL,
  genes = NULL,
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
  metadata_columns = NULL,
  outline_by = NULL
) {
  obs <- x@meta.data
  if (!is.data.frame(obs)) {
    stop("Seurat object does not have a usable meta.data frame.")
  }

  cell_names <- rownames(obs)
  if (is.null(cell_names) || length(cell_names) == 0) {
    stop("Seurat meta.data must have row names matching cell names.")
  }

  if (length(x@images) == 0) {
    stop("Seurat object has no spatial image data in @images.")
  }

  coord_frames <- list()
  x_col <- NULL
  y_col <- NULL
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

    x_col <- x_col %||% image_x_col
    y_col <- y_col %||% image_y_col
    coord_frames[[image_name]] <- image_coords
  }

  if (length(coord_frames) == 0) {
    stop("Could not resolve coordinates from any Seurat image slot.")
  }

  coord_df <- do.call(rbind, unname(coord_frames))
  coord_df <- coord_df[!duplicated(rownames(coord_df)), , drop = FALSE]
  missing_cells <- setdiff(cell_names, rownames(coord_df))
  if (length(missing_cells) > 0) {
    stop(
      "Seurat image coordinates do not cover all cells in meta.data. Missing ",
      length(missing_cells),
      " cells."
    )
  }
  coord_df <- coord_df[cell_names, , drop = FALSE]

  coords <- as.matrix(coord_df[, c(x_col, y_col), drop = FALSE])
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

  assay_name <- if ("Spatial" %in% names(x@assays)) "Spatial" else x@active.assay
  assay_obj <- x@assays[[assay_name]]
  expression <- NULL
  if ("data" %in% slotNames(assay_obj) && nrow(assay_obj@data) > 0) {
    expression <- assay_obj@data
  } else if ("counts" %in% slotNames(assay_obj) && nrow(assay_obj@counts) > 0) {
    expression <- assay_obj@counts
  }

  normalize_list_source(
    x = list(
      obs = obs,
      coordinates = coords,
      umap = umap,
      expression = expression,
      gene_names = if (!is.null(expression)) rownames(expression) else character()
    ),
    groupby = groupby,
    initial_color = initial_color,
    additional_colors = additional_colors,
    genes = genes,
    metadata_columns = metadata_columns,
    outline_by = outline_by
  )
}

normalize_single_cell_experiment <- function(
  x,
  groupby,
  initial_color,
  additional_colors = NULL,
  genes = NULL,
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

  assay_names <- SummarizedExperiment::assayNames(x)
  assay_name <- if ("logcounts" %in% assay_names) "logcounts" else assay_names[[1]]
  expression <- SummarizedExperiment::assay(x, assay_name)

  normalize_list_source(
    x = list(
      obs = obs,
      coordinates = spatial,
      umap = umap,
      expression = expression,
      gene_names = rownames(x)
    ),
    groupby = groupby,
    initial_color = initial_color,
    additional_colors = additional_colors,
    genes = genes,
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

  assay_names <- SummarizedExperiment::assayNames(x)
  assay_name <- if ("logcounts" %in% assay_names) "logcounts" else assay_names[[1]]
  expression <- SummarizedExperiment::assay(x, assay_name)

  normalize_list_source(
    x = list(
      obs = obs,
      coordinates = spatial,
      umap = umap,
      expression = expression,
      gene_names = rownames(x)
    ),
    groupby = groupby,
    initial_color = initial_color,
    additional_colors = additional_colors,
    genes = genes,
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
    n_cells = nrow(obs)
  )

  list(
    obs = obs,
    coords = coords,
    umap = umap,
    expression = expression_info$expression,
    gene_names = expression_info$gene_names,
    selected_genes = resolve_selected_genes(
      genes = genes,
      gene_names = expression_info$gene_names
    ),
    groupby = groupby,
    initial_color = initial_color,
    additional_colors = unique(c(initial_color, additional_colors)),
    metadata_columns = resolve_metadata_columns(
      obs = obs,
      groupby = groupby,
      metadata_columns = metadata_columns
    ),
    outline_by = outline_by
  )
}

normalize_expression <- function(expression, gene_names = NULL, n_cells) {
  if (is.null(expression)) {
    return(list(expression = NULL, gene_names = character()))
  }

  expr <- as_plain_matrix(expression)
  dims <- dim(expr)
  if (length(dims) != 2) {
    stop("expression must be two-dimensional.")
  }

  if (dims[[2]] == n_cells) {
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

resolve_selected_genes <- function(genes, gene_names) {
  gene_names <- as.character(gene_names %||% character())
  if (length(gene_names) == 0) {
    return(character())
  }

  if (is.null(genes) || length(genes) == 0) {
    return(gene_names[seq_len(min(20, length(gene_names)))])
  }

  selected <- intersect(as.character(genes), gene_names)
  if (length(selected) == 0) {
    stop("None of the requested genes were found in the expression matrix.")
  }
  selected
}

resolve_metadata_columns <- function(obs, groupby, metadata_columns = NULL) {
  if (!is.null(metadata_columns) && length(metadata_columns) > 0) {
    missing <- setdiff(metadata_columns, names(obs))
    if (length(missing) > 0) {
      stop("metadata_columns not found in obs: ", paste(missing, collapse = ", "))
    }
    return(as.character(metadata_columns))
  }

  is_metadata <- vapply(
    obs,
    function(column) is.factor(column) || is.character(column) || is.logical(column),
    logical(1)
  )
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
