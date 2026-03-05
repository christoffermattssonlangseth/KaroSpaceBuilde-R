build_viewer_payload <- function(
  input,
  groupby,
  initial_color,
  additional_colors = NULL,
  genes = NULL,
  assay = NULL,
  metadata_input = NULL,
  metadata_input_columns = NULL,
  metadata_prefix = NULL,
  metadata_columns = NULL,
  outline_by = NULL
) {
  source <- normalize_input_source(
    x = prepare_karospace_input(
      input = input,
      metadata_input = metadata_input,
      metadata_input_columns = metadata_input_columns,
      metadata_prefix = metadata_prefix
    ),
    groupby = groupby,
    initial_color = initial_color,
    additional_colors = additional_colors,
    genes = genes,
    assay = assay,
    metadata_columns = metadata_columns,
    outline_by = outline_by
  )

  build_payload_from_normalized(source)
}

build_payload_from_normalized <- function(source) {
  obs <- source$obs
  coords <- source$coords
  umap <- source$umap
  expression <- source$expression
  group_values <- as.character(obs[[source$groupby]])
  section_ids <- unique(group_values)

  color_columns <- unique(source$additional_colors)
  color_data <- lapply(color_columns, function(col) build_color_column(obs[[col]]))
  names(color_data) <- color_columns

  gene_data <- build_gene_data(
    expression = expression,
    gene_names = source$gene_names,
    selected_genes = source$selected_genes
  )

  metadata_filters <- build_metadata_filters(
    obs = obs,
    metadata_columns = source$metadata_columns
  )

  sections <- vector("list", length(section_ids))
  all_indices <- seq_len(nrow(obs)) - 1L
  has_umap <- !is.null(umap)

  for (i in seq_along(section_ids)) {
    section_id <- section_ids[[i]]
    idx <- which(group_values == section_id)
    section_obs <- obs[idx, , drop = FALSE]
    section_coords <- coords[idx, , drop = FALSE]

    section_colors <- lapply(color_data, function(info) {
      as_json_array(as.numeric(info$values[idx]))
    })

    section_genes <- lapply(gene_data$values, function(values) {
      as_json_array(as.numeric(values[idx]))
    })

    section_entry <- list(
      id = section_id,
      metadata = build_section_metadata(section_obs, source$metadata_columns),
      n_cells = length(idx),
      x = as_json_array(as.numeric(section_coords[, 1])),
      y = as_json_array(as.numeric(section_coords[, 2])),
      xb64 = NULL,
      yb64 = NULL,
      obs_idx = as_json_array(as.integer(all_indices[idx])),
      obs_idxb64 = NULL,
      colors = section_colors,
      colors_b64 = empty_named_list(),
      genes = section_genes,
      genes_sparse = empty_named_list(),
      bounds = list(
        xmin = min(section_coords[, 1]),
        xmax = max(section_coords[, 1]),
        ymin = min(section_coords[, 2]),
        ymax = max(section_coords[, 2])
      ),
      umap_x = NULL,
      umap_y = NULL,
      umap_xb64 = NULL,
      umap_yb64 = NULL,
      edges = list(),
      edges_b64 = NULL
    )

    if (has_umap) {
      section_entry$umap_x <- as_json_array(as.numeric(umap[idx, 1]))
      section_entry$umap_y <- as_json_array(as.numeric(umap[idx, 2]))
    }

    sections[[i]] <- section_entry
  }

  colors_meta <- lapply(color_data, function(info) info$meta)
  genes_meta <- lapply(gene_data$meta, identity)
  gene_encodings <- stats::setNames(
    as.list(rep("dense", length(gene_data$meta))),
    names(gene_data$meta)
  )

  payload <- list(
    schema_version = "1.0.0",
    initial_color = source$initial_color,
    groupby = source$groupby,
    colors_meta = colors_meta,
    genes_meta = genes_meta,
    gene_encodings = gene_encodings,
    metadata_filters = metadata_filters,
    n_sections = length(sections),
    total_cells = nrow(obs),
    sections = sections,
    available_colors = as_json_array(color_columns),
    available_genes = as_json_array(names(gene_data$meta)),
    marker_genes = empty_named_list(),
    has_umap = has_umap,
    umap_bounds = build_umap_bounds(umap),
    has_neighbors = FALSE,
    neighbors_key = NULL,
    neighbor_stats = empty_named_list(),
    interaction_markers = empty_named_list()
  )

  if (!is.null(source$outline_by)) {
    payload$outline_by <- source$outline_by
  }

  payload
}

build_color_column <- function(column) {
  if (is.numeric(column)) {
    values <- as.numeric(column)
    finite <- values[is.finite(values)]
    meta <- list(
      is_continuous = TRUE,
      categories = NULL,
      vmin = if (length(finite) > 0) min(finite) else 0,
      vmax = if (length(finite) > 0) max(finite) else 1
    )
    return(list(values = values, meta = meta))
  }

  factor_column <- if (is.factor(column)) column else factor(as.character(column))
  values <- as.numeric(factor_column) - 1
  values[is.na(values)] <- NA_real_
  meta <- list(
    is_continuous = FALSE,
    categories = as_json_array(as.character(levels(factor_column))),
    vmin = 0,
    vmax = max(0, length(levels(factor_column)) - 1)
  )

  list(values = values, meta = meta)
}

build_gene_data <- function(expression, gene_names, selected_genes) {
  if (is.null(expression) || length(selected_genes) == 0) {
    return(list(values = empty_named_list(), meta = empty_named_list()))
  }

  gene_names <- as.character(gene_names)
  gene_index <- match(selected_genes, gene_names)

  values <- vector("list", length(selected_genes))
  meta <- vector("list", length(selected_genes))
  names(values) <- selected_genes
  names(meta) <- selected_genes

  for (i in seq_along(selected_genes)) {
    gene <- selected_genes[[i]]
    idx <- gene_index[[i]]
    vector <- extract_gene_vector(expression, idx)
    finite <- vector[is.finite(vector)]

    values[[gene]] <- vector
    meta[[gene]] <- list(
      vmin = if (length(finite) > 0) min(finite) else 0,
      vmax = if (length(finite) > 0) max(finite) else 1
    )
  }

  list(values = values, meta = meta)
}

extract_gene_vector <- function(expression, idx) {
  if (methods::is(expression, "Matrix")) {
    return(as.numeric(expression[idx, , drop = TRUE]))
  }
  as.numeric(expression[idx, , drop = TRUE])
}

build_metadata_filters <- function(obs, metadata_columns) {
  filters <- empty_named_list()
  for (column_name in metadata_columns %||% character()) {
    values <- unique_non_missing(obs[[column_name]])
    if (length(values) == 0) {
      next
    }
    filters[[column_name]] <- as_json_array(sort(values))
  }
  filters
}

build_section_metadata <- function(obs, metadata_columns) {
  metadata <- empty_named_list()
  for (column_name in metadata_columns %||% character()) {
    value <- first_non_missing(obs[[column_name]])
    if (!is.null(value)) {
      metadata[[column_name]] <- value
    }
  }
  metadata
}

build_umap_bounds <- function(umap) {
  if (is.null(umap)) {
    return(NULL)
  }

  list(
    xmin = min(umap[, 1]),
    xmax = max(umap[, 1]),
    ymin = min(umap[, 2]),
    ymax = max(umap[, 2])
  )
}
