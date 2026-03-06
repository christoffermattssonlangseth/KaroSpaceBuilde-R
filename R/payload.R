build_viewer_payload <- function(
  input,
  groupby,
  initial_color,
  additional_colors = NULL,
  genes = NULL,
  top_genes_n = NULL,
  assay = NULL,
  metadata_input = NULL,
  metadata_input_columns = NULL,
  metadata_prefix = NULL,
  neighbor_mode = "spatial",
  neighbor_graph = NULL,
  neighbor_k = 6L,
  metadata_columns = NULL,
  outline_by = NULL,
  lightweight = FALSE,
  pack_arrays = TRUE,
  pack_arrays_min_len = 1000L,
  gene_sparse_zero_threshold = 0.8,
  gene_sparse_pack = TRUE,
  gene_sparse_pack_min_nnz = 32L,
  marker_genes_groupby = NULL,
  marker_genes_top_n = 20L,
  interaction_markers_groupby = NULL,
  interaction_markers_top_targets = 8L,
  interaction_markers_top_genes = 12L,
  interaction_markers_min_cells = 30L,
  interaction_markers_min_neighbors = 1L,
  marker_test = "mean_diff",
  neighbor_stats_permutations = 0L,
  neighbor_stats_seed = 42L
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
    top_genes_n = top_genes_n,
    assay = assay,
    neighbor_mode = neighbor_mode,
    neighbor_graph = neighbor_graph,
    neighbor_k = neighbor_k,
    metadata_columns = metadata_columns,
    outline_by = outline_by
  )

  build_payload_from_normalized(
    source = source,
    lightweight = lightweight,
    pack_arrays = pack_arrays,
    pack_arrays_min_len = pack_arrays_min_len,
    gene_sparse_zero_threshold = gene_sparse_zero_threshold,
    gene_sparse_pack = gene_sparse_pack,
    gene_sparse_pack_min_nnz = gene_sparse_pack_min_nnz,
    marker_genes_groupby = marker_genes_groupby,
    marker_genes_top_n = marker_genes_top_n,
    interaction_markers_groupby = interaction_markers_groupby,
    interaction_markers_top_targets = interaction_markers_top_targets,
    interaction_markers_top_genes = interaction_markers_top_genes,
    interaction_markers_min_cells = interaction_markers_min_cells,
    interaction_markers_min_neighbors = interaction_markers_min_neighbors,
    marker_test = marker_test,
    neighbor_stats_permutations = neighbor_stats_permutations,
    neighbor_stats_seed = neighbor_stats_seed
  )
}

build_payload_from_normalized <- function(
  source,
  lightweight = FALSE,
  pack_arrays = TRUE,
  pack_arrays_min_len = 1000L,
  gene_sparse_zero_threshold = 0.8,
  gene_sparse_pack = TRUE,
  gene_sparse_pack_min_nnz = 32L,
  marker_genes_groupby = NULL,
  marker_genes_top_n = 20L,
  interaction_markers_groupby = NULL,
  interaction_markers_top_targets = 8L,
  interaction_markers_top_genes = 12L,
  interaction_markers_min_cells = 30L,
  interaction_markers_min_neighbors = 1L,
  marker_test = "mean_diff",
  neighbor_stats_permutations = 0L,
  neighbor_stats_seed = 42L
) {
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
    selected_genes = source$selected_genes,
    sparse_zero_threshold = gene_sparse_zero_threshold
  )

  metadata_filters <- build_metadata_filters(
    obs = obs,
    metadata_columns = source$metadata_columns
  )

  sections <- vector("list", length(section_ids))
  section_index_map <- stats::setNames(vector("list", length(section_ids)), section_ids)
  all_indices <- seq_len(nrow(obs)) - 1L
  has_umap <- !is.null(umap)
  section_edges <- source$section_edges %||% empty_named_list()

  for (i in seq_along(section_ids)) {
    section_id <- section_ids[[i]]
    idx <- which(group_values == section_id)
    section_index_map[[section_id]] <- idx
    section_obs <- obs[idx, , drop = FALSE]
    section_coords <- coords[idx, , drop = FALSE]
    section_edge_pairs <- as.integer(section_edges[[section_id]] %||% integer())
    should_pack_arrays <- isTRUE(pack_arrays) && length(idx) >= as.integer(pack_arrays_min_len %||% 1000L)

    section_colors <- empty_named_list()
    section_colors_b64 <- empty_named_list()
    for (color_name in names(color_data)) {
      section_values <- as.numeric(color_data[[color_name]]$values[idx])
      if (should_pack_arrays) {
        section_colors_b64[[color_name]] <- pack_float32_base64(section_values)
      } else {
        section_colors[[color_name]] <- as_json_array(section_values)
      }
    }

    section_genes <- empty_named_list()
    section_genes_sparse <- empty_named_list()
    for (gene_name in names(gene_data$values)) {
      section_values <- as.numeric(gene_data$values[[gene_name]][idx])
      mode <- gene_data$encodings[[gene_name]] %||% "dense"
      if (identical(mode, "sparse")) {
        sparse_entry <- encode_sparse_gene_values(
          values = section_values,
          pack_sparse = gene_sparse_pack,
          pack_min_nnz = gene_sparse_pack_min_nnz
        )
        if (!is.null(sparse_entry)) {
          section_genes_sparse[[gene_name]] <- sparse_entry
        }
      } else {
        section_genes[[gene_name]] <- as_json_array(section_values)
      }
    }

    packed_edges <- pack_uint32_base64(section_edge_pairs)

    section_entry <- list(
      id = section_id,
      metadata = build_section_metadata(section_obs, source$metadata_columns),
      n_cells = length(idx),
      x = if (should_pack_arrays) NULL else as_json_array(as.numeric(section_coords[, 1])),
      y = if (should_pack_arrays) NULL else as_json_array(as.numeric(section_coords[, 2])),
      xb64 = if (should_pack_arrays) pack_float32_base64(section_coords[, 1]) else NULL,
      yb64 = if (should_pack_arrays) pack_float32_base64(section_coords[, 2]) else NULL,
      obs_idx = if (should_pack_arrays) NULL else as_json_array(as.integer(all_indices[idx])),
      obs_idxb64 = if (should_pack_arrays) pack_uint32_base64(as.integer(all_indices[idx])) else NULL,
      colors = section_colors,
      colors_b64 = section_colors_b64,
      genes = section_genes,
      genes_sparse = section_genes_sparse,
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
      edges = if (is.null(packed_edges)) list() else NULL,
      edges_b64 = packed_edges
    )

    if (has_umap) {
      section_entry$umap_x <- if (should_pack_arrays) NULL else as_json_array(as.numeric(umap[idx, 1]))
      section_entry$umap_y <- if (should_pack_arrays) NULL else as_json_array(as.numeric(umap[idx, 2]))
      section_entry$umap_xb64 <- if (should_pack_arrays) pack_float32_base64(umap[idx, 1]) else NULL
      section_entry$umap_yb64 <- if (should_pack_arrays) pack_float32_base64(umap[idx, 2]) else NULL
    }

    sections[[i]] <- section_entry
  }

  colors_meta <- lapply(color_data, function(info) info$meta)
  genes_meta <- lapply(gene_data$meta, identity)
  gene_encodings <- gene_data$encodings
  neighbor_stats <- build_neighbor_stats(
    color_data = color_data,
    section_index_map = section_index_map,
    section_edges = section_edges,
    n_perms = neighbor_stats_permutations,
    seed = neighbor_stats_seed
  )
  has_neighbors <- has_neighbor_edges(section_edges)
  if (isTRUE(lightweight)) {
    if (length(marker_genes_groupby %||% character()) == 0L) {
      marker_genes_groupby <- "none"
    }
    if (length(interaction_markers_groupby %||% character()) == 0L) {
      interaction_markers_groupby <- "none"
    }
  }
  marker_columns <- resolve_phase4_groupby_columns(
    requested_groupby = unique(c(marker_genes_groupby, interaction_markers_groupby)),
    color_data = color_data
  )
  marker_genes <- build_marker_genes(
    expression = expression,
    gene_names = source$gene_names,
    color_data = color_data,
    requested_groupby = marker_columns,
    top_n = marker_genes_top_n,
    marker_test = marker_test
  )
  interaction_markers <- build_interaction_markers(
    expression = expression,
    gene_names = source$gene_names,
    color_data = color_data,
    requested_groupby = interaction_markers_groupby,
    neighbor_stats = neighbor_stats,
    section_index_map = section_index_map,
    section_edges = section_edges,
    top_targets = interaction_markers_top_targets,
    top_genes = interaction_markers_top_genes,
    min_cells = interaction_markers_min_cells,
    min_neighbors = interaction_markers_min_neighbors,
    marker_test = marker_test
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
    marker_genes = marker_genes,
    has_umap = has_umap,
    umap_bounds = build_umap_bounds(umap),
    has_neighbors = has_neighbors,
    neighbors_key = if (has_neighbors) source$neighbors_key %||% "neighbors" else NULL,
    neighbor_stats = neighbor_stats,
    interaction_markers = interaction_markers
  )

  if (!is.null(source$outline_by)) {
    payload$outline_by <- source$outline_by
  }

  payload
}

build_neighbor_stats <- function(color_data, section_index_map, section_edges, n_perms = 0L, seed = 42L) {
  stats <- empty_named_list()
  if (!has_neighbor_edges(section_edges) || length(color_data) == 0L) {
    return(stats)
  }

  n_perms <- max(0L, as.integer(n_perms %||% 0L))

  for (color_name in names(color_data)) {
    info <- color_data[[color_name]]
    if (isTRUE(info$meta$is_continuous)) {
      next
    }

    categories <- unlist(info$meta$categories, use.names = FALSE)
    n_categories <- length(categories)
    if (n_categories == 0L) {
      next
    }

    encoded_all <- as.integer(info$values) + 1L
    n_cells <- tabulate(encoded_all, nbins = n_categories)
    counts <- matrix(0, nrow = n_categories, ncol = n_categories)
    degree_sum <- numeric(n_categories)
    n_cats_sq <- n_categories^2L

    for (section_id in names(section_index_map)) {
      idx <- section_index_map[[section_id]]
      if (length(idx) == 0L) {
        next
      }
      packed_edges <- as.integer(section_edges[[section_id]] %||% integer())
      if (length(packed_edges) == 0L) {
        next
      }
      edge_matrix <- matrix(packed_edges, ncol = 2L, byrow = TRUE)
      section_values <- encoded_all[idx]
      source_cat <- section_values[edge_matrix[, 1] + 1L]
      target_cat <- section_values[edge_matrix[, 2] + 1L]

      valid <- is.finite(source_cat) & is.finite(target_cat) &
               source_cat >= 1L & source_cat <= n_categories &
               target_cat >= 1L & target_cat <= n_categories
      if (!any(valid)) {
        next
      }
      sc <- as.integer(source_cat[valid])
      tc <- as.integer(target_cat[valid])

      combined_fwd <- (sc - 1L) * n_categories + tc
      combined_rev <- (tc - 1L) * n_categories + sc
      counts <- counts +
        matrix(tabulate(combined_fwd, nbins = n_cats_sq), n_categories, n_categories) +
        matrix(tabulate(combined_rev, nbins = n_cats_sq), n_categories, n_categories)
      degree_sum <- degree_sum +
        tabulate(sc, nbins = n_categories) +
        tabulate(tc, nbins = n_categories)
    }

    z <- if (n_perms > 0L) {
      compute_neighbor_zscores(
        encoded_all = encoded_all,
        n_categories = n_categories,
        section_index_map = section_index_map,
        section_edges = section_edges,
        observed_counts = counts,
        n_perms = n_perms,
        seed = seed
      )
    } else {
      NULL
    }

    mean_degree <- ifelse(n_cells > 0, degree_sum / n_cells, 0)
    stats[[color_name]] <- list(
      categories = as_json_array(categories),
      counts = lapply(seq_len(n_categories), function(i) as_json_array(as.numeric(counts[i, ]))),
      n_cells = as_json_array(as.integer(n_cells)),
      mean_degree = as_json_array(as.numeric(mean_degree)),
      zscore = if (!is.null(z)) {
        lapply(seq_len(n_categories), function(i) as_json_array(as.numeric(z[i, ])))
      } else {
        NULL
      },
      perm_n = as.integer(n_perms)
    )
  }

  stats
}

# Estimate neighbor enrichment z-scores via permutation testing.
#
# The approach: for each permutation, cell-type labels are shuffled
# independently within each section (preserving per-section composition),
# and a symmetric co-occurrence count matrix [n_categories x n_categories]
# is tallied over all edges. The mean and variance of this null distribution
# are tracked online using Welford's algorithm (perm_mean / perm_M2).
#
# Final z-score = (observed - null_mean) / null_sd
#   - positive z: category pair co-occurs more than expected by chance
#   - negative z: category pair avoids each other
#   - NA: null_sd == 0 (e.g. only one category present in every section)
#
# Args:
#   encoded_all        — integer vector of category codes (1..n_categories) for all cells
#   n_categories       — number of distinct categories
#   section_index_map  — named list: section_id → global cell indices in that section
#   section_edges      — named list: section_id → flat edge index pairs (0-based local)
#   observed_counts    — [n_categories x n_categories] observed co-occurrence matrix
#   n_perms            — number of label-shuffle permutations
#   seed               — RNG seed (caller's RNG state is preserved on exit)
compute_neighbor_zscores <- function(
  encoded_all, n_categories, section_index_map, section_edges,
  observed_counts, n_perms = 200L, seed = 42L
) {
  # Unpack section edges into global-index edge lists; also record which cells
  # belong to each section (needed for section-preserving shuffles below).
  section_idx_list <- list()
  all_src <- integer(0)
  all_dst <- integer(0)
  for (section_id in names(section_index_map)) {
    idx <- section_index_map[[section_id]]
    packed <- as.integer(section_edges[[section_id]] %||% integer())
    if (!length(packed) || !length(idx)) next
    # packed is a flat interleaved vector: [src0, dst0, src1, dst1, ...]
    # indices are 0-based local; +1 converts to 1-based for idx lookup.
    em <- matrix(packed, ncol = 2L, byrow = TRUE) + 1L
    all_src <- c(all_src, idx[em[, 1]])
    all_dst <- c(all_dst, idx[em[, 2]])
    section_idx_list[[section_id]] <- idx
  }
  if (!length(all_src)) return(NULL)

  n_cats_sq <- n_categories^2L

  # Preserve caller's RNG state so permutation testing is side-effect-free.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  # Welford online accumulators for the null distribution.
  perm_mean <- matrix(0, n_categories, n_categories)
  perm_M2   <- matrix(0, n_categories, n_categories)

  for (k in seq_len(n_perms)) {
    # Section-preserving shuffle: permute labels within each section independently.
    # This holds section composition fixed — the correct null for a block-diagonal graph.
    pl <- encoded_all
    for (idx in section_idx_list) {
      pl[idx] <- sample(encoded_all[idx])
    }
    # Symmetrised co-occurrence: count each undirected edge in both directions
    # so the result matrix is symmetric (same as observed_counts).
    sc <- pl[all_src]
    dc <- pl[all_dst]
    valid <- is.finite(sc) & is.finite(dc) &
             sc >= 1L & sc <= n_categories & dc >= 1L & dc <= n_categories
    sc <- sc[valid]
    dc <- dc[valid]
    mat <- matrix(
      tabulate((sc - 1L) * n_categories + dc, nbins = n_cats_sq) +
      tabulate((dc - 1L) * n_categories + sc, nbins = n_cats_sq),
      n_categories, n_categories
    )
    # Welford update: accumulate running mean and sum of squared deviations.
    delta <- mat - perm_mean
    perm_mean <- perm_mean + delta / k
    perm_M2   <- perm_M2 + delta * (mat - perm_mean)
  }

  sd_perm <- sqrt(pmax(0, perm_M2 / max(1L, n_perms - 1L)))
  z <- matrix(NA_real_, n_categories, n_categories)
  nz <- sd_perm > 0
  z[nz] <- (observed_counts[nz] - perm_mean[nz]) / sd_perm[nz]
  z
}

resolve_phase4_groupby_columns <- function(requested_groupby, color_data, max_categories = 50L) {
  if (length(color_data) == 0L) {
    return(character())
  }

  categorical_columns <- names(color_data)[!vapply(color_data, function(info) {
    isTRUE(info$meta$is_continuous)
  }, logical(1))]

  requested <- trimws(as.character(requested_groupby %||% character()))
  requested <- requested[nzchar(requested)]
  if (length(requested) == 0L) {
    return(character())
  }

  requested_lower <- tolower(requested)
  if ("none" %in% requested_lower) {
    return(character())
  }

  columns <- character()
  if ("auto" %in% requested_lower) {
    auto_columns <- categorical_columns[vapply(categorical_columns, function(column_name) {
      categories <- unlist(color_data[[column_name]]$meta$categories %||% list(), use.names = FALSE)
      length(categories) >= 2L && length(categories) <= max_categories
    }, logical(1))]
    columns <- c(columns, auto_columns)
  }

  explicit_columns <- requested[!(requested_lower %in% c("auto", "none"))]
  columns <- c(columns, explicit_columns)
  columns <- unique(intersect(columns, categorical_columns))
  columns
}

build_marker_genes <- function(
  expression,
  gene_names,
  color_data,
  requested_groupby = NULL,
  top_n = 20L,
  min_cells = 3L,
  marker_test = "mean_diff"
) {
  markers <- empty_named_list()
  groupby_columns <- resolve_phase4_groupby_columns(
    requested_groupby = requested_groupby,
    color_data = color_data
  )
  if (length(groupby_columns) == 0L || is.null(expression) || length(gene_names %||% character()) == 0L) {
    return(markers)
  }

  top_n <- max(1L, as.integer(top_n %||% 20L))
  min_cells <- max(2L, as.integer(min_cells %||% 3L))

  for (color_name in groupby_columns) {
    info <- color_data[[color_name]]
    categories <- as.character(unlist(info$meta$categories %||% list(), use.names = FALSE))
    codes <- as.integer(info$values) + 1L
    valid_codes <- is.finite(codes) & codes > 0L
    group_markers <- empty_named_list()

    for (category_idx in seq_along(categories)) {
      category_name <- categories[[category_idx]]
      if (identical(category_name, "(missing)")) {
        next
      }

      pos_idx <- which(codes == category_idx)
      neg_idx <- which(valid_codes & codes != category_idx)
      if (length(pos_idx) < min_cells || length(neg_idx) < min_cells) {
        next
      }

      ranked <- rank_marker_genes_between_groups(
        expression = expression,
        gene_names = gene_names,
        pos_idx = pos_idx,
        neg_idx = neg_idx,
        top_n = top_n,
        marker_test = marker_test
      )
      if (length(ranked$genes) == 0L) {
        next
      }

      group_markers[[category_name]] <- as_json_array(ranked$genes)
    }

    if (length(group_markers) > 0L) {
      markers[[color_name]] <- group_markers
    }
  }

  markers
}

# Compute contact-conditioned marker genes for all source→target category pairs.
#
# For each categorical color column and each (source, target) cell-type pair,
# cells of the source type are split into two groups based on whether they
# physically contact target cells in the neighbor graph:
#   - "contact" (pos_idx):     source cells with >= min_neighbors edges to a target cell
#   - "non-contact" (neg_idx): source cells with zero edges to any target cell
#
# Marker genes are then ranked between these two groups, identifying genes
# whose expression is modulated by proximity to the target cell type.
#
# Outer loop structure:
#   color_name  → categorical column (e.g. "cell_type")
#     source_name → source category (e.g. "Neuron")
#       target_name → target category ranked by enrichment z-score, then edge count
#
# Filtering:
#   - Skips "(missing)" categories
#   - Only considers targets that share at least one edge with the source
#   - Skips source→target pairs where either contact or non-contact group is
#     smaller than min_cells (result is stored with available=FALSE)
#   - Limits to top_targets targets per source (by descending z-score, then edge count)
build_interaction_markers <- function(
  expression,
  gene_names,
  color_data,
  requested_groupby = NULL,
  neighbor_stats,
  section_index_map,
  section_edges,
  top_targets = 8L,
  top_genes = 12L,
  min_cells = 30L,
  min_neighbors = 1L,
  marker_test = "mean_diff"
) {
  interactions <- empty_named_list()
  if (is.null(expression) || !has_neighbor_edges(section_edges) || length(gene_names %||% character()) == 0L) {
    return(interactions)
  }

  groupby_columns <- resolve_phase4_groupby_columns(
    requested_groupby = requested_groupby,
    color_data = color_data
  )
  if (length(groupby_columns) == 0L) {
    return(interactions)
  }

  top_targets <- max(1L, as.integer(top_targets %||% 8L))
  top_genes <- max(1L, as.integer(top_genes %||% 12L))
  min_cells <- max(2L, as.integer(min_cells %||% 30L))
  min_neighbors <- max(1L, as.integer(min_neighbors %||% 1L))

  for (color_name in groupby_columns) {
    stats_entry <- neighbor_stats[[color_name]] %||% NULL
    if (is.null(stats_entry)) {
      next
    }

    categories <- as.character(unlist(stats_entry$categories %||% list(), use.names = FALSE))
    if (length(categories) == 0L) {
      next
    }

    # counts_matrix[i, j]: total edges between category i and category j (symmetric)
    # zscore_matrix[i, j]: enrichment z-score from permutation testing (NA if unavailable)
    counts_matrix <- neighbor_counts_to_matrix(
      counts = stats_entry$counts,
      n_categories = length(categories)
    )
    zscore_matrix <- neighbor_counts_to_matrix(
      counts = stats_entry$zscore,
      n_categories = length(categories),
      default = NA_real_
    )

    info <- color_data[[color_name]]
    codes <- as.integer(info$values) + 1L
    # neighbor_category_counts[cell_i, cat_j]: how many neighbors of cell_i belong to category j
    neighbor_category_counts <- build_neighbor_category_counts(
      labels = codes,
      n_categories = length(categories),
      section_index_map = section_index_map,
      section_edges = section_edges
    )

    color_interactions <- empty_named_list()
    for (source_idx in seq_along(categories)) {
      source_name <- categories[[source_idx]]
      if (identical(source_name, "(missing)")) {
        next
      }

      # source_mask: logical vector identifying all cells of this source type
      source_mask <- codes == source_idx
      if (!any(source_mask, na.rm = TRUE)) {
        next
      }

      # Filter to targets that actually share edges with this source, excluding
      # self-interactions and missing labels. Rank by descending z-score (most
      # enriched first), then by raw edge count as a tie-breaker.
      row_counts <- counts_matrix[source_idx, ]
      candidate_targets <- which(row_counts > 0)
      candidate_targets <- setdiff(candidate_targets, source_idx)
      candidate_targets <- candidate_targets[categories[candidate_targets] != "(missing)"]
      if (length(candidate_targets) == 0L) {
        next
      }

      ranked_targets <- candidate_targets[order(
        ifelse(is.finite(zscore_matrix[source_idx, candidate_targets]), -zscore_matrix[source_idx, candidate_targets], Inf),
        -row_counts[candidate_targets],
        categories[candidate_targets]
      )]
      ranked_targets <- utils::head(ranked_targets, top_targets)

      source_result <- empty_named_list()
      for (target_idx in ranked_targets) {
        target_name <- categories[[target_idx]]
        # Per-cell count of how many neighbors belong to the target category.
        target_neighbor_counts <- neighbor_category_counts[, target_idx]
        # Contact group: source cells with at least min_neighbors target neighbors.
        # Non-contact group: source cells with zero target neighbors.
        # Cells between 1..(min_neighbors-1) neighbors are excluded from both groups.
        pos_idx <- which(source_mask & target_neighbor_counts >= min_neighbors)
        neg_idx <- which(source_mask & target_neighbor_counts == 0L)
        n_pos <- length(pos_idx)
        n_neg <- length(neg_idx)

        entry <- list(
          available = n_pos >= min_cells && n_neg >= min_cells,
          genes = list(),
          logfoldchanges = list(),
          pvals_adj = list(),
          n_contact = n_pos,
          n_non_contact = n_neg,
          pct_contact = 100 * n_pos / max(1L, n_pos + n_neg),
          mean_target_neighbors_contact = if (n_pos > 0L) mean(target_neighbor_counts[pos_idx]) else 0,
          mean_target_neighbors_non_contact = if (n_neg > 0L) mean(target_neighbor_counts[neg_idx]) else 0,
          target_edge_count = as.numeric(row_counts[[target_idx]] %||% 0),
          target_zscore = if (is.finite(zscore_matrix[source_idx, target_idx])) {
            as.numeric(zscore_matrix[source_idx, target_idx])
          } else {
            NULL
          }
        )

        if (!isTRUE(entry$available)) {
          entry$reason <- "insufficient_cells"
          entry$min_cells_required <- min_cells
          source_result[[target_name]] <- entry
          next
        }

        ranked <- rank_marker_genes_between_groups(
          expression = expression,
          gene_names = gene_names,
          pos_idx = pos_idx,
          neg_idx = neg_idx,
          top_n = top_genes,
          marker_test = marker_test
        )
        entry$genes <- as_json_array(ranked$genes)
        entry$logfoldchanges <- as_json_array(as.numeric(ranked$effect))
        entry$pvals_adj <- if (!is.null(ranked$pval_adj)) as_json_array(as.numeric(ranked$pval_adj)) else list()
        source_result[[target_name]] <- entry
      }

      if (length(source_result) > 0L) {
        color_interactions[[source_name]] <- source_result
      }
    }

    if (length(color_interactions) > 0L) {
      interactions[[color_name]] <- color_interactions
    }
  }

  interactions
}

rank_marker_genes_between_groups <- function(expression, gene_names, pos_idx, neg_idx, top_n = 20L, marker_test = "mean_diff") {
  if (identical(tolower(marker_test %||% "mean_diff"), "wilcoxon")) {
    return(rank_marker_genes_wilcoxon(
      expression = expression,
      gene_names = gene_names,
      pos_idx = pos_idx,
      neg_idx = neg_idx,
      top_n = top_n
    ))
  }

  top_n <- max(1L, as.integer(top_n %||% 20L))
  if (length(pos_idx) == 0L || length(neg_idx) == 0L) {
    return(list(genes = character(), effect = numeric()))
  }

  mean_pos <- matrix_rowmeans(expression[, pos_idx, drop = FALSE])
  mean_neg <- matrix_rowmeans(expression[, neg_idx, drop = FALSE])
  pct_pos <- matrix_detection_rate(expression[, pos_idx, drop = FALSE])
  pct_neg <- matrix_detection_rate(expression[, neg_idx, drop = FALSE])

  effect <- mean_pos - mean_neg
  pct_diff <- pct_pos - pct_neg
  valid <- is.finite(effect) &
    is.finite(pct_diff) &
    is.finite(mean_pos) &
    is.character(gene_names) &
    !is.na(gene_names) &
    nzchar(gene_names)

  preferred <- valid & (effect > 0 | (effect == 0 & pct_diff > 0))
  candidate_idx <- if (any(preferred)) {
    which(preferred)
  } else {
    which(valid & mean_pos > 0)
  }
  if (length(candidate_idx) == 0L) {
    return(list(genes = character(), effect = numeric()))
  }

  ord <- order(
    -effect[candidate_idx],
    -pct_diff[candidate_idx],
    -mean_pos[candidate_idx],
    gene_names[candidate_idx]
  )
  top_idx <- candidate_idx[utils::head(ord, top_n)]
  list(
    genes = as.character(gene_names[top_idx]),
    effect = as.numeric(effect[top_idx])
  )
}

rank_marker_genes_wilcoxon <- function(
  expression, gene_names, pos_idx, neg_idx, top_n = 20L, prefilter_n = NULL
) {
  prefilter_n <- max(top_n, as.integer(prefilter_n %||% (top_n * 5L)))
  if (!length(pos_idx) || !length(neg_idx)) {
    return(list(genes = character(), effect = numeric(), pval_adj = numeric()))
  }

  mean_pos <- matrix_rowmeans(expression[, pos_idx, drop = FALSE])
  mean_neg <- matrix_rowmeans(expression[, neg_idx, drop = FALSE])
  pct_pos  <- matrix_detection_rate(expression[, pos_idx, drop = FALSE])
  pct_neg  <- matrix_detection_rate(expression[, neg_idx, drop = FALSE])
  effect   <- mean_pos - mean_neg
  pct_diff <- pct_pos - pct_neg

  valid <- is.finite(effect) & is.character(gene_names) & !is.na(gene_names) & nzchar(gene_names)
  preferred <- valid & (effect > 0 | (effect == 0 & pct_diff > 0))
  cands <- if (any(preferred)) which(preferred) else which(valid & mean_pos > 0)
  if (!length(cands)) {
    return(list(genes = character(), effect = numeric(), pval_adj = numeric()))
  }

  cands <- cands[utils::head(order(-effect[cands], -pct_diff[cands]), prefilter_n)]

  pos_m <- expression[cands, pos_idx, drop = FALSE]
  neg_m <- expression[cands, neg_idx, drop = FALSE]
  pvals <- vapply(seq_along(cands), function(i) {
    x <- as.numeric(pos_m[i, ])
    y <- as.numeric(neg_m[i, ])
    if (stats::var(c(x, y)) == 0) return(1)
    tryCatch(
      stats::wilcox.test(x, y, alternative = "greater", exact = FALSE)$p.value,
      error = function(e) 1
    )
  }, numeric(1))

  padj <- stats::p.adjust(pvals, method = "BH")
  ord  <- order(padj, -effect[cands])
  top  <- cands[utils::head(ord, top_n)]
  padj_top <- padj[utils::head(ord, top_n)]

  list(
    genes    = as.character(gene_names[top]),
    effect   = as.numeric(effect[top]),
    pval_adj = as.numeric(padj_top)
  )
}

matrix_detection_rate <- function(x) {
  if (is.null(x) || ncol(x) == 0L) {
    return(rep.int(0, nrow(x %||% matrix(numeric(), nrow = 0L))))
  }
  if (methods::is(x, "Matrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Sparse marker computation requires the Matrix package.")
    }
    return(as.numeric(Matrix::rowMeans(x > 0)))
  }
  as.numeric(rowMeans(x > 0))
}

neighbor_counts_to_matrix <- function(counts, n_categories, default = 0) {
  if (is.null(counts) || length(counts) == 0L || n_categories < 1L) {
    return(matrix(default, nrow = n_categories, ncol = n_categories))
  }

  out <- matrix(default, nrow = n_categories, ncol = n_categories)
  limit <- min(length(counts), n_categories)
  for (row_idx in seq_len(limit)) {
    row_values <- as.numeric(unlist(counts[[row_idx]] %||% list(), use.names = FALSE))
    if (length(row_values) == 0L) {
      next
    }
    row_limit <- min(length(row_values), n_categories)
    out[row_idx, seq_len(row_limit)] <- row_values[seq_len(row_limit)]
  }
  out
}

build_neighbor_category_counts <- function(labels, n_categories, section_index_map, section_edges) {
  n_cells <- length(labels)
  counts <- matrix(0L, nrow = n_cells, ncol = n_categories)
  if (n_cells == 0L || n_categories == 0L || !has_neighbor_edges(section_edges)) {
    return(counts)
  }

  labels <- as.integer(labels)
  for (section_id in names(section_index_map)) {
    idx <- section_index_map[[section_id]]
    edge_values <- as.integer(section_edges[[section_id]] %||% integer())
    if (length(idx) == 0L || length(edge_values) == 0L) {
      next
    }

    edge_matrix <- matrix(edge_values, ncol = 2L, byrow = TRUE) + 1L
    src_local <- edge_matrix[, 1]
    dst_local <- edge_matrix[, 2]
    src_global <- idx[src_local]
    dst_global <- idx[dst_local]
    src_labels <- labels[src_global]
    dst_labels <- labels[dst_global]

    valid_dst <- is.finite(dst_labels) & dst_labels >= 1L & dst_labels <= n_categories
    if (any(valid_dst)) {
      for (pair_idx in which(valid_dst)) {
        counts[src_global[[pair_idx]], dst_labels[[pair_idx]]] <- counts[src_global[[pair_idx]], dst_labels[[pair_idx]]] + 1L
      }
    }

    valid_src <- is.finite(src_labels) & src_labels >= 1L & src_labels <= n_categories
    if (any(valid_src)) {
      for (pair_idx in which(valid_src)) {
        counts[dst_global[[pair_idx]], src_labels[[pair_idx]]] <- counts[dst_global[[pair_idx]], src_labels[[pair_idx]]] + 1L
      }
    }
  }

  counts
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

  missing_label <- "(missing)"
  factor_column <- if (is.factor(column)) column else factor(as.character(column))
  if (any(is.na(factor_column))) {
    factor_column <- addNA(factor_column, ifany = TRUE)
    levels(factor_column)[is.na(levels(factor_column))] <- missing_label
  }
  values <- as.numeric(factor_column) - 1
  meta <- list(
    is_continuous = FALSE,
    categories = as_json_array(as.character(levels(factor_column))),
    vmin = 0,
    vmax = max(0, length(levels(factor_column)) - 1)
  )

  list(values = values, meta = meta)
}

build_gene_data <- function(expression, gene_names, selected_genes, sparse_zero_threshold = 0.8) {
  if (is.null(expression) || length(selected_genes) == 0) {
    return(list(
      values = empty_named_list(),
      meta = empty_named_list(),
      encodings = empty_named_list()
    ))
  }

  gene_names <- as.character(gene_names)
  gene_index <- match(selected_genes, gene_names)

  values <- vector("list", length(selected_genes))
  meta <- vector("list", length(selected_genes))
  encodings <- vector("list", length(selected_genes))
  names(values) <- selected_genes
  names(meta) <- selected_genes
  names(encodings) <- selected_genes

  for (i in seq_along(selected_genes)) {
    gene <- selected_genes[[i]]
    idx <- gene_index[[i]]
    vector <- extract_gene_vector(expression, idx)
    finite <- vector[is.finite(vector)]
    zero_fraction <- if (length(finite) == 0L) 1 else mean(finite == 0)

    values[[gene]] <- vector
    meta[[gene]] <- list(
      vmin = if (length(finite) > 0) min(finite) else 0,
      vmax = if (length(finite) > 0) max(finite) else 1
    )
    encodings[[gene]] <- if (isTRUE(zero_fraction >= sparse_zero_threshold)) "sparse" else "dense"
  }

  list(values = values, meta = meta, encodings = encodings)
}

extract_gene_vector <- function(expression, idx) {
  if (methods::is(expression, "Matrix")) {
    return(as.numeric(expression[idx, , drop = TRUE]))
  }
  as.numeric(expression[idx, , drop = TRUE])
}

encode_sparse_gene_values <- function(values, pack_sparse = TRUE, pack_min_nnz = 32L) {
  finite <- is.finite(values)
  nonzero <- finite & (values != 0)
  if (!any(nonzero)) {
    nan_idx <- which(is.na(values)) - 1L
    out <- list(i = list(), v = list())
    if (length(nan_idx) > 0L) {
      out$nan <- as_json_array(as.integer(nan_idx))
    }
    return(out)
  }

  nz_idx <- which(nonzero) - 1L
  nz_vals <- as.numeric(values[nonzero])
  should_pack_sparse <- isTRUE(pack_sparse) && length(nz_idx) >= as.integer(pack_min_nnz %||% 32L)
  out <- if (should_pack_sparse) {
    list(
      ib64 = pack_uint32_base64(as.integer(nz_idx)),
      vb64 = pack_float32_base64(nz_vals)
    )
  } else {
    list(
      i = as_json_array(as.integer(nz_idx)),
      v = as_json_array(nz_vals)
    )
  }

  nan_idx <- which(is.na(values)) - 1L
  if (length(nan_idx) > 0L) {
    out$nan <- as_json_array(as.integer(nan_idx))
  }
  out
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
