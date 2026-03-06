karospace_default_palette <- function() {
  c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    "#393b79", "#5254a3", "#6b6ecf", "#9c9ede", "#637939",
    "#8ca252", "#b5cf6b", "#cedb9c", "#8c6d31", "#bd9e39",
    "#e7ba52", "#e7cb94", "#843c39", "#ad494a", "#d6616b",
    "#e7969c", "#7b4173", "#a55194", "#ce6dbd", "#de9ed6"
  )
}

karospace_default_metadata_labels <- function() {
  list(
    last_score = "disease score",
    last_day = "day of sacrifice"
  )
}

karospace_default_viewer_info <- function() {
  paste0(
    '<div class="info-block">',
    '<div class="info-title">Viewer</div>',
    '<div class="info-text">KaroSpace interactive spatial viewer for exploring ',
    'sections, cell types, and gene expression.</div>',
    '</div>',
    '<div class="info-block">',
    '<div class="info-title">Source</div>',
    '<div class="info-text">Built from an RDS input via the standalone KaroSpaceR exporter.</div>',
    '</div>'
  )
}

karospace_theme_context <- function(theme = "light") {
  theme <- tolower(theme %||% "light")
  if (identical(theme, "dark")) {
    return(list(
      background = "#1a1a1a",
      text_color = "#e0e0e0",
      header_bg = "#2a2a2a",
      panel_bg = "#2a2a2a",
      border_color = "#404040",
      input_bg = "#333333",
      muted_color = "#888888",
      hover_bg = "#3a3a3a",
      graph_color = "rgba(255, 255, 255, 0.12)",
      theme_icon = "\u2600\ufe0f",
      initial_theme = "dark"
    ))
  }

  list(
    background = "#f5f5f5",
    text_color = "#1a1a1a",
    header_bg = "#ffffff",
    panel_bg = "#ffffff",
    border_color = "#e0e0e0",
    input_bg = "#ffffff",
    muted_color = "#666666",
    hover_bg = "#f0f0f0",
    graph_color = "rgba(0, 0, 0, 0.12)",
    theme_icon = "\ud83c\udf19",
    initial_theme = "light"
  )
}

render_viewer_html <- function(
  payload,
  output_path,
  title = "KaroSpace",
  theme = "light",
  min_panel_size = 150,
  spot_size = 2,
  outline_by = NULL,
  viewer_info_html = NULL,
  viewer_shell_path = NULL
) {
  shell_path <- resolve_asset_path(
    rel_path = file.path("inst", "viewer", "karospace_viewer_shell.html"),
    explicit_path = viewer_shell_path
  )
  html <- paste(readLines(shell_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

  theme_ctx <- karospace_theme_context(theme)
  viewer_info_html <- viewer_info_html %||% karospace_default_viewer_info()
  payload_json <- gsub("</", "<\\\\/", compact_json(payload), fixed = TRUE)

  replacements <- c(
    "__KAROSPACE_TITLE__" = escape_html(title),
    "__KAROSPACE_MIN_PANEL_SIZE__" = as.character(as.integer(min_panel_size)),
    "__KAROSPACE_MAX_PANEL_SIZE__" = as.character(as.integer(min_panel_size * 2L)),
    "__KAROSPACE_SPOT_SIZE__" = format(as.numeric(spot_size), scientific = FALSE, trim = TRUE),
    "__KAROSPACE_DATA_JSON__" = payload_json,
    "__KAROSPACE_PALETTE_JSON__" = compact_json(unname(karospace_default_palette())),
    "__KAROSPACE_METADATA_LABELS_JSON__" = compact_json(karospace_default_metadata_labels()),
    "__KAROSPACE_OUTLINE_BY_JSON__" = compact_json(outline_by),
    "__KAROSPACE_VIEWER_INFO_HTML_JSON__" = compact_json(viewer_info_html),
    "__KAROSPACE_VIEWER_INFO_HTML__" = viewer_info_html,
    "__KAROSPACE_THEME_ICON__" = theme_ctx$theme_icon,
    "__KAROSPACE_INITIAL_THEME__" = theme_ctx$initial_theme,
    "__KAROSPACE_FAVICON_LINK__" = "",
    "__KAROSPACE_FOOTER_LOGO__" = paste0(
      '<div class="footer-logo">',
      '<span>KaroSpace</span>',
      '<span class="footer-link">Standalone R Export</span>',
      '</div>'
    ),
    "__KAROSPACE_BACKGROUND__" = theme_ctx$background,
    "__KAROSPACE_TEXT_COLOR__" = theme_ctx$text_color,
    "__KAROSPACE_HEADER_BG__" = theme_ctx$header_bg,
    "__KAROSPACE_PANEL_BG__" = theme_ctx$panel_bg,
    "__KAROSPACE_BORDER_COLOR__" = theme_ctx$border_color,
    "__KAROSPACE_INPUT_BG__" = theme_ctx$input_bg,
    "__KAROSPACE_MUTED_COLOR__" = theme_ctx$muted_color,
    "__KAROSPACE_HOVER_BG__" = theme_ctx$hover_bg,
    "__KAROSPACE_GRAPH_COLOR__" = theme_ctx$graph_color
  )

  for (token in names(replacements)) {
    html <- gsub(token, replacements[[token]], html, fixed = TRUE)
  }

  output_path <- resolve_output_path(output_path)
  writeLines(html, con = output_path, useBytes = TRUE)
  invisible(output_path)
}

export_karospace_viewer <- function(
  input,
  output_path,
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
  neighbor_stats_seed = 42L,
  title = "KaroSpace",
  theme = "light",
  min_panel_size = 150,
  spot_size = 2,
  viewer_info_html = NULL,
  viewer_shell_path = NULL
) {
  payload <- build_viewer_payload(
    input = input,
    groupby = groupby,
    initial_color = initial_color,
    additional_colors = additional_colors,
    genes = genes,
    top_genes_n = top_genes_n,
    assay = assay,
    metadata_input = metadata_input,
    metadata_input_columns = metadata_input_columns,
    metadata_prefix = metadata_prefix,
    neighbor_mode = neighbor_mode,
    neighbor_graph = neighbor_graph,
    neighbor_k = neighbor_k,
    metadata_columns = metadata_columns,
    outline_by = outline_by,
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

  render_viewer_html(
    payload = payload,
    output_path = output_path,
    title = title,
    theme = theme,
    min_panel_size = min_panel_size,
    spot_size = spot_size,
    outline_by = outline_by %||% payload$outline_by,
    viewer_info_html = viewer_info_html,
    viewer_shell_path = viewer_shell_path
  )
}

karospace_build_config_keys <- function() {
  c(
    "version",
    "input",
    "output",
    "groupby",
    "initial_color",
    "additional_colors",
    "genes",
    "top_genes_n",
    "assay",
    "metadata_input",
    "metadata_input_columns",
    "metadata_prefix",
    "neighbor_mode",
    "neighbor_graph",
    "neighbor_k",
    "metadata_columns",
    "outline_by",
    "lightweight",
    "marker_genes_groupby",
    "marker_genes_top_n",
    "interaction_markers_groupby",
    "interaction_markers_top_targets",
    "interaction_markers_top_genes",
    "interaction_markers_min_cells",
    "interaction_markers_min_neighbors",
    "marker_test",
    "neighbor_stats_permutations",
    "neighbor_stats_seed",
    "title",
    "theme",
    "min_panel_size",
    "spot_size"
  )
}

karospace_config_scalar_character <- function(x, field, allow_null = TRUE) {
  if (is.null(x) || length(x) == 0L) {
    if (allow_null) {
      return(NULL)
    }
    stop("Missing required config field: ", field)
  }
  value <- as.character(x)[[1L]]
  if (!nzchar(value)) {
    if (allow_null) {
      return(NULL)
    }
    stop("Config field must be a non-empty string: ", field)
  }
  value
}

karospace_config_character_vector <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NULL)
  }
  values <- as.character(x)
  values <- values[nzchar(values)]
  if (length(values) == 0L) {
    return(NULL)
  }
  unname(values)
}

karospace_config_scalar_integer <- function(x, field, allow_null = TRUE, min_value = NULL) {
  if (is.null(x) || length(x) == 0L) {
    if (allow_null) {
      return(NULL)
    }
    stop("Missing required config field: ", field)
  }
  value <- suppressWarnings(as.integer(x[[1L]]))
  if (is.na(value)) {
    stop("Config field must be an integer: ", field)
  }
  if (!is.null(min_value) && value < min_value) {
    stop("Config field must be >= ", min_value, ": ", field)
  }
  value
}

karospace_config_scalar_numeric <- function(x, field, allow_null = TRUE, min_value = NULL) {
  if (is.null(x) || length(x) == 0L) {
    if (allow_null) {
      return(NULL)
    }
    stop("Missing required config field: ", field)
  }
  value <- suppressWarnings(as.numeric(x[[1L]]))
  if (is.na(value)) {
    stop("Config field must be numeric: ", field)
  }
  if (!is.null(min_value) && value < min_value) {
    stop("Config field must be >= ", min_value, ": ", field)
  }
  value
}

karospace_config_scalar_logical <- function(x, field, allow_null = TRUE) {
  if (is.null(x) || length(x) == 0L) {
    if (allow_null) {
      return(NULL)
    }
    stop("Missing required config field: ", field)
  }

  value <- x[[1L]]
  if (is.logical(value)) {
    return(isTRUE(value))
  }

  normalized <- tolower(trimws(as.character(value)))
  if (normalized %in% c("true", "t", "1", "yes", "y", "on")) {
    return(TRUE)
  }
  if (normalized %in% c("false", "f", "0", "no", "n", "off")) {
    return(FALSE)
  }

  stop("Config field must be boolean: ", field)
}

resolve_config_path_value <- function(path_value, base_dir = NULL, must_exist = FALSE) {
  if (is.null(path_value) || !nzchar(path_value)) {
    return(NULL)
  }

  resolved <- path_value
  if (!grepl("^(/|[A-Za-z]:[/\\\\]|~)", resolved) && !is.null(base_dir)) {
    resolved <- file.path(base_dir, resolved)
  }

  if (must_exist) {
    return(normalizePath(resolved, mustWork = TRUE))
  }

  normalizePath(resolved, mustWork = FALSE)
}

normalize_karospace_build_config <- function(
  config,
  config_path = NULL,
  require_required = FALSE,
  resolve_paths = TRUE
) {
  if (is.character(config) && length(config) == 1L && file.exists(config)) {
    return(read_karospace_build_config(config))
  }
  if (!is.list(config)) {
    stop("Build config must be a named list or a JSON config path.")
  }

  config_names <- names(config) %||% character()
  if (length(config_names) == 0L) {
    stop("Build config must be a named list.")
  }

  unknown_fields <- setdiff(config_names, karospace_build_config_keys())
  if (length(unknown_fields) > 0L) {
    stop("Unsupported build config fields: ", paste(sort(unknown_fields), collapse = ", "))
  }

  base_dir <- NULL
  if (!is.null(config_path)) {
    base_dir <- dirname(normalizePath(config_path, mustWork = FALSE))
  }

  normalized <- list(
    version = karospace_config_scalar_integer(config$version %||% 1L, "version", allow_null = FALSE, min_value = 1L),
    input = karospace_config_scalar_character(config$input, "input", allow_null = !require_required),
    output = karospace_config_scalar_character(config$output, "output", allow_null = !require_required),
    groupby = karospace_config_scalar_character(config$groupby, "groupby", allow_null = !require_required),
    initial_color = karospace_config_scalar_character(config$initial_color, "initial_color", allow_null = !require_required),
    additional_colors = karospace_config_character_vector(config$additional_colors),
    genes = karospace_config_character_vector(config$genes),
    top_genes_n = karospace_config_scalar_integer(config$top_genes_n, "top_genes_n", allow_null = TRUE, min_value = 1L),
    assay = karospace_config_scalar_character(config$assay, "assay", allow_null = TRUE),
    metadata_input = karospace_config_scalar_character(config$metadata_input, "metadata_input", allow_null = TRUE),
    metadata_input_columns = karospace_config_character_vector(config$metadata_input_columns),
    metadata_prefix = karospace_config_scalar_character(config$metadata_prefix, "metadata_prefix", allow_null = TRUE),
    neighbor_mode = karospace_config_scalar_character(config$neighbor_mode %||% "spatial", "neighbor_mode", allow_null = FALSE),
    neighbor_graph = karospace_config_scalar_character(config$neighbor_graph, "neighbor_graph", allow_null = TRUE),
    neighbor_k = karospace_config_scalar_integer(config$neighbor_k %||% 6L, "neighbor_k", allow_null = FALSE, min_value = 1L),
    metadata_columns = karospace_config_character_vector(config$metadata_columns),
    outline_by = karospace_config_scalar_character(config$outline_by, "outline_by", allow_null = TRUE),
    lightweight = karospace_config_scalar_logical(config$lightweight %||% FALSE, "lightweight", allow_null = FALSE),
    marker_genes_groupby = karospace_config_character_vector(config$marker_genes_groupby),
    marker_genes_top_n = karospace_config_scalar_integer(config$marker_genes_top_n %||% 20L, "marker_genes_top_n", allow_null = FALSE, min_value = 1L),
    interaction_markers_groupby = karospace_config_character_vector(config$interaction_markers_groupby),
    interaction_markers_top_targets = karospace_config_scalar_integer(config$interaction_markers_top_targets %||% 8L, "interaction_markers_top_targets", allow_null = FALSE, min_value = 1L),
    interaction_markers_top_genes = karospace_config_scalar_integer(config$interaction_markers_top_genes %||% 12L, "interaction_markers_top_genes", allow_null = FALSE, min_value = 1L),
    interaction_markers_min_cells = karospace_config_scalar_integer(config$interaction_markers_min_cells %||% 30L, "interaction_markers_min_cells", allow_null = FALSE, min_value = 2L),
    interaction_markers_min_neighbors = karospace_config_scalar_integer(config$interaction_markers_min_neighbors %||% 1L, "interaction_markers_min_neighbors", allow_null = FALSE, min_value = 1L),
    marker_test = tolower(karospace_config_scalar_character(config$marker_test %||% "mean_diff", "marker_test", allow_null = FALSE)),
    neighbor_stats_permutations = karospace_config_scalar_integer(config$neighbor_stats_permutations %||% 0L, "neighbor_stats_permutations", allow_null = FALSE, min_value = 0L),
    neighbor_stats_seed = karospace_config_scalar_integer(config$neighbor_stats_seed %||% 42L, "neighbor_stats_seed", allow_null = FALSE),
    title = karospace_config_scalar_character(config$title, "title", allow_null = TRUE),
    theme = karospace_config_scalar_character(config$theme %||% "light", "theme", allow_null = FALSE),
    min_panel_size = karospace_config_scalar_numeric(config$min_panel_size %||% 150, "min_panel_size", allow_null = FALSE, min_value = 1),
    spot_size = karospace_config_scalar_numeric(config$spot_size %||% 2, "spot_size", allow_null = FALSE, min_value = 0)
  )

  if (normalized$version != 1L) {
    stop("Unsupported build config version: ", normalized$version, ". Expected version 1.")
  }
  if (!normalized$neighbor_mode %in% KAROSPACE_NEIGHBOR_MODES) {
    stop("Unsupported neighbor_mode in build config: ", normalized$neighbor_mode)
  }
  if (!normalized$marker_test %in% KAROSPACE_MARKER_TESTS) {
    stop("Unsupported marker_test in build config: ", normalized$marker_test)
  }
  if (!normalized$theme %in% KAROSPACE_THEMES) {
    stop("Unsupported theme in build config: ", normalized$theme)
  }
  if (!is.null(normalized$genes) && !is.null(normalized$top_genes_n)) {
    warning("top_genes_n is ignored because genes was provided explicitly.", call. = FALSE)
    normalized$top_genes_n <- NULL
  }

  if (isTRUE(resolve_paths)) {
    normalized$input <- resolve_config_path_value(normalized$input, base_dir = base_dir, must_exist = !is.null(normalized$input))
    normalized$output <- resolve_config_path_value(normalized$output, base_dir = base_dir, must_exist = FALSE)
    normalized$metadata_input <- resolve_config_path_value(
      normalized$metadata_input,
      base_dir = base_dir,
      must_exist = !is.null(normalized$metadata_input)
    )
  }

  normalized
}

read_karospace_build_config <- function(path) {
  config_path <- normalizePath(path, mustWork = TRUE)
  config <- jsonlite::fromJSON(config_path, simplifyVector = TRUE)
  normalize_karospace_build_config(
    config = config,
    config_path = config_path,
    require_required = FALSE,
    resolve_paths = TRUE
  )
}

write_karospace_build_config <- function(config, path) {
  config_path <- resolve_output_path(path)
  normalized <- normalize_karospace_build_config(
    config = config,
    require_required = FALSE,
    resolve_paths = FALSE
  )
  json_ready <- karospace_build_config_for_json(normalized)
  json <- jsonlite::toJSON(
    json_ready,
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    digits = NA,
    pretty = TRUE
  )
  writeLines(json, con = config_path, useBytes = TRUE)
  invisible(config_path)
}

karospace_build_config_for_json <- function(config) {
  normalized <- normalize_karospace_build_config(
    config = config,
    require_required = FALSE,
    resolve_paths = FALSE
  )
  array_fields <- c(
    "additional_colors",
    "genes",
    "metadata_input_columns",
    "metadata_columns",
    "marker_genes_groupby",
    "interaction_markers_groupby"
  )
  json_ready <- normalized
  for (field in array_fields) {
    if (!is.null(json_ready[[field]])) {
      json_ready[[field]] <- unname(as.list(json_ready[[field]]))
    }
  }
  json_ready
}

export_karospace_viewer_from_config <- function(config, config_path = NULL) {
  normalized <- if (is.character(config) && length(config) == 1L && file.exists(config)) {
    read_karospace_build_config(config)
  } else {
    normalize_karospace_build_config(
      config = config,
      config_path = config_path,
      require_required = TRUE,
      resolve_paths = TRUE
    )
  }

  export_karospace_viewer(
    input = normalized$input,
    output_path = normalized$output,
    groupby = normalized$groupby,
    initial_color = normalized$initial_color,
    additional_colors = normalized$additional_colors,
    genes = normalized$genes,
    top_genes_n = normalized$top_genes_n,
    assay = normalized$assay,
    metadata_input = normalized$metadata_input,
    metadata_input_columns = normalized$metadata_input_columns,
    metadata_prefix = normalized$metadata_prefix,
    neighbor_mode = normalized$neighbor_mode,
    neighbor_graph = normalized$neighbor_graph,
    neighbor_k = normalized$neighbor_k,
    metadata_columns = normalized$metadata_columns,
    outline_by = normalized$outline_by,
    lightweight = normalized$lightweight,
    marker_genes_groupby = normalized$marker_genes_groupby,
    marker_genes_top_n = normalized$marker_genes_top_n,
    interaction_markers_groupby = normalized$interaction_markers_groupby,
    interaction_markers_top_targets = normalized$interaction_markers_top_targets,
    interaction_markers_top_genes = normalized$interaction_markers_top_genes,
    interaction_markers_min_cells = normalized$interaction_markers_min_cells,
    interaction_markers_min_neighbors = normalized$interaction_markers_min_neighbors,
    marker_test = normalized$marker_test,
    neighbor_stats_permutations = normalized$neighbor_stats_permutations,
    neighbor_stats_seed = normalized$neighbor_stats_seed,
    title = normalized$title %||% "KaroSpace",
    theme = normalized$theme,
    min_panel_size = normalized$min_panel_size,
    spot_size = normalized$spot_size
  )

  invisible(normalized$output)
}
