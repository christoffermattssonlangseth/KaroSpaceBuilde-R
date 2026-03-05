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
    '<div class="info-text">Built from an RDS input via the standalone KaroSpaceBuildR exporter.</div>',
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
  assay = NULL,
  metadata_input = NULL,
  metadata_input_columns = NULL,
  metadata_prefix = NULL,
  metadata_columns = NULL,
  outline_by = NULL,
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
    assay = assay,
    metadata_input = metadata_input,
    metadata_input_columns = metadata_input_columns,
    metadata_prefix = metadata_prefix,
    metadata_columns = metadata_columns,
    outline_by = outline_by
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
