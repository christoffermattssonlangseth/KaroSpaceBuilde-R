source("R/helpers.R")
source("R/source.R")
source("R/payload.R")
source("R/export.R")
source("R/inspect.R")

set.seed(7)
n_cells <- 60
n_genes <- 12

obs <- data.frame(
  sample_id = rep(c("section_a", "section_b", "section_c"), each = 20),
  cell_type = c(
    rep("A", 10), rep("B", 10),
    rep("A", 10), rep("C", 10),
    rep("B", 10), rep("C", 10)
  ),
  course = rep(c("naive", "peak", "recovery"), each = 20),
  stringsAsFactors = FALSE
)
rownames(obs) <- sprintf("cell_%03d", seq_len(n_cells))

coords <- cbind(
  x = stats::rnorm(n_cells),
  y = stats::rnorm(n_cells)
)

umap <- cbind(
  x = stats::rnorm(n_cells),
  y = stats::rnorm(n_cells)
)

expr <- matrix(
  rexp(n_cells * n_genes, rate = 1),
  nrow = n_genes,
  ncol = n_cells,
  dimnames = list(sprintf("Gene%02d", seq_len(n_genes)), NULL)
)
expr["Gene01", obs$cell_type == "A"] <- expr["Gene01", obs$cell_type == "A"] + 8
expr["Gene02", obs$cell_type == "B"] <- expr["Gene02", obs$cell_type == "B"] + 8
expr["Gene03", obs$cell_type == "C"] <- expr["Gene03", obs$cell_type == "C"] + 8
expr["Gene04", obs$sample_id == "section_a" & obs$cell_type == "A"] <-
  expr["Gene04", obs$sample_id == "section_a" & obs$cell_type == "A"] + 10
expr["Gene05", obs$sample_id == "section_b" & obs$cell_type == "A"] <-
  expr["Gene05", obs$sample_id == "section_b" & obs$cell_type == "A"] + 10
expr["Gene06", ] <- 0
expr["Gene06", c(1, 25)] <- c(4, 6)

section_edges <- list(
  section_a = as.integer(c(
    0, 10,
    1, 11,
    2, 12,
    3, 13,
    4, 14
  )),
  section_b = as.integer(c(
    0, 10,
    1, 11,
    2, 12
  )),
  section_c = as.integer(c(
    0, 10,
    1, 11,
    2, 12
  ))
)

toy <- list(
  obs = obs,
  coordinates = coords,
  umap = umap,
  expression = expr,
  neighbor_edges_by_section = section_edges,
  neighbors_key = "toy_neighbors"
)

cluster_meta <- data.frame(
  seurat_clusters = rep(c("0", "1", "2", NA), length.out = n_cells),
  stringsAsFactors = FALSE,
  row.names = rownames(obs)
)
merged_toy <- merge_external_metadata(
  primary = toy,
  secondary = cluster_meta[seq_len(n_cells / 2), , drop = FALSE],
  columns = c("seurat_clusters")
)
stopifnot("seurat_clusters" %in% names(merged_toy$obs))
stopifnot(identical(merged_toy$obs$seurat_clusters[1:3], c("0", "1", "2")))
merge_report <- extract_metadata_merge_report(merged_toy)
stopifnot(merge_report$overlap_rows == n_cells / 2)
stopifnot(merge_report$primary_rows == n_cells)
stopifnot(identical(merge_report$added_columns, c("seurat_clusters")))
reported_coverage <- merge_report$column_coverage$coverage[
  merge_report$column_coverage$column == "seurat_clusters"
]
expected_coverage <- count_non_missing_values(merged_toy$obs$seurat_clusters) / n_cells
stopifnot(isTRUE(all.equal(reported_coverage, expected_coverage)))

payload <- build_viewer_payload(
  input = merged_toy,
  groupby = "sample_id",
  initial_color = "cell_type",
  additional_colors = c("course", "seurat_clusters"),
  genes = c("Gene01", "Gene02"),
  neighbor_mode = "existing",
  marker_genes_groupby = "auto",
  marker_genes_top_n = 4L,
  interaction_markers_groupby = "cell_type",
  interaction_markers_top_targets = 2L,
  interaction_markers_top_genes = 3L,
  interaction_markers_min_cells = 3L,
  interaction_markers_min_neighbors = 1L,
  metadata_columns = c("course")
)

stopifnot(payload$n_sections == 3L)
stopifnot(payload$total_cells == n_cells)
stopifnot(identical(unlist(payload$available_colors, use.names = FALSE), c("cell_type", "course", "seurat_clusters")))
stopifnot(identical(unlist(payload$available_genes, use.names = FALSE), c("Gene01", "Gene02")))
stopifnot(isTRUE(payload$has_umap))
stopifnot(isTRUE(payload$has_neighbors))
stopifnot(identical(payload$neighbors_key, "toy_neighbors"))
stopifnot(identical(names(payload$genes_meta), c("Gene01", "Gene02")))
stopifnot(length(payload$sections[[1]]$genes$Gene01) == 20L)
stopifnot(identical(payload$gene_encodings$Gene01, "dense"))
stopifnot("(missing)" %in% unlist(payload$colors_meta$seurat_clusters$categories, use.names = FALSE))
stopifnot(is.character(payload$sections[[1]]$edges_b64) && nzchar(payload$sections[[1]]$edges_b64))
stopifnot(is.null(payload$sections[[1]]$edges))
stopifnot("cell_type" %in% names(payload$neighbor_stats))
stopifnot(length(payload$neighbor_stats$cell_type$categories) == 3L)
stopifnot(length(payload$neighbor_stats$cell_type$counts) == 3L)
stopifnot("cell_type" %in% names(payload$marker_genes))
stopifnot("Gene01" %in% unlist(payload$marker_genes$cell_type$A, use.names = FALSE))
stopifnot("Gene02" %in% unlist(payload$marker_genes$cell_type$B, use.names = FALSE))
stopifnot("Gene03" %in% unlist(payload$marker_genes$cell_type$C, use.names = FALSE))
stopifnot("cell_type" %in% names(payload$interaction_markers))
stopifnot(isTRUE(payload$interaction_markers$cell_type$A$B$available))
stopifnot("Gene04" %in% unlist(payload$interaction_markers$cell_type$A$B$genes, use.names = FALSE))
stopifnot(payload$interaction_markers$cell_type$A$B$n_contact >= 5L)
stopifnot(payload$interaction_markers$cell_type$A$B$n_non_contact >= 5L)

packed_payload <- build_viewer_payload(
  input = merged_toy,
  groupby = "sample_id",
  initial_color = "cell_type",
  additional_colors = c("course"),
  genes = c("Gene06"),
  neighbor_mode = "existing",
  pack_arrays = TRUE,
  pack_arrays_min_len = 1L,
  gene_sparse_zero_threshold = 0.9,
  gene_sparse_pack = TRUE,
  gene_sparse_pack_min_nnz = 1L,
  lightweight = TRUE
)
stopifnot(is.null(packed_payload$sections[[1]]$x))
stopifnot(is.character(packed_payload$sections[[1]]$xb64) && nzchar(packed_payload$sections[[1]]$xb64))
stopifnot(is.character(packed_payload$sections[[1]]$colors_b64$cell_type) && nzchar(packed_payload$sections[[1]]$colors_b64$cell_type))
stopifnot(identical(packed_payload$gene_encodings$Gene06, "sparse"))
stopifnot(is.character(packed_payload$sections[[1]]$genes_sparse$Gene06$ib64) && nzchar(packed_payload$sections[[1]]$genes_sparse$Gene06$ib64))
stopifnot(is.character(packed_payload$sections[[1]]$genes_sparse$Gene06$vb64) && nzchar(packed_payload$sections[[1]]$genes_sparse$Gene06$vb64))
stopifnot(length(packed_payload$marker_genes) == 0L)
stopifnot(length(packed_payload$interaction_markers) == 0L)

gene_info <- resolve_input_gene_names(toy)
stopifnot(identical(unlist(gene_info$gene_names[1:3], use.names = FALSE), c("Gene01", "Gene02", "Gene03")))
top_gene_selection <- resolve_selected_genes(
  genes = NULL,
  gene_names = rownames(expr),
  expression = expr,
  top_genes_n = 3L
)
stopifnot(length(top_gene_selection) == 3L)
stopifnot(setequal(top_gene_selection, c("Gene01", "Gene02", "Gene03")))

fake_staffli <- structure(
  list(),
  meta_data = data.frame(
    barcode = c("cell_b", "cell_a"),
    pxl_col_in_fullres = c(11, 22),
    pxl_row_in_fullres = c(33, 44),
    sampleID = c("1", "1"),
    stringsAsFactors = FALSE
  ),
  imgs = c(
    "/tmp/demo_section_a/spatial/tissue_hires_image.png"
  ),
  image_info = data.frame(
    sampleID = "1",
    stringsAsFactors = FALSE
  ),
  class = structure("Staffli", package = "semla")
)
staffli_coords <- extract_staffli_coordinate_df(fake_staffli)
stopifnot(identical(rownames(staffli_coords), c("cell_b", "cell_a")))
stopifnot(identical(unname(as.numeric(staffli_coords["cell_a", ])), c(22, 44)))
staffli_obs <- extract_staffli_obs_df(fake_staffli)
stopifnot(identical(staffli_obs$sampleID, c("1", "1")))
stopifnot(identical(staffli_obs$sample_name, c("demo_section_a", "demo_section_a")))

if (requireNamespace("Matrix", quietly = TRUE)) {
  graph <- Matrix::Matrix(
    c(
      0, 1, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 1, 0
    ),
    nrow = 4,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(graph) <- colnames(graph) <- rownames(obs)[1:4]
  graph_edges <- extract_sparse_graph_edges_by_section(
    graph = graph,
    cell_names = rownames(obs)[1:4],
    group_values = c("section_a", "section_a", "section_b", "section_b")
  )
  stopifnot(identical(graph_edges$section_a, c(0L, 1L)))
  stopifnot(identical(graph_edges$section_b, c(0L, 1L)))
}

wx_payload <- build_viewer_payload(
  input = merged_toy,
  groupby = "sample_id",
  initial_color = "cell_type",
  genes = c("Gene01", "Gene02"),
  neighbor_mode = "existing",
  marker_genes_groupby = "auto",
  marker_test = "wilcoxon"
)
stopifnot("Gene01" %in% unlist(wx_payload$marker_genes$cell_type$A, use.names = FALSE))

zs_payload <- build_viewer_payload(
  input = merged_toy,
  groupby = "sample_id",
  initial_color = "cell_type",
  neighbor_mode = "existing",
  neighbor_stats_permutations = 50L,
  neighbor_stats_seed = 1L
)
stopifnot(!is.null(zs_payload$neighbor_stats$cell_type$zscore))
stopifnot(zs_payload$neighbor_stats$cell_type$perm_n == 50L)

output_path <- tempfile(fileext = ".html")
render_viewer_html(
  payload = payload,
  output_path = output_path,
  title = "Smoke Viewer",
  viewer_shell_path = "inst/viewer/karospace_viewer_shell.html"
)

html <- paste(readLines(output_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
stopifnot(grepl("Smoke Viewer", html, fixed = TRUE))
stopifnot(grepl("section_a", html, fixed = TRUE))
stopifnot(grepl("Gene01", html, fixed = TRUE))
stopifnot(grepl("marker_genes", html, fixed = TRUE))
stopifnot(grepl("interaction_markers", html, fixed = TRUE))

config_dir <- tempfile("karospacer-config-")
dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
config_input_path <- file.path(config_dir, "toy_input.rds")
config_output_path <- file.path(config_dir, "toy_from_config.html")
config_path <- file.path(config_dir, "build.json")
saveRDS(merged_toy, config_input_path)

write_karospace_build_config(
  config = list(
    version = 1L,
    input = config_input_path,
    output = config_output_path,
    groupby = "sample_id",
    initial_color = "cell_type",
    additional_colors = c("course", "seurat_clusters"),
    genes = c("Gene01", "Gene02"),
    neighbor_mode = "existing",
    marker_genes_groupby = "auto",
    marker_test = "wilcoxon",
    neighbor_stats_permutations = 5L,
    neighbor_stats_seed = 1L,
    title = "Config Viewer"
  ),
  path = config_path
)
read_config <- read_karospace_build_config(config_path)
stopifnot(identical(read_config$groupby, "sample_id"))
stopifnot(identical(read_config$initial_color, "cell_type"))
stopifnot(identical(read_config$output, normalizePath(config_output_path, mustWork = FALSE)))

export_karospace_viewer_from_config(config_path)
config_html <- paste(readLines(config_output_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
stopifnot(grepl("Config Viewer", config_html, fixed = TRUE))
stopifnot(grepl("Gene01", config_html, fixed = TRUE))
stopifnot(grepl("neighbor_stats", config_html, fixed = TRUE))

inspect_summary <- inspect_karospace_input(
  input = config_input_path,
  gene_query = "Gene01",
  gene_limit = 2L
)
stopifnot(identical(inspect_summary$defaults$groupby, "sample_id"))
stopifnot(identical(inspect_summary$defaults$initial_color, "cell_type"))
stopifnot(identical(inspect_summary$defaults$neighbor_mode, "spatial"))
stopifnot(isTRUE(inspect_summary$genes$matched_count == 1L))
stopifnot(identical(inspect_summary$default_config$groupby, "sample_id"))
stopifnot(is.list(inspect_summary$default_config$marker_genes_groupby))
stopifnot(identical(unlist(inspect_summary$default_config$marker_genes_groupby, use.names = FALSE), "auto"))

cat("smoke_export: ok\n")
