source("R/helpers.R")
source("R/source.R")
source("R/payload.R")
source("R/export.R")

set.seed(7)
n_cells <- 60
n_genes <- 12

obs <- data.frame(
  sample_id = rep(c("section_a", "section_b", "section_c"), each = 20),
  cell_type = sample(c("A", "B", "C"), size = n_cells, replace = TRUE),
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

toy <- list(
  obs = obs,
  coordinates = coords,
  umap = umap,
  expression = expr
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
  metadata_columns = c("course")
)

stopifnot(payload$n_sections == 3L)
stopifnot(payload$total_cells == n_cells)
stopifnot(identical(unlist(payload$available_colors, use.names = FALSE), c("cell_type", "course", "seurat_clusters")))
stopifnot(identical(unlist(payload$available_genes, use.names = FALSE), c("Gene01", "Gene02")))
stopifnot(isTRUE(payload$has_umap))
stopifnot(identical(names(payload$genes_meta), c("Gene01", "Gene02")))
stopifnot(length(payload$sections[[1]]$genes$Gene01) == 20L)
stopifnot(identical(payload$gene_encodings$Gene01, "dense"))
stopifnot("(missing)" %in% unlist(payload$colors_meta$seurat_clusters$categories, use.names = FALSE))

gene_info <- resolve_input_gene_names(toy)
stopifnot(identical(unlist(gene_info$gene_names[1:3], use.names = FALSE), c("Gene01", "Gene02", "Gene03")))

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

cat("smoke_export: ok\n")
