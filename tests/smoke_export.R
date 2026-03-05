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

payload <- build_viewer_payload(
  input = toy,
  groupby = "sample_id",
  initial_color = "cell_type",
  additional_colors = c("course"),
  genes = c("Gene01", "Gene02"),
  metadata_columns = c("course")
)

stopifnot(payload$n_sections == 3L)
stopifnot(payload$total_cells == n_cells)
stopifnot(identical(unlist(payload$available_colors, use.names = FALSE), c("cell_type", "course")))
stopifnot(identical(unlist(payload$available_genes, use.names = FALSE), c("Gene01", "Gene02")))
stopifnot(isTRUE(payload$has_umap))

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
