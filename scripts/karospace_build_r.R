#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }
    out[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

options <- parse_args(args)

required <- c("input", "output", "groupby", "initial-color")
missing <- required[!required %in% names(options)]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo_root, "R", "helpers.R"))
source(file.path(repo_root, "R", "source.R"))
source(file.path(repo_root, "R", "payload.R"))
source(file.path(repo_root, "R", "export.R"))

split_csv <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return(NULL)
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

export_karospace_viewer(
  input = options[["input"]],
  output_path = options[["output"]],
  groupby = options[["groupby"]],
  initial_color = options[["initial-color"]],
  additional_colors = split_csv(options[["additional-colors"]]),
  genes = split_csv(options[["genes"]]),
  metadata_columns = split_csv(options[["metadata-columns"]]),
  outline_by = options[["outline-by"]],
  title = options[["title"]] %||% "KaroSpace",
  theme = options[["theme"]] %||% "light",
  min_panel_size = as.numeric(options[["min-panel-size"]] %||% 150),
  spot_size = as.numeric(options[["spot-size"]] %||% 2)
)

cat("Wrote viewer to ", normalizePath(options[["output"]], mustWork = FALSE), "\n", sep = "")
