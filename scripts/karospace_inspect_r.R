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

    next_is_value <- i < length(args) && !startsWith(args[[i + 1L]], "--")
    if (next_is_value) {
      out[[substring(key, 3L)]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      out[[substring(key, 3L)]] <- TRUE
      i <- i + 1L
    }
  }
  out
}

split_csv <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return(NULL)
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

options <- parse_args(args)

if (isTRUE(options$help) || is.null(options$input)) {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/karospace_inspect_r.R --input path/to/object.rds",
      "[--metadata-input metadata.rds] [--metadata-input-columns col1,col2] [--metadata-prefix ext_]",
      "[--assay SCT] [--gene-query COL] [--gene-limit 50] [--output inspect.json]",
      sep = "\n"
    ),
    "\n"
  )
  quit(save = "no", status = if (isTRUE(options$help)) 0L else 1L)
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(repo_root, "R", "helpers.R"))
source(file.path(repo_root, "R", "source.R"))
source(file.path(repo_root, "R", "payload.R"))
source(file.path(repo_root, "R", "export.R"))
source(file.path(repo_root, "R", "inspect.R"))

result <- inspect_karospace_input(
  input = normalizePath(options$input, mustWork = TRUE),
  metadata_input = if (!is.null(options[["metadata-input"]])) normalizePath(options[["metadata-input"]], mustWork = TRUE) else NULL,
  metadata_input_columns = split_csv(options[["metadata-input-columns"]]),
  metadata_prefix = options[["metadata-prefix"]],
  assay = options[["assay"]],
  gene_query = options[["gene-query"]],
  gene_limit = suppressWarnings(as.integer(options[["gene-limit"]] %||% 50L))
)

json <- jsonlite::toJSON(
  result,
  auto_unbox = TRUE,
  null = "null",
  na = "null",
  digits = NA,
  pretty = TRUE
)

if (!is.null(options$output)) {
  output_path <- resolve_output_path(options$output)
  writeLines(json, con = output_path, useBytes = TRUE)
  cat("Wrote inspect JSON to ", output_path, "\n", sep = "")
} else {
  cat(json, "\n", sep = "")
}
