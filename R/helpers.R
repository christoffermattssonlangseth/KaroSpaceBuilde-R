`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

empty_named_list <- function() {
  stats::setNames(vector("list", 0), character(0))
}

compact_json <- function(x) {
  jsonlite::toJSON(
    x,
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    digits = NA,
    POSIXt = "ISO8601",
    force = TRUE
  )
}

as_json_array <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(list())
  }
  unname(as.list(x))
}

escape_html <- function(x) {
  if (is.null(x)) {
    return("")
  }

  value <- enc2utf8(as.character(x))
  value <- gsub("&", "&amp;", value, fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  value <- gsub(">", "&gt;", value, fixed = TRUE)
  value <- gsub("\"", "&quot;", value, fixed = TRUE)
  value
}

find_package_root <- function() {
  installed <- system.file(package = "KaroSpaceBuildR")
  if (nzchar(installed)) {
    return(installed)
  }

  cwd_candidate <- normalizePath(getwd(), mustWork = FALSE)
  if (file.exists(file.path(cwd_candidate, "DESCRIPTION"))) {
    return(cwd_candidate)
  }

  stop(
    "Could not resolve the package root. Install KaroSpaceBuildR or run from the repo root."
  )
}

resolve_asset_path <- function(rel_path, explicit_path = NULL) {
  if (!is.null(explicit_path)) {
    return(normalizePath(explicit_path, mustWork = TRUE))
  }

  installed <- system.file(rel_path, package = "KaroSpaceBuildR")
  if (nzchar(installed)) {
    return(installed)
  }

  repo_path <- file.path(find_package_root(), rel_path)
  normalizePath(repo_path, mustWork = TRUE)
}

resolve_output_path <- function(path) {
  outdir <- dirname(path)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, mustWork = FALSE)
}

unique_non_missing <- function(x) {
  values <- as.character(x)
  values <- values[!is.na(values) & nzchar(values)]
  unique(values)
}

first_non_missing <- function(x) {
  values <- unique_non_missing(x)
  if (length(values) == 0) {
    NULL
  } else {
    values[[1]]
  }
}

as_plain_matrix <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.matrix(x)) {
    return(x)
  }
  if (is.data.frame(x)) {
    return(as.matrix(x))
  }
  if (methods::is(x, "Matrix")) {
    return(x)
  }
  stop("Unsupported matrix-like object: ", paste(class(x), collapse = ", "))
}

coalesce_character <- function(x, y = character()) {
  out <- x %||% y
  as.character(out)
}
