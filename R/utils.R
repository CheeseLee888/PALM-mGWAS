#' Read a tabular file and set the first column as row names
#'
#' Convenience helper used throughout the package for reading abundance and
#' covariate tables that store sample IDs in the first column.
#'
#' @param file Path to a delimited text file readable by `data.table::fread()`.
#'
#' @return A `data.frame` with row names taken from the first column and that
#'   column removed.
#' @export
read_firstcol_as_rownames <- function(file) {
  if (missing(file) || !nzchar(file)) {
    stop("'file' must be provided.")
  }
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  df <- data.table::fread(
    file = file,
    data.table = FALSE,
    check.names = FALSE
  )

  if (ncol(df) < 1) {
    stop("Input file has no columns: ", file)
  }

  id <- as.character(df[[1]])
  rownames(df) <- id
  df[[1]] <- NULL

  df
}

# Read a named vector from a tabular file whose first column stores IDs.
# If `column` is NULL, use the first non-ID column.
read_named_vector <- function(file, column = NULL, numeric = FALSE) {
  if (missing(file) || !nzchar(file)) {
    stop("'file' must be provided.")
  }
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  df <- data.table::fread(
    file = file,
    data.table = FALSE,
    check.names = FALSE
  )

  if (ncol(df) < 2) {
    stop("Expected at least two columns in: ", file)
  }

  ids <- as.character(df[[1]])
  if (anyDuplicated(ids)) {
    stop("Duplicate IDs detected in first column of ", file)
  }

  if (is.null(column) || !nzchar(as.character(column))) {
    value_idx <- 2L
  } else if (is.numeric(column)) {
    value_idx <- as.integer(column)
  } else if (grepl("^[0-9]+$", as.character(column))) {
    value_idx <- as.integer(column)
  } else if (column %in% colnames(df)) {
    value_idx <- match(column, colnames(df))
  } else {
    stop("Column '", column, "' not found in ", file)
  }

  if (value_idx < 2L || value_idx > ncol(df)) {
    stop("Selected column must refer to a non-ID column in ", file)
  }

  values <- df[[value_idx]]
  names(values) <- ids

  if (numeric) {
    values_num <- suppressWarnings(as.numeric(values))
    if (any(is.na(values_num) & !is.na(values))) {
      stop("Selected column in ", file, " cannot be safely converted to numeric.")
    }
    values <- values_num
    names(values) <- ids
  }

  values
}
