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
    data.table   = FALSE,
    check.names  = FALSE
  )

  if (ncol(df) < 1) {
    stop("Input file has no columns: ", file)
  }

  id <- as.character(df[[1]])
  rownames(df) <- id
  df[[1]] <- NULL

  df
}
