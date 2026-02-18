read_firstcol_as_rownames <- function(file) {
  df <- data.table::fread(
    file = file,
    data.table   = FALSE,
    check.names  = FALSE
  )
  id <- as.character(df[[1]])
  rownames(df) <- id
  df[[1]] <- NULL
  return(df)
}
