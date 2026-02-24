#' Fit a PALM null model
#'
#' Reads abundance and optional covariate tables, fits `PALM::palm.null.model`,
#' saves the fitted object, and invisibly returns it.
#'
#' @param abdFile Path to abundance table with sample IDs in the first column.
#' @param covFile Optional path to covariate table with matching sample IDs.
#'   Use `NULL` (default) to fit without covariates.
#' @param outputPrefix File prefix (without extension) for the saved `.rda` file.
#'
#' @return Invisibly returns the fitted PALM null model object.
#' @export
fitNULL <- function(abdFile,
                    covFile = NULL,
                    outputPrefix) {
  if (!requireNamespace("PALM", quietly = TRUE)) {
    stop("Package 'PALM' is required but not installed.")
  }

  if (missing(abdFile) || !nzchar(abdFile)) {
    stop("'abdFile' must be provided.")
  }
  if (!file.exists(abdFile)) {
    stop("'abdFile' does not exist: ", abdFile)
  }
  if (missing(outputPrefix) || !nzchar(outputPrefix)) {
    stop("'outputPrefix' must be provided.")
  }

  abd <- read_firstcol_as_rownames(abdFile)
  abd <- as.matrix(abd)

  if (is.null(covFile)) {
    message("Fitting PALM null model without covariates.")
    modglmm <- PALM::palm.null.model(
      rel.abd = abd,
      prev.filter = 0
    )
  } else {
    if (!file.exists(covFile)) {
      stop("'covFile' does not exist: ", covFile)
    }
    cov <- read_firstcol_as_rownames(covFile)
    modglmm <- PALM::palm.null.model(
      rel.abd = abd,
      covariate.adjust = cov,
      prev.filter = 0
    )
  }

  dir.create(dirname(outputPrefix), recursive = TRUE, showWarnings = FALSE)
  save(modglmm, file = paste0(outputPrefix, ".rda"))
  message("Done. PALM null model saved to ", paste0(outputPrefix, ".rda"))
  invisible(modglmm)
}
