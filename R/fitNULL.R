#' Fit a PALM null model
#'
#' Reads abundance and optional covariate tables, fits `PALM::palm.null.model`,
#' saves the fitted object, and invisibly returns it.
#'
#' @param abdFile Path to abundance table with sample IDs in the first column.
#' @param covFile Optional path to covariate table with matching sample IDs.
#'   Use `NULL` (default) to fit without covariates.
#' @param depthFile Optional path to a sequencing-depth table with sample IDs in
#'   the first column and depth values in the second column.
#' @param depth.filter Passed to `PALM::palm.null.model()`; samples with depth
#'   less than or equal to this threshold are removed. Defaults to `0`.
#' @param prev.filter Passed to `PALM::palm.null.model()`; features with
#'   prevalence less than or equal to this threshold are removed. Defaults to `0.1`.
#' @param NULLmodelFile Full output path for the saved `.rda` file.
#'
#' @return Invisibly returns the fitted PALM null model object.
#' @export
fitNULL <- function(abdFile,
                    covFile = NULL,
                    depthFile = NULL,
                    depth.filter = 0,
                    prev.filter = 0.1,
                    NULLmodelFile) {
  if (!requireNamespace("PALM", quietly = TRUE)) {
    stop("Package 'PALM' is required but not installed.")
  }

  if (missing(abdFile) || !nzchar(abdFile)) {
    stop("'abdFile' must be provided.")
  }
  if (!file.exists(abdFile)) {
    stop("'abdFile' does not exist: ", abdFile)
  }
  if (missing(NULLmodelFile) || !nzchar(NULLmodelFile)) {
    stop("'NULLmodelFile' must be provided.")
  }
  if (!grepl("\\.rda$", NULLmodelFile, ignore.case = TRUE)) {
    stop("'NULLmodelFile' must include the .rda suffix.")
  }

  abd <- read_firstcol_as_rownames(abdFile)
  abd <- as.matrix(abd)

  depth <- NULL
  if (!is.null(depthFile)) {
    depth <- read_named_vector(depthFile, numeric = TRUE)
    missing_ids <- setdiff(rownames(abd), names(depth))
    if (length(missing_ids) > 0) {
      stop("Missing depth values for sample IDs: ", paste(utils::head(missing_ids, 5), collapse = ", "))
    }
    depth <- depth[rownames(abd)]
  }

  if (is.null(covFile)) {
    message("Fitting PALM null model without covariates.")
    modglmm <- PALM::palm.null.model(
      rel.abd = abd,
      depth = depth,
      depth.filter = depth.filter,
      prev.filter = prev.filter
    )
  } else {
    if (!file.exists(covFile)) {
      stop("'covFile' does not exist: ", covFile)
    }
    cov <- read_firstcol_as_rownames(covFile)
    modglmm <- PALM::palm.null.model(
      rel.abd = abd,
      covariate.adjust = cov,
      depth = depth,
      depth.filter = depth.filter,
      prev.filter = prev.filter
    )
  }

  dir.create(dirname(NULLmodelFile), recursive = TRUE, showWarnings = FALSE)
  save(modglmm, file = NULLmodelFile)
  message("Done. PALM null model saved to ", NULLmodelFile)
  invisible(modglmm)
}
