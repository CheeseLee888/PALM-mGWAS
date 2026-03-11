#' Fit a PALM null model
#'
#' Reads abundance and optional covariate tables, fits `PALM::palm.null.model`,
#' saves the fitted object, and invisibly returns it.
#'
#' @param abdFile Path to abundance table with sample IDs in the first column.
#' @param phenoColList Optional phenotype column names to keep from `abdFile`.
#'   Accepts either a character vector or a single comma-separated string.
#'   By default, all non-ID columns in `abdFile` are used.
#' @param covFile Optional path to covariate table with matching sample IDs.
#'   Use `NULL` (default) to fit without covariates.
#' @param covarColList Optional covariate column names to keep from `covFile`.
#'   Accepts either a character vector or a single comma-separated string.
#'   By default, all non-ID columns in `covFile` are used.
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
                    phenoColList = NULL,
                    covFile = NULL,
                    covarColList = NULL,
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
  if (is.null(covFile) && !is.null(covarColList)) {
    stop("'covarColList' requires a non-NULL 'covFile'.")
  }
  if (missing(NULLmodelFile) || !nzchar(NULLmodelFile)) {
    stop("'NULLmodelFile' must be provided.")
  }
  if (!grepl("\\.rda$", NULLmodelFile, ignore.case = TRUE)) {
    stop("'NULLmodelFile' must include the .rda suffix.")
  }

  normalize_col_list <- function(x, arg_name) {
    if (is.null(x)) {
      return(NULL)
    }
    if (length(x) == 1L) {
      x <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
    } else {
      x <- trimws(as.character(x))
    }
    x <- x[nzchar(x)]
    if (!length(x)) {
      stop("'", arg_name, "' did not contain any valid column names.")
    }
    x
  }

  abd <- read_firstcol_as_rownames(abdFile)
  phenoColList <- normalize_col_list(phenoColList, "phenoColList")
  if (!is.null(phenoColList)) {
    missing_cols <- setdiff(phenoColList, colnames(abd))
    if (length(missing_cols) > 0) {
      stop(
        "Requested phenotype column(s) not found in 'abdFile': ",
        paste(missing_cols, collapse = ", ")
      )
    }
    abd <- abd[, phenoColList, drop = FALSE]
  }
  abd <- as.matrix(abd)
  message(
    "Input abundance matrix: ", nrow(abd), " samples x ", ncol(abd), " features."
  )

  depth <- NULL
  if (!is.null(depthFile)) {
    message("depthFile provided: using the second column as sequencing depth from ", depthFile)
    depth <- read_named_vector(depthFile, numeric = TRUE)
    missing_ids <- setdiff(rownames(abd), names(depth))
    if (length(missing_ids) > 0) {
      stop("Missing depth values for sample IDs: ", paste(utils::head(missing_ids, 5), collapse = ", "))
    }
    depth <- depth[rownames(abd)]
    message(
      "Depth filtering enabled. depth.filter=", depth.filter,
      "; depth summary min/median/max=",
      min(depth), "/", stats::median(depth), "/", max(depth)
    )
  } else {
    message("No depthFile provided: PALM will use row sums of abundance as depth.")
  }
  message("Prevalence filter setting: prev.filter=", prev.filter)

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
    covarColList <- normalize_col_list(covarColList, "covarColList")
    if (!is.null(covarColList)) {
      missing_cols <- setdiff(covarColList, colnames(cov))
      if (length(missing_cols) > 0) {
        stop(
          "Requested covariate column(s) not found in 'covFile': ",
          paste(missing_cols, collapse = ", ")
        )
      }
      cov <- cov[, covarColList, drop = FALSE]
    }
    message(
      "covFile provided: fitting PALM null model with ", ncol(cov), " covariate column(s) from ", covFile
    )
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
