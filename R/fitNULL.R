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
#'   By default, all non-ID columns in `covFile` are used. If `depthCol` is
#'   provided, it is excluded from the covariate-adjustment matrix.
#' @param depthCol Optional column name in `covFile` to use as sequencing depth.
#'   If not provided, `depth = NULL` is passed to PALM so sequencing depth is
#'   computed from row sums of `abdFile`.
#' @param prev.filter Passed to `PALM::palm.null.model()`; features with
#'   prevalence less than or equal to this threshold are removed. Defaults to `0.1`.
#' @param FeatureInfoFile Optional output path for feature prevalence and
#'   average proportion computed from the Step1 abundance input before
#'   prevalence filtering. Use `NULL` (default) to skip writing this file.
#' @param FeatureNameListFile Optional output path for the modeled feature IDs
#'   retained in the fitted Step1 null model after prevalence filtering. Use
#'   `NULL` (default) to skip writing this file.
#' @param NULLObjPrefix Output prefix for the saved null model object. The
#'   function writes `<NULLObjPrefix>.rda`.
#'
#' @return Invisibly returns the fitted PALM null model object.
#' @export
fitNULL <- function(abdFile,
                    phenoColList = NULL,
                    covFile = NULL,
                    covarColList = NULL,
                    depthCol = NULL,
                    prev.filter = 0.1,
                    FeatureInfoFile = NULL,
                    FeatureNameListFile = NULL,
                    NULLObjPrefix) {
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
  if (is.null(covFile) && !is.null(depthCol)) {
    stop("'depthCol' requires a non-NULL 'covFile'.")
  }
  if (missing(NULLObjPrefix) || !nzchar(NULLObjPrefix)) {
    stop("'NULLObjPrefix' must be provided.")
  }
  null_model_file <- NULLObjPrefix
  if (grepl("\\.rda$", null_model_file, ignore.case = TRUE)) {
    null_model_file <- sub("\\.rda$", "", null_model_file, ignore.case = TRUE)
    message("NULLObjPrefix should not include .rda; normalizing to prefix: ", null_model_file)
  }
  null_model_file <- paste0(null_model_file, ".rda")

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
  storage.mode(abd) <- "numeric"
  message(
    "Input abundance matrix: ", nrow(abd), " samples x ", ncol(abd), " features."
  )
  if (!is.null(FeatureInfoFile) && nzchar(FeatureInfoFile)) {
    message(
      "Generating FeatureInfo from the Step1 input abundance matrix before prev.filter. ",
      "Output path: ", FeatureInfoFile
    )
    feature_stats <- feature_info_from_matrix(abd, feature_ids = colnames(abd))
    dir.create(dirname(FeatureInfoFile), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(
      feature_stats,
      file = FeatureInfoFile,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )
    message("FeatureInfo finished: ", nrow(feature_stats), " feature(s) written to ", FeatureInfoFile)
  } else {
    message("FeatureInfo skipped: FeatureInfoFile is NULL.")
  }

  cov <- NULL
  if (!is.null(covFile)) {
    if (!file.exists(covFile)) {
      stop("'covFile' does not exist: ", covFile)
    }
    cov <- read_firstcol_as_rownames(covFile)
  }

  depthCol <- normalize_col_list(depthCol, "depthCol")
  if (!is.null(depthCol) && length(depthCol) != 1L) {
    stop("'depthCol' must specify exactly one column name.")
  }

  depth <- NULL
  if (!is.null(depthCol)) {
    if (!(depthCol %in% colnames(cov))) {
      stop("Requested depth column not found in 'covFile': ", depthCol)
    }
    depth <- suppressWarnings(as.numeric(cov[[depthCol]]))
    if (any(is.na(depth) & !is.na(cov[[depthCol]]))) {
      stop("Requested depth column in 'covFile' cannot be safely converted to numeric: ", depthCol)
    }
    names(depth) <- rownames(cov)
    missing_ids <- setdiff(rownames(abd), names(depth))
    if (length(missing_ids) > 0) {
      stop("Missing depth values for sample IDs: ", paste(utils::head(missing_ids, 5), collapse = ", "))
    }
    depth <- depth[rownames(abd)]
    message("depthCol provided: using '", depthCol, "' from ", covFile, " as sequencing depth.")
    message("Depth summary min/median/max=", min(depth), "/", stats::median(depth), "/", max(depth))
  } else {
    message("No depthCol provided: PALM will use row sums of abundance as depth.")
  }
  message("Prevalence filter setting: prev.filter=", prev.filter)

  if (is.null(cov)) {
    message("Fitting PALM null model without covariates.")
    modglmm <- PALM::palm.null.model(
      rel.abd = abd,
      depth = depth,
      prev.filter = prev.filter
    )
  } else {
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
    if (!is.null(depthCol) && depthCol %in% colnames(cov)) {
      cov <- cov[, setdiff(colnames(cov), depthCol), drop = FALSE]
      message("Excluding depthCol '", depthCol, "' from covariate.adjust.")
    }
    if (ncol(cov) == 0L) {
      message("covFile provided, but no covariate columns remain after excluding depthCol.")
      modglmm <- PALM::palm.null.model(
        rel.abd = abd,
        depth = depth,
        prev.filter = prev.filter
      )
    } else {
      message(
        "covFile provided: fitting PALM null model with ", ncol(cov), " covariate column(s) from ", covFile
      )
      modglmm <- PALM::palm.null.model(
        rel.abd = abd,
        covariate.adjust = cov,
        depth = depth,
        prev.filter = prev.filter
      )
    }
  }

  dir.create(dirname(null_model_file), recursive = TRUE, showWarnings = FALSE)
  save(modglmm, file = null_model_file)
  message("Done. PALM null model saved to ", null_model_file)

  if (!is.null(FeatureNameListFile) && nzchar(FeatureNameListFile)) {
    feature_ids <- unique(unlist(lapply(modglmm, function(x) colnames(x$Y_I)), use.names = FALSE))
    if (!length(feature_ids)) {
      stop("No modeled features found in fitted null model; cannot write FeatureNameListFile.")
    }
    dir.create(dirname(FeatureNameListFile), recursive = TRUE, showWarnings = FALSE)
    writeLines(feature_ids, FeatureNameListFile, useBytes = TRUE)
    message(
      "FeatureList finished: ", length(feature_ids),
      " feature(s) written to ", FeatureNameListFile
    )
  } else {
    message("FeatureList skipped: FeatureNameListFile is NULL.")
  }
  invisible(modglmm)
}
