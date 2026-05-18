#' Fit a PALM null model
#'
#' Reads abundance and optional covariate tables, fits `PALM::palm.null.model`,
#' saves the fitted object, and invisibly returns it.
#'
#' @param abdFile Path to abundance table with subject IDs in the first column.
#' @param phenoColList Optional phenotype column names to keep from `abdFile`.
#'   Accepts either a character vector or a single comma-separated string.
#'   By default, all non-ID columns in `abdFile` are used.
#' @param covFile Optional path to covariate table with matching subject IDs.
#'   Use `NULL` (default) to fit without covariates.
#' @param covarColList Optional covariate column names to keep from `covFile`.
#'   Accepts either a character vector or a single comma-separated string.
#'   By default, all non-ID columns in `covFile` are used. If `depthCol` is
#'   provided, it is excluded from the covariate-adjustment matrix.
#' @param depthCol Optional column name in `covFile` to use as sequencing depth.
#'   If not provided, `depth = NULL` is passed to PALM so sequencing depth is
#'   computed from row sums of `abdFile`.
#' @param clusterCol Optional column name in `covFile` to use as the
#'   family/pedigree cluster ID for Step2.1. This column is stored in the saved
#'   null model and excluded from `covariate.adjust`.
#' @param prev.filter Passed to `PALM::palm.null.model()`; features with
#'   prevalence less than or equal to this threshold are removed. Defaults to `0.1`.
#' @param featureInfoFile Optional output path for feature prevalence and
#'   average proportion computed from the final Step1 modeled feature set after
#'   prevalence filtering. Use `NULL` (default) to skip writing this file.
#' @param nullModelPrefix Output prefix for the saved null model object. The
#'   function writes `<nullModelPrefix>.rda`.
#'
#' @return Invisibly returns the fitted PALM null model object.
#' @export
fitNULL <- function(abdFile,
                    phenoColList = NULL,
                    covFile = NULL,
                    covarColList = NULL,
                    depthCol = NULL,
                    clusterCol = NULL,
                    prev.filter = 0.1,
                    featureInfoFile = NULL,
                    nullModelPrefix) {
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
  if (is.null(covFile) && !is.null(clusterCol)) {
    stop("'clusterCol' requires a non-NULL 'covFile'.")
  }
  if (missing(nullModelPrefix) || !nzchar(nullModelPrefix)) {
    stop("'nullModelPrefix' must be provided.")
  }
  null_model_file <- nullModelPrefix
  if (grepl("\\.rda$", null_model_file, ignore.case = TRUE)) {
    null_model_file <- sub("\\.rda$", "", null_model_file, ignore.case = TRUE)
    message("nullModelPrefix should not include .rda; normalizing to prefix: ", null_model_file)
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
  format_name_list <- function(x, max_show = 20L) {
    if (!length(x)) {
      return("<none>")
    }
    shown <- utils::head(x, max_show)
    suffix <- if (length(x) > max_show) paste0(", ... (+", length(x) - max_show, " more)") else ""
    paste0(paste(shown, collapse = ", "), suffix)
  }
  read_firstcol_as_data <- function(file) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package 'data.table' is required but not installed.")
    }
    df <- data.table::fread(file = file, data.table = FALSE, check.names = FALSE)
    if (ncol(df) < 1L) {
      stop("Input file has no columns: ", file)
    }
    ids <- as.character(df[[1]])
    if (anyNA(ids) || any(!nzchar(ids))) {
      stop("Missing/empty subject ID detected in first column of ", file)
    }
    df[[1]] <- NULL
    list(data = df, ids = ids)
  }

  abd_input <- read_firstcol_as_data(abdFile)
  abd <- abd_input$data
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
  subject_ids <- abd_input$ids
  model_row_ids <- make.unique(subject_ids, sep = "__rep")
  abd <- as.matrix(abd)
  storage.mode(abd) <- "numeric"
  rownames(abd) <- model_row_ids
  message(
    "Input abundance matrix: ", nrow(abd), " rows x ", ncol(abd), " features."
  )
  cov <- NULL
  if (!is.null(covFile)) {
    if (!file.exists(covFile)) {
      stop("'covFile' does not exist: ", covFile)
    }
    cov_input <- read_firstcol_as_data(covFile)
    cov <- cov_input$data
    rownames(cov) <- make.unique(cov_input$ids)
    message(
      "Input covariate table: ", nrow(cov), " rows x ", ncol(cov),
      " covariate column(s)."
    )
    message("Covariate columns in covFile: ", format_name_list(colnames(cov)))
    missing_cov_ids <- setdiff(unique(subject_ids), unique(cov_input$ids))
    extra_cov_ids <- setdiff(unique(cov_input$ids), unique(subject_ids))
    same_sample_order <- identical(subject_ids, cov_input$ids)
    message(
      "Covariate subject match: matched rows=", sum(subject_ids == cov_input$ids),
      "/", nrow(abd),
      ", missing=", length(missing_cov_ids),
      ", extra=", length(extra_cov_ids),
      ", same_order=", same_sample_order
    )
    if (length(missing_cov_ids) > 0L) {
      message("Covariate missing subject IDs: ", format_name_list(missing_cov_ids, max_show = 5L))
    }
    if (length(extra_cov_ids) > 0L) {
      message("Covariate extra subject IDs: ", format_name_list(extra_cov_ids, max_show = 5L))
    }
    if (!same_sample_order) {
      stop("Step1 requires covFile rows to be aligned to abdFile rows by subject ID. Run Step0 first.")
    }
  }

  depthCol <- normalize_col_list(depthCol, "depthCol")
  if (!is.null(depthCol) && length(depthCol) != 1L) {
    stop("'depthCol' must specify exactly one column name.")
  }

  depth <- NULL
  cluster_ids <- NULL
  if (!is.null(depthCol)) {
    if (!(depthCol %in% colnames(cov))) {
      stop("Requested depth column not found in 'covFile': ", depthCol)
    }
    depth <- suppressWarnings(as.numeric(cov[[depthCol]]))
    if (any(is.na(depth) & !is.na(cov[[depthCol]]))) {
      stop("Requested depth column in 'covFile' cannot be safely converted to numeric: ", depthCol)
    }
    if (length(depth) != nrow(abd)) {
      stop("Depth column length does not match abundance row count. Run Step0 first.")
    }
    names(depth) <- subject_ids
    message("depthCol provided: using '", depthCol, "' from ", covFile, " as sequencing depth.")
    message("Depth summary min/median/max=", min(depth), "/", stats::median(depth), "/", max(depth))
  } else {
    message("No depthCol provided: PALM will use row sums of abundance as depth.")
  }
  clusterCol <- normalize_col_list(clusterCol, "clusterCol")
  if (!is.null(clusterCol) && length(clusterCol) != 1L) {
    stop("'clusterCol' must specify exactly one column name.")
  }
  if (!is.null(clusterCol)) {
    if (!(clusterCol %in% colnames(cov))) {
      stop("Requested cluster column not found in 'covFile': ", clusterCol)
    }
    cluster_ids <- as.character(cov[[clusterCol]])
    if (anyNA(cluster_ids) || any(!nzchar(trimws(cluster_ids)))) {
      stop("Requested cluster column in 'covFile' contains missing or empty values: ", clusterCol)
    }
    names(cluster_ids) <- model_row_ids
    message("clusterCol provided: using '", clusterCol, "' from ", covFile, " as Step2.1 cluster ID.")
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
      message("Requested covariate columns from covarColList: ", format_name_list(covarColList))
      missing_cols <- setdiff(covarColList, colnames(cov))
      if (length(missing_cols) > 0) {
        stop(
          "Requested covariate column(s) not found in 'covFile': ",
          paste(missing_cols, collapse = ", ")
        )
      }
      cov <- cov[, covarColList, drop = FALSE]
      message("Covariate columns after covarColList selection: ", format_name_list(colnames(cov)))
    } else {
      message("covarColList is NULL: starting from all covariate columns in covFile.")
    }
    if (!is.null(depthCol) && depthCol %in% colnames(cov)) {
      cov <- cov[, setdiff(colnames(cov), depthCol), drop = FALSE]
      message("Excluding depthCol '", depthCol, "' from covariate.adjust.")
    } else if (!is.null(depthCol)) {
      message("depthCol '", depthCol, "' is not in selected covariate columns; no covariate exclusion needed.")
    } else {
      message("depthCol is NULL: no covariate column is excluded as sequencing depth.")
    }
    if (!is.null(clusterCol) && clusterCol %in% colnames(cov)) {
      cov <- cov[, setdiff(colnames(cov), clusterCol), drop = FALSE]
      message("Excluding clusterCol '", clusterCol, "' from covariate.adjust.")
    }
    message("Final covariate.adjust columns: ", format_name_list(colnames(cov)))
    if (ncol(cov) == 0L) {
      message("covFile provided, but no covariate columns remain after excluding depthCol.")
      modglmm <- PALM::palm.null.model(
        rel.abd = abd,
        depth = depth,
        prev.filter = prev.filter
      )
    } else {
      rownames(cov) <- model_row_ids
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
  attr(modglmm, "subject_ids") <- subject_ids
  attr(modglmm, "model_row_ids") <- model_row_ids
  if (!is.null(cluster_ids)) {
    attr(modglmm, "cluster_ids") <- cluster_ids
    attr(modglmm, "cluster_col") <- clusterCol
  }
  save(modglmm, file = null_model_file)
  message("Done. PALM null model saved to ", null_model_file)

  feature_ids <- unique(unlist(lapply(modglmm, function(x) colnames(x$Y_I)), use.names = FALSE))
  if (!length(feature_ids)) {
    stop("No modeled features found in fitted null model.")
  }

  if (!is.null(featureInfoFile) && nzchar(featureInfoFile)) {
    missing_modeled_features <- setdiff(feature_ids, colnames(abd))
    if (length(missing_modeled_features) > 0L) {
      stop(
        "Modeled feature(s) not found in Step1 abundance input: ",
        paste(utils::head(missing_modeled_features, 5), collapse = ", ")
      )
    }
    message(
      "Generating feature info from the final Step1 modeled feature set after prev.filter. ",
      "Output path: ", featureInfoFile
    )
    feature_stats <- feature_info_from_matrix(abd[, feature_ids, drop = FALSE], feature_ids = feature_ids)
    dir.create(dirname(featureInfoFile), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(
      feature_stats,
      file = featureInfoFile,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )
    message("Feature info finished: ", nrow(feature_stats), " modeled feature(s) written to ", featureInfoFile)
  } else {
    message("Feature info skipped: featureInfoFile is NULL.")
  }

  invisible(modglmm)
}
