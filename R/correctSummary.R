#' Apply PALM compositional correction to Step2 result files
#'
#' @description
#' Reads per-feature Step2 result files produced by [getSummary()] or
#' `extdata/step2_1_summary.R`, then applies the compositional correction logic
#' from `thirdParty/palm.get.summary.R` across features for each SNP. This is
#' intended for a post-processing "Step2.3" stage after the Step2.2
#' chromosome-merge step and after all feature-level Step2 jobs have completed.
#'
#' @param inputPrefix Step2 file prefix. Files are expected at
#'   `<inputPrefix>_<feature>.txt`.
#' @param outputPrefix Output prefix for corrected files. Defaults to
#'   `inputPrefix`, which overwrites the original Step2 files.
#' @param NULLmodelFile Optional Step1 NULL-model `.rda` file containing
#'   `modglmm`. When provided, the function mirrors the PALM behavior of forcing
#'   correction to be skipped if any study has `abd = TRUE`.
#'
#' @return Invisibly returns the corrected file paths.
#' @export
correctSummary <- function(inputPrefix,
                           outputPrefix = inputPrefix,
                           NULLmodelFile = NULL) {
  if (!requireNamespace("PALM", quietly = TRUE)) {
    stop("Package 'PALM' is required but not installed.")
  }

  if (missing(inputPrefix) || !nzchar(inputPrefix)) {
    stop("'inputPrefix' must be provided.")
  }
  if (missing(outputPrefix) || !nzchar(outputPrefix)) {
    stop("'outputPrefix' must be provided.")
  }

  escape_regex <- function(x) {
    gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
  }

  step2_dir <- dirname(inputPrefix)
  step2_base <- basename(inputPrefix)
  pattern <- paste0("^", escape_regex(step2_base), "_.*[.]txt$")
  files <- list.files(step2_dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) {
    stop("No Step2 files found for inputPrefix: ", inputPrefix)
  }

  extract_feature <- function(path) {
    base <- basename(path)
    sub(paste0("^", escape_regex(step2_base), "_(.*)[.]txt$"), "\\1", base)
  }
  features <- vapply(files, extract_feature, character(1))
  if (anyDuplicated(features)) {
    dup <- unique(features[duplicated(features)])
    stop("Duplicated feature names detected for inputPrefix: ", paste(dup, collapse = ", "))
  }
  names(files) <- features

  copy_step2_outputs <- function() {
    out_dir <- dirname(outputPrefix)
    if (!out_dir %in% c("", ".")) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    out_paths <- character(length(files))
    names(out_paths) <- names(files)
    for (feat in names(files)) {
      src <- files[[feat]]
      dst <- paste0(outputPrefix, "_", feat, ".txt")
      if (!identical(normalizePath(src, winslash = "/", mustWork = FALSE),
                     normalizePath(dst, winslash = "/", mustWork = FALSE))) {
        ok <- file.copy(src, dst, overwrite = TRUE)
        if (!isTRUE(ok)) {
          stop("Failed to copy uncorrected Step2 file to ", dst)
        }
      }
      out_paths[[feat]] <- dst
    }
    invisible(unname(out_paths))
  }

  if (!is.null(NULLmodelFile)) {
    if (!file.exists(NULLmodelFile)) {
      stop("NULL model file not found: ", NULLmodelFile)
    }
    env <- new.env()
    load(NULLmodelFile, envir = env)
    if (!exists("modglmm", envir = env)) {
      stop("Object 'modglmm' not found in ", NULLmodelFile)
    }
    modglmm <- env$modglmm

    if (any(vapply(modglmm, function(d) isTRUE(d$abd), logical(1)))) {
      warning(
        "At least one study includes features with an average proportion larger than 90%. ",
        "Skipping correction to match PALM behavior."
      )
      return(copy_step2_outputs())
    }

  }

  read_step2_one <- function(path) {
    dat <- utils::read.table(
      path,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    need <- c("SNP", "CHR", "POS", "est", "stderr", "pval")
    miss <- setdiff(need, colnames(dat))
    if (length(miss) > 0L) {
      stop("Missing required columns in ", path, ": ", paste(miss, collapse = ", "))
    }
    dat <- dat[, need, drop = FALSE]
    dat$SNP <- as.character(dat$SNP)
    dat$CHR <- suppressWarnings(as.integer(dat$CHR))
    dat$POS <- suppressWarnings(as.integer(dat$POS))
    dat$est <- suppressWarnings(as.numeric(dat$est))
    dat$stderr <- suppressWarnings(as.numeric(dat$stderr))
    dat$pval <- suppressWarnings(as.numeric(dat$pval))
    dat
  }

  feature_data <- lapply(files, read_step2_one)
  feature_ids <- names(feature_data)
  snp_ids <- unique(unlist(lapply(feature_data, `[[`, "SNP"), use.names = FALSE))

  est_mat <- matrix(
    NA_real_,
    nrow = length(feature_ids),
    ncol = length(snp_ids),
    dimnames = list(feature_ids, snp_ids)
  )
  stderr_mat <- est_mat

  for (feat in feature_ids) {
    dat <- feature_data[[feat]]
    idx <- match(dat$SNP, snp_ids)
    est_mat[feat, idx] <- dat$est
    stderr_mat[feat, idx] <- dat$stderr
  }

  delta <- apply(est_mat, 2, function(x) stats::median(-x, na.rm = TRUE))

  adjust_part <- vapply(snp_ids, function(snp) {
    se <- stderr_mat[, snp]
    non_na <- !is.na(se)
    if (!any(non_na)) {
      return(NA_real_)
    }
    sum(non_na) / (2 * sum(1 / sqrt(2 * pi) / se[non_na]))^2
  }, numeric(1))
  names(adjust_part) <- snp_ids

  out_dir <- dirname(outputPrefix)
  if (!out_dir %in% c("", ".")) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  out_paths <- character(length(feature_ids))
  names(out_paths) <- feature_ids
  for (feat in feature_ids) {
    dat <- feature_data[[feat]]
    idx <- match(dat$SNP, snp_ids)
    ok <- !is.na(dat$est) & !is.na(dat$stderr) & !is.na(delta[idx]) & !is.na(adjust_part[idx])
    dat$est[ok] <- dat$est[ok] + delta[idx][ok]
    dat$stderr[ok] <- sqrt(dat$stderr[ok]^2 + adjust_part[idx][ok])
    dat$pval <- 1 - stats::pchisq((dat$est / dat$stderr)^2, df = 1)

    out_file <- paste0(outputPrefix, "_", feat, ".txt")
    utils::write.table(
      dat,
      file = out_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
    out_paths[[feat]] <- out_file
  }

  message(
    "Step2 correction finished: wrote ", length(out_paths),
    " corrected file(s) with prefix ", outputPrefix
  )
  invisible(unname(out_paths))
}
