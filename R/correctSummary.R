#' Apply PALM compositional correction to Step2 result files
#'
#' @description
#' Reads per-feature Step2 result files produced by [getSummary()] or
#' `extdata/step2_1_summary.R`, then applies the compositional correction logic
#' from `thirdParty/palm.get.summary.R` across features for each SNP. This is
#' intended for a post-processing "Step2.2" stage after all feature-level Step2
#' jobs have completed. It can operate on either one chromosome suffix
#' (`_chr1`..`_chr22`) or the merged `_allchr` files.
#'
#' @param inputPrefix Shared Step2 base prefix. Files are expected at
#'   `<inputPrefix>_chr<chrom>_<feature>.txt` when `chrom` is set, or
#'   `<inputPrefix>_allchr_<feature>.txt` when `chrom` is `NULL`. An optional
#'   trailing underscore is ignored.
#' @param chrom Optional chromosome selector. Use `NULL` to correct the merged
#'   `_allchr` files. Use `1`..`22` or strings like `"chr1"` to correct one
#'   chromosome-specific shard across all features.
#' @param overwriteOutput Logical; if `TRUE` overwrite the original Step2 files.
#'   If `FALSE`, keep the originals and write new files as
#'   `<inputPrefix>_corrected_<allchr|chrN>_<feature>.txt`.
#' @param correct Correction method. Supported values are `"median"` and
#'   `"tune"`.
#' @param NULLmodelFile Optional Step1 NULL model `.rda` containing `modglmm`.
#'   Required when `correct = "tune"`; the sample size is inferred from
#'   `nrow(modglmm[[1]]$Y_I)`.
#'
#' @return Invisibly returns the corrected file paths.
#' @export
correctSummary <- function(inputPrefix,
                           chrom = NULL,
                           overwriteOutput = TRUE,
                           correct = c("median", "tune"),
                           NULLmodelFile = NULL) {
  if (!requireNamespace("PALM", quietly = TRUE)) {
    stop("Package 'PALM' is required but not installed.")
  }

  correct <- match.arg(correct)
  if (missing(inputPrefix) || !nzchar(inputPrefix)) {
    stop("'inputPrefix' must be provided.")
  }
  if (!is.null(chrom) && length(chrom) != 1L) {
    stop("'chrom' must be NULL or a single chromosome value.")
  }
  if (!is.logical(overwriteOutput) || length(overwriteOutput) != 1L || is.na(overwriteOutput)) {
    stop("'overwriteOutput' must be a single TRUE/FALSE value.")
  }
  tuneN <- NULL
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
    if (identical(correct, "tune")) {
      tuneN <- nrow(modglmm[[1]]$Y_I)
      message("Inferred tuneN from NULL model: ", tuneN)
    }
  }
  if (identical(correct, "tune")) {
    if (is.null(tuneN)) {
      stop("'NULLmodelFile' must be supplied when correct='tune' so tuneN can be inferred.")
    }
  }

  escape_regex <- function(x) {
    gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
  }

  normalize_scope <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    x <- trimws(as.character(x))
    if (!nzchar(x) || toupper(x) == "NULL") {
      return(NULL)
    }
    x <- sub("^chr", "", x, ignore.case = TRUE)
    if (!grepl("^([1-9]|1[0-9]|2[0-2])$", x)) {
      stop("'chrom' must be NULL or one of 1..22 (optionally prefixed with 'chr'). Received: ", x)
    }
    paste0("chr", as.integer(x))
  }

  step2_dir <- dirname(inputPrefix)
  step2_base <- sub("_+$", "", basename(inputPrefix))
  if (grepl("_(allchr|chr([1-9]|1[0-9]|2[0-2]))$", step2_base)) {
    stop(
      "'inputPrefix' must be the shared Step2 base prefix without '_allchr' or '_chrN'. ",
      "Use inputPrefix='", sub("_(allchr|chr([1-9]|1[0-9]|2[0-2]))$", "", step2_base),
      "' together with --chrom=NULL or --chrom=1..22."
    )
  }

  requested_scope <- normalize_scope(chrom)
  if (is.null(requested_scope)) {
    requested_scope <- "allchr"
  }

  pattern <- paste0("^", escape_regex(step2_base), "_", escape_regex(requested_scope), "_.*[.]txt$")
  files <- list.files(step2_dir, pattern = pattern, full.names = TRUE)
  files <- files[!grepl("_corrected[.]txt$", basename(files))]
  if (!length(files)) {
    expected <- paste0(step2_base, "_", requested_scope, "_<feature>.txt")
    if (identical(requested_scope, "allchr")) {
      stop(
        "No Step2 files found for correction scope 'allchr' under base prefix: ", inputPrefix,
        ". Expected to see files like ", expected, "."
      )
    }
    stop(
      "No Step2 files found for correction scope '", requested_scope,
      "' under base prefix: ", inputPrefix,
      ". Expected to see files like ", expected, "."
    )
  }

  extract_feature <- function(path) {
    base <- basename(path)
    sub(paste0("^", escape_regex(step2_base), "_", escape_regex(requested_scope), "_(.*)[.]txt$"), "\\1", base)
  }
  features <- vapply(files, extract_feature, character(1))
  if (anyDuplicated(features)) {
    dup <- unique(features[duplicated(features)])
    stop("Duplicated feature names detected for inputPrefix: ", paste(dup, collapse = ", "))
  }
  names(files) <- features

  build_output_path <- function(src, feat) {
    if (isTRUE(overwriteOutput)) {
      return(src)
    }
    file.path(step2_dir, paste0(step2_base, "_corrected_", requested_scope, "_", feat, ".txt"))
  }

  copy_step2_outputs <- function() {
    out_paths <- character(length(files))
    names(out_paths) <- names(files)
    for (feat in names(files)) {
      src <- files[[feat]]
      dst <- build_output_path(src, feat)
      out_dir <- dirname(dst)
      if (!out_dir %in% c("", ".")) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      }
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

  if (identical(correct, "median")) {
    delta <- apply(est_mat, 2, function(x) stats::median(-x, na.rm = TRUE))
  } else {
    palm_tune <- utils::getFromNamespace("palm_tune", "PALM")
    delta <- rep(NA_real_, length(snp_ids))
    names(delta) <- snp_ids
    for (snp in snp_ids) {
      tune_res <- palm_tune(
        summary.stats = list(
          Study1 = list(
            est = est_mat[, snp, drop = FALSE],
            stderr = stderr_mat[, snp, drop = FALSE],
            n = tuneN
          )
        ),
        output.best.one = TRUE,
        verbose = FALSE
      )
      delta[[snp]] <- unname(tune_res[[1]]$delta[[1]])
    }
  }

  adjust_part <- vapply(snp_ids, function(snp) {
    se <- stderr_mat[, snp]
    non_na <- !is.na(se)
    if (!any(non_na)) {
      return(NA_real_)
    }
    sum(non_na) / (2 * sum(1 / sqrt(2 * pi) / se[non_na]))^2
  }, numeric(1))
  names(adjust_part) <- snp_ids

  out_paths <- character(length(feature_ids))
  names(out_paths) <- feature_ids
  for (feat in feature_ids) {
    dat <- feature_data[[feat]]
    idx <- match(dat$SNP, snp_ids)
    ok <- !is.na(dat$est) & !is.na(dat$stderr) & !is.na(delta[idx]) & !is.na(adjust_part[idx])
    dat$est[ok] <- dat$est[ok] + delta[idx][ok]
    dat$stderr[ok] <- sqrt(dat$stderr[ok]^2 + adjust_part[idx][ok])
    dat$pval <- 1 - stats::pchisq((dat$est / dat$stderr)^2, df = 1)

    out_file <- build_output_path(files[[feat]], feat)
    out_dir <- dirname(out_file)
    if (!out_dir %in% c("", ".")) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
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
    "Step2 ", correct, " correction finished for scope ", requested_scope,
    ": wrote ", length(out_paths), " corrected file(s)"
  )
  invisible(unname(out_paths))
}
