#' Meta-analyze step2 results across studies (meta columns prefixed with `meta_`)
#'
#' @param study_dirs Named character vector/list. Names are study IDs, values are directories.
#' @param in_prefix Input file prefix in each study dir, e.g. "palm1_step2_allchr_"
#' @param in_suffix Input file suffix, default ".txt"
#' @param features Optional feature names (without prefix/suffix). If NULL, infer from first study dir.
#' @param out_dir If not NULL, write per-feature meta files to this directory.
#' @param out_prefix Output meta file prefix, e.g. "meta_step2_allchr_"
#' @param out_suffix Output file suffix, default ".txt"
#' @param keep_het If TRUE and multi-study, keep pval.het column; if FALSE, drop it to match 6-column step2 format exactly.
#' @param meta.method (deprecated) no longer used; meta-analysis now uses fixed-effect inverse-variance weighting.
#'
#' @return Named list: each element is a data.frame/tibble in step2 format.
#' @import dplyr
#' @export
metaSummary <- function(study_dirs,
                        in_prefix,
                        in_suffix = ".txt",
                        features = NULL,
                        out_dir = NULL,
                        out_prefix = "meta_step2_allchr_",
                        out_suffix = ".txt",
                        keep_het = TRUE,
                        meta.method = NULL) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }

  if (is.null(names(study_dirs)) || any(names(study_dirs) == "")) {
    stop("study_dirs must be a named vector/list: names are study IDs.")
  }
  study.ID <- names(study_dirs)

  # info: how many studies were provided
  message(sprintf(
    "metaSummary: reading %d study(ies): %s", length(study.ID),
    paste(study.ID, collapse = ", ")
  ))
  message("metaSummary: study directories: ",
    paste(sprintf("%s -> %s", study.ID, study_dirs), collapse = "; ")
  )

  # ensure study directories exist
  missing_dir <- study_dirs[!dir.exists(study_dirs)]
  if (length(missing_dir) > 0) {
    stop(sprintf(
      "metaSummary: %d study directory(ies) not found: %s",
      length(missing_dir),
      paste(names(missing_dir), missing_dir, sep = " -> ", collapse = "; ")
    ))
  }

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # feature (phenotype) availability summary across studies + file matching echo
  feature_scan <- lapply(names(study_dirs), function(sid) {
    d <- study_dirs[[sid]]
    ff <- list.files(
      d,
      pattern = paste0("^", in_prefix, ".*", gsub("\\.", "\\\\.", in_suffix), "$"),
      full.names = FALSE
    )
    feats <- sub(paste0("^", in_prefix), "", ff)
    feats <- sub(paste0(gsub("\\.", "\\\\.", in_suffix), "$"), "", feats)
    feats <- unique(feats)

    # verbose but concise logging: matched files and extracted features
    if (length(ff) == 0) {
      message("metaSummary: study ", sid, " matched 0 files with pattern ", in_prefix, "*", in_suffix)
    } else {
      show_files <- if (length(ff) > 6) c(ff[1:6], "...") else ff
      show_feats <- if (length(feats) > 10) c(feats[1:10], "...") else feats
      message("metaSummary: study ", sid, " matched files: ", paste(show_files, collapse = ", "))
      message("metaSummary: study ", sid, " extracted features: ", paste(show_feats, collapse = ", "))
    }

    list(files = ff, features = feats)
  })

  feature_lists <- lapply(feature_scan, `[[`, "features")
  names(feature_lists) <- names(study_dirs)

  feat_union <- sort(unique(unlist(feature_lists)))
  feat_inter <- if (length(feature_lists) > 1) Reduce(intersect, feature_lists) else feat_union

  # if user does not specify, use union
  if (is.null(features)) {
    features <- feat_union
  }

  considered_feats <- features
  miss_counts <- vapply(feature_lists, function(x) sum(!considered_feats %in% x), integer(1))

  feature_counts <- vapply(feature_lists, length, integer(1))
  message("metaSummary: per-study feature counts: ",
    paste(sprintf("%s=%d", names(feature_counts), feature_counts), collapse = "; ")
  )
  if (length(feature_lists) > 1) {
    message(sprintf("metaSummary: feature intersection size=%d", length(feat_inter)))
  }
  message("metaSummary: per-study missing features (relative to requested set): ",
    paste(sprintf("%s missing=%d", names(miss_counts), miss_counts), collapse = "; ")
  )

  # existence check for logging only (do not abort; missing files => study skipped for that feature)
  expected <- expand.grid(study = names(study_dirs), feature = considered_feats, stringsAsFactors = FALSE)
  expected$path <- file.path(study_dirs[expected$study], paste0(in_prefix, expected$feature, in_suffix))
  missing_mask <- !file.exists(expected$path)
  if (any(missing_mask)) {
    miss <- expected[missing_mask, , drop = FALSE]
    msg_lines <- paste0(miss$study, ":", miss$feature)
    message(sprintf(
      "metaSummary: %d missing input file(s); those study-feature pairs will be skipped. Examples: %s",
      nrow(miss), paste(utils::head(msg_lines, 10), collapse = "; ")
    ))
  } else {
    message("metaSummary: all requested input files exist across studies.")
  }

  .read_step2 <- function(path) {
    if (!file.exists(path)) {
      return(NULL)
    }
    dat <- tryCatch(
      read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = ""),
      error = function(e) NULL
    )
    if (is.null(dat) || nrow(dat) == 0) {
      return(NULL)
    }

    need <- c("SNP", "CHR", "POS", "est", "stderr", "pval")
    miss <- setdiff(need, colnames(dat))
    if (length(miss) > 0) {
      stop(
        "File missing required columns: ", paste(miss, collapse = ", "),
        "\nFile: ", path
      )
    }
    dat <- dat[, need, drop = FALSE]

    dat$CHR <- suppressWarnings(as.integer(dat$CHR))
    dat$POS <- suppressWarnings(as.integer(dat$POS))
    dat$est <- suppressWarnings(as.numeric(dat$est))
    dat$stderr <- suppressWarnings(as.numeric(dat$stderr))
    dat$pval <- suppressWarnings(as.numeric(dat$pval))

    dat <- dat[!duplicated(dat$SNP), , drop = FALSE]
    dat
  }

  .meta_one_feature <- function(feat) {
    per_study <- setNames(vector("list", length(study.ID)), study.ID)
    for (d in study.ID) {
      fpath <- file.path(study_dirs[[d]], paste0(in_prefix, feat, in_suffix))
      per_study[[d]] <- .read_step2(fpath)
    }

    has <- vapply(per_study, function(x) !is.null(x) && nrow(x) > 0, logical(1))
    if (!any(has)) {
      return(NULL)
    }

    used_studies <- names(per_study)[has]
    per_study <- per_study[used_studies]
    # info: how many studies contribute to this feature
    n_used <- length(used_studies)
    message(sprintf(
      "metaSummary: feature '%s' - using %d study(ies): %s",
      feat, n_used, paste(used_studies, collapse = ", ")
    ))

    snp.ID <- unique(unlist(lapply(per_study, `[[`, "SNP")))
    snp.ID <- as.character(snp.ID)

    AA.est <- matrix(NA_real_,
      nrow = length(snp.ID), ncol = length(used_studies),
      dimnames = list(snp.ID, used_studies)
    )
    AA.var <- matrix(NA_real_,
      nrow = length(snp.ID), ncol = length(used_studies),
      dimnames = list(snp.ID, used_studies)
    )

    # CHR/POS: take first observed
    CHR <- rep(NA_integer_, length(snp.ID))
    names(CHR) <- snp.ID
    POS <- rep(NA_integer_, length(snp.ID))
    names(POS) <- snp.ID

    for (d in used_studies) {
      dat <- per_study[[d]]
      idx <- match(dat$SNP, snp.ID)

      AA.est[idx, d] <- dat$est
      AA.var[idx, d] <- (dat$stderr)^2

      miss_chr <- is.na(CHR[idx]) & !is.na(dat$CHR)
      if (any(miss_chr)) CHR[idx[miss_chr]] <- dat$CHR[miss_chr]
      miss_pos <- is.na(POS[idx]) & !is.na(dat$POS)
      if (any(miss_pos)) POS[idx[miss_pos]] <- dat$POS[miss_pos]
    }

    # # progress info: how many SNPs overlap
    # if (length(used_studies) > 1) {
    #   union_n <- length(snp.ID)
    #   per_snp_non_na <- rowSums(!is.na(AA.est))
    #   inter_n <- sum(per_snp_non_na == length(used_studies))
    #   dropped_if_intersect <- union_n - inter_n

    #   # per-study missing counts
    #   per_study_missing <- colSums(is.na(AA.est))
    #   per_study_present <- colSums(!is.na(AA.est))
    #   msg1 <- sprintf(
    #     "metaSummary: feature '%s' SNP union=%d, intersection=%d, would-drop-if-intersect=%d",
    #     feat, union_n, inter_n, dropped_if_intersect
    #   )
    #   msg2 <- paste(sprintf("%s missing=%d present=%d", names(per_study_missing), per_study_missing, per_study_present), collapse = "; ")
    #   message(msg1)
    #   message("metaSummary: per-study SNP counts: ", msg2)
    # }

    if (length(used_studies) > 1) {
      # fixed-effect inverse-variance weighting in matrix form (whole-column ops)
      w <- 1 / AA.var
      # rows with all NA will have wsum=0; guard to avoid inf
      wsum <- rowSums(w, na.rm = TRUE)
      wsum[wsum == 0] <- NA_real_

      meta_est <- rowSums(w * AA.est, na.rm = TRUE) / wsum
      meta_stderr <- sqrt(1 / wsum)
      z <- meta_est / meta_stderr
      meta_pval <- 2 * stats::pnorm(-abs(z))

      # heterogeneity Q (still vectorized)
      centered <- AA.est - meta_est
      Q <- rowSums(w * centered * centered, na.rm = TRUE)
      df <- rowSums(!is.na(AA.est)) - 1
      meta_pval.het <- stats::pchisq(Q, df = df, lower.tail = FALSE)

      meta_fits <- data.frame(
        est = meta_est,
        stderr = meta_stderr,
        pval = meta_pval,
        pval.het = meta_pval.het,
        stringsAsFactors = FALSE
      )

      out <- dplyr::tibble(
        SNP = snp.ID,
        CHR = unname(CHR[snp.ID]),
        POS = unname(POS[snp.ID]),
        meta_est = meta_fits$est,
        meta_stderr = meta_fits$stderr,
        meta_pval = meta_fits$pval
      )

      if (isTRUE(keep_het)) {
        out$meta_pval.het <- meta_fits$`pval.het`
      }

      # add per-study est/stderr columns (SNP-aligned)
      for (d in used_studies) {
        dat <- per_study[[d]]
        m2 <- match(out$SNP, dat$SNP)
        out[[paste0(d, "_est")]] <- dat$est[m2]
        out[[paste0(d, "_stderr")]] <- dat$stderr[m2]
      }

      # optional: keep SNP ordering stable (by CHR/POS if available)
      if (all(!is.na(out$CHR)) && all(!is.na(out$POS))) {
        out <- dplyr::arrange(out, CHR, POS)
      }

      return(out)
    } else {
      d <- used_studies[1]
      dat <- per_study[[d]]

      beta.coef <- dat$est
      std.coef <- dat$stderr
      pval <- 1 - pchisq((beta.coef / std.coef)^2, df = 1)

      out <- dplyr::tibble(
        SNP = dat$SNP,
        CHR = dat$CHR,
        POS = dat$POS,
        meta_est = beta.coef,
        meta_stderr = std.coef,
        meta_pval = pval
      )

      if (all(!is.na(out$CHR)) && all(!is.na(out$POS))) {
        out <- dplyr::arrange(out, CHR, POS)
      }
      return(out)
    }
  }

  res <- setNames(vector("list", length(features)), features)
  for (feat in features) {
    out <- .meta_one_feature(feat)
    res[[feat]] <- out

    # write to disk if requested
    if (!is.null(out_dir) && !is.null(out)) {
      out_path <- file.path(out_dir, paste0(out_prefix, feat, out_suffix))
      write.table(out,
        file = out_path, sep = "\t",
        quote = FALSE, row.names = FALSE, col.names = TRUE
      )
    }
  }

  res <- res[!vapply(res, is.null, logical(1))]
  return(invisible(res))
}
