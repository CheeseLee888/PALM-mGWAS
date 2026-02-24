#' Meta-analyze step2 results across studies (output in step2 file format)
#'
#' @param study_dirs Named character vector/list. Names are study IDs, values are directories.
#' @param in_prefix Input file prefix in each study dir, e.g. "palm1_step2_allchr_"
#' @param in_suffix Input file suffix, default ".txt"
#' @param features Optional feature names (without prefix/suffix). If NULL, infer from first study dir.
#' @param out_dir If not NULL, write per-feature meta files to this directory.
#' @param out_prefix Output meta file prefix, e.g. "meta_step2_allchr_"
#' @param out_suffix Output file suffix, default ".txt"
#' @param keep_het If TRUE and multi-study, keep pval.het column; if FALSE, drop it to match 6-column step2 format exactly.
#' @param meta.method Passed to metafor::rma(method=...), consistent with your old code.
#'
#' @return Named list: each element is a data.frame/tibble in step2 format.
#' @import dplyr metafor
#' @export
metaSummary <- function(study_dirs,
                        in_prefix,
                        in_suffix = ".txt",
                        features = NULL,
                        out_dir = NULL,
                        out_prefix = "meta_step2_allchr_",
                        out_suffix = ".txt",
                        keep_het = TRUE,
                        meta.method = "EE") {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required but not installed.")
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

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # infer features from first study dir
  if (is.null(features)) {
    d0 <- study_dirs[[1]]
    ff <- list.files(
      d0,
      pattern = paste0("^", in_prefix, ".*", gsub("\\.", "\\\\.", in_suffix), "$"),
      full.names = FALSE
    )
    features <- sub(paste0("^", in_prefix), "", ff)
    features <- sub(paste0(gsub("\\.", "\\\\.", in_suffix), "$"), "", features)
    features <- unique(features)
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

    if (length(used_studies) > 1) {
      # ---- CORE CALC (kept consistent with your old code) ----
      meta_fits <- sapply(seq_len(nrow(AA.est)), function(i) {
        non.id <- !is.na(AA.est[i, ])
        m <- try(
          metafor::rma(yi = AA.est[i, non.id], vi = AA.var[i, non.id], method = meta.method),
          silent = TRUE
        )
        if (class(m)[1] != "try-error") {
          return(c(est = m$beta, stderr = m$se, pval = m$QMp, pval.het = m$QEp))
        } else {
          return(c(est = NA, stderr = NA, pval = NA, pval.het = NA))
        }
      })
      meta_fits <- data.frame(t(meta_fits), stringsAsFactors = FALSE)
      # -------------------------------------------------------

      out <- dplyr::tibble(
        SNP = snp.ID,
        CHR = unname(CHR[snp.ID]),
        POS = unname(POS[snp.ID]),
        est = meta_fits$est,
        stderr = meta_fits$stderr,
        pval = meta_fits$pval
      )

      if (isTRUE(keep_het)) {
        out$pval.het <- meta_fits$`pval.het`
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
        est = beta.coef,
        stderr = std.coef,
        pval = pval
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
  return(res)
}
