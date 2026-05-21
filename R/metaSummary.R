#' Meta-analyze step2 results across studies (meta columns prefixed with `meta_`)
#'
#' @param study_dirs Named character vector/list. Names are study IDs, values are directories.
#' @param inputPrefix Shared Step2 base prefix, or a named vector/list of
#'   study-specific Step2 prefixes. Files are expected at
#'   `<prefix>_chr<chrom>_<feature>.txt` when `chrom` is set, or
#'   `<prefix>_allchr_<feature>.txt` when `chrom` is `NULL`. An optional
#'   trailing underscore is ignored.
#' @param chrom Optional chromosome selector. Use `NULL` to meta-analyze
#'   `_allchr` files. Use `1`..`22` or strings like `"chr1"` to
#'   meta-analyze one chromosome-specific shard across all studies.
#' @param featureList Optional feature names (without prefix/suffix). If
#'   `NULL`, infer all features from the selected Step2 scope.
#' @param out_dir If not NULL, write per-feature meta files to this directory.
#' @param out_prefix Output meta file prefix, e.g. "step3_meta". A trailing underscore is ignored.
#' @param out_suffix Output file suffix, default ".txt"
#' @param keep_het If TRUE and multi-study, keep pval.het column; if FALSE, drop it to match 6-column step2 format exactly.
#' @param meta.method Meta-analysis method passed to `metafor::rma.uni()`.
#'   Defaults to `"EE"`.
#'
#' @return Named list: each element is a data.frame/tibble in step2 format.
#' @import dplyr
#' @export
metaSummary <- function(study_dirs,
                        inputPrefix,
                        chrom = NULL,
                        featureList = NULL,
                        out_dir = NULL,
                        out_prefix = "step3_meta",
                        out_suffix = ".txt",
                        keep_het = TRUE,
                        meta.method = "EE") {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (is.null(names(study_dirs)) || any(names(study_dirs) == "")) {
    stop("study_dirs must be a named vector/list: names are study IDs.")
  }
  if (missing(inputPrefix) || is.null(inputPrefix) || length(inputPrefix) == 0L) {
    stop("'inputPrefix' must be provided.")
  }
  if (!is.null(chrom) && length(chrom) != 1L) {
    stop("'chrom' must be NULL or a single chromosome value.")
  }
  if (is.null(meta.method) || length(meta.method) != 1L || !nzchar(trimws(meta.method))) {
    stop("'meta.method' must be a single non-empty method name for metafor::rma.uni().")
  }
  meta.method <- trimws(meta.method)
  fast_fixed_effect <- identical(toupper(meta.method), "EE")
  if (!fast_fixed_effect && !requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required for meta.method='", meta.method, "' but is not installed.")
  }
  study.ID <- names(study_dirs)
  if (length(inputPrefix) == 1L) {
    input_prefixes <- stats::setNames(file.path(as.character(study_dirs), basename(as.character(inputPrefix))), study.ID)
  } else {
    input_prefixes <- stats::setNames(as.character(inputPrefix), names(inputPrefix))
    if (is.null(names(input_prefixes)) || any(names(input_prefixes) == "")) {
      stop("'inputPrefix' must be a single shared prefix or a named vector/list with study IDs as names.")
    }
    missing_prefix <- setdiff(study.ID, names(input_prefixes))
    if (length(missing_prefix) > 0L) {
      stop("'inputPrefix' is missing prefix value(s) for study ID(s): ", paste(missing_prefix, collapse = ", "))
    }
    input_prefixes <- input_prefixes[study.ID]
  }
  input_prefixes <- sub("_+$", "", input_prefixes)
  if (any(is.na(input_prefixes) | !nzchar(input_prefixes))) {
    stop("'inputPrefix' contains empty prefix value(s).")
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

  step2_base <- basename(input_prefixes)
  bad_scope_prefix <- grepl("_(allchr|chr([1-9]|1[0-9]|2[0-2]))$", step2_base)
  if (any(bad_scope_prefix)) {
    stop(
      "'inputPrefix' must be the shared Step2 base prefix without '_allchr' or '_chrN'. ",
      "Use inputPrefix='", sub("_(allchr|chr([1-9]|1[0-9]|2[0-2]))$", "", step2_base[bad_scope_prefix][[1L]]),
      "' together with --chrom=NULL or --chrom=1..22."
    )
  }
  requested_scope <- normalize_scope(chrom)
  if (is.null(requested_scope)) {
    requested_scope <- "allchr"
  }
  file_pattern_by_study <- stats::setNames(
    paste0("^", vapply(step2_base, escape_regex, character(1)), "_", escape_regex(requested_scope), "_(.*)[.]txt$"),
    study.ID
  )

  # info: how many studies were provided
  message(sprintf(
    "metaSummary: reading %d study(ies): %s", length(study.ID),
    paste(study.ID, collapse = ", ")
  ))
  message("metaSummary: input prefixes: ",
    paste(sprintf("%s -> %s", study.ID, input_prefixes), collapse = "; ")
  )

  study_dirs <- stats::setNames(dirname(input_prefixes), study.ID)
  study_dirs[study_dirs %in% c("", ".")] <- "."

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
  out_prefix <- sub("_+$", "", out_prefix)

  # high-level run context (output directory printed near end)
  message("metaSummary: meta method = ", meta.method)
  message("metaSummary: chromosome scope = ", requested_scope)

  extract_feature_names <- function(files, sid) {
    base_names <- basename(files)
    sub(
      file_pattern_by_study[[sid]],
      "\\1",
      base_names
    )
  }

  # feature (phenotype) availability summary across studies + file matching echo
  feature_scan <- lapply(names(study_dirs), function(sid) {
    d <- study_dirs[[sid]]
    ff <- list.files(
      d,
      pattern = file_pattern_by_study[[sid]],
      full.names = FALSE
    )
    feats <- extract_feature_names(ff, sid)

    dup_feats <- unique(feats[duplicated(feats)])
    if (length(dup_feats) > 0) {
      stop(
        "metaSummary: duplicated extracted feature names in study ", sid, ": ",
        paste(dup_feats, collapse = ", "),
        ". The Step2 scope should contain exactly one file per feature."
      )
    }

    file_map <- stats::setNames(ff, feats)

    # logging: only counts, suppress listing file names or phenotype names
    if (length(ff) == 0) {
      message("metaSummary: study ", sid, " matched 0 files for scope ", requested_scope)
    } else {
      message("metaSummary: study ", sid, " matched ", length(ff), " file(s) for scope ", requested_scope)
      message("metaSummary: study ", sid, " extracted ", length(unique(feats)), " feature(s)")
    }

    list(files = ff, features = unique(feats), file_map = file_map)
  })
  names(feature_scan) <- names(study_dirs)

  feature_lists <- lapply(feature_scan, `[[`, "features")
  names(feature_lists) <- names(study_dirs)

  feat_union <- sort(unique(unlist(feature_lists)))
  feat_inter <- if (length(feature_lists) > 1) Reduce(intersect, feature_lists) else feat_union

  if (length(feat_union) == 0L) {
    expected <- paste0(step2_base[[1L]], "_", requested_scope, "_<feature>.txt")
    if (identical(requested_scope, "allchr")) {
      stop(
        "No Step2 files found for meta-analysis scope 'allchr' under the provided input prefix(es). ",
        "Expected to see files like ", expected,
        ". If only chromosome-split files exist, run Step3 once per chromosome with --chrom=1..22."
      )
    }
    stop(
      "No Step2 files found for meta-analysis scope '", requested_scope,
      "' under the provided input prefix(es)",
      ". Expected to see files like ", expected, "."
    )
  }

  if (is.null(featureList)) {
    features <- feat_union
  } else {
    features <- unique(as.character(featureList))
    features <- trimws(features)
    features <- features[nzchar(features)]
    if (length(features) == 0L) {
      stop("'featureList' must be NULL or contain at least one feature name.")
    }
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
  expected$path <- mapply(function(study, feature) {
    fmap <- feature_scan[[study]][["file_map"]]
    if (is.null(fmap) || !(feature %in% names(fmap))) {
      return(NA_character_)
    }
    rel <- unname(fmap[[feature]])
    if (is.null(rel) || is.na(rel) || !nzchar(rel)) return(NA_character_)
    file.path(study_dirs[[study]], rel)
  }, expected$study, expected$feature, USE.NAMES = FALSE)
  missing_mask <- is.na(expected$path) | !file.exists(expected$path)
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
    if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
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
      fmap <- feature_scan[[d]][["file_map"]]
      if (is.null(fmap) || !(feat %in% names(fmap))) {
        rel <- NA_character_
      } else {
        rel <- unname(fmap[[feat]])
      }
      fpath <- if (is.null(rel) || is.na(rel) || !nzchar(rel)) NA_character_ else file.path(study_dirs[[d]], rel)
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
      if (fast_fixed_effect) {
        valid <- !is.na(AA.est) & !is.na(AA.var) & AA.var > 0
        weights <- matrix(0, nrow = nrow(AA.var), ncol = ncol(AA.var), dimnames = dimnames(AA.var))
        weights[valid] <- 1 / AA.var[valid]

        est_values <- AA.est
        est_values[!valid] <- 0

        sum_w <- rowSums(weights)
        meta_est <- rep(NA_real_, length(snp.ID))
        meta_stderr <- rep(NA_real_, length(snp.ID))
        meta_pval <- rep(NA_real_, length(snp.ID))
        meta_pval_het <- rep(NA_real_, length(snp.ID))

        has_weight <- sum_w > 0
        meta_est[has_weight] <- rowSums(weights * est_values)[has_weight] / sum_w[has_weight]
        meta_stderr[has_weight] <- sqrt(1 / sum_w[has_weight])
        meta_pval[has_weight] <- stats::pchisq(
          (meta_est[has_weight] / meta_stderr[has_weight])^2,
          df = 1,
          lower.tail = FALSE
        )

        n_per_snp <- rowSums(valid)
        centered <- sweep(AA.est, 1, meta_est, "-")
        centered[!valid] <- 0
        q_stat <- rowSums(weights * centered^2)
        has_het <- n_per_snp > 1
        meta_pval_het[has_het] <- stats::pchisq(q_stat[has_het], df = n_per_snp[has_het] - 1, lower.tail = FALSE)

        meta_fits <- data.frame(
          est = meta_est,
          stderr = meta_stderr,
          pval = meta_pval,
          pval.het = meta_pval_het,
          stringsAsFactors = FALSE
        )
      } else {
        meta_fits <- data.frame(
          est = rep(NA_real_, length(snp.ID)),
          stderr = rep(NA_real_, length(snp.ID)),
          pval = rep(NA_real_, length(snp.ID)),
          pval.het = rep(NA_real_, length(snp.ID)),
          stringsAsFactors = FALSE
        )
        for (i in seq_along(snp.ID)) {
          keep <- !is.na(AA.est[i, ]) & !is.na(AA.var[i, ]) & AA.var[i, ] > 0
          if (!any(keep)) {
            next
          }
          if (sum(keep) == 1L) {
            beta.coef <- AA.est[i, keep][[1L]]
            std.coef <- sqrt(AA.var[i, keep][[1L]])
            meta_fits$est[i] <- beta.coef
            meta_fits$stderr[i] <- std.coef
            meta_fits$pval[i] <- stats::pchisq((beta.coef / std.coef)^2, df = 1, lower.tail = FALSE)
            next
          }
          fit <- tryCatch(
            metafor::rma.uni(
              yi = as.numeric(AA.est[i, keep]),
              vi = as.numeric(AA.var[i, keep]),
              method = meta.method
            ),
            error = function(e) {
              stop(
                "metafor::rma.uni() failed for feature '", feat,
                "', SNP '", snp.ID[[i]], "' with meta.method='", meta.method,
                "': ", conditionMessage(e),
                call. = FALSE
              )
            }
          )
          meta_fits$est[i] <- as.numeric(fit$b)
          meta_fits$stderr[i] <- fit$se
          meta_fits$pval[i] <- fit$pval
          meta_fits$pval.het[i] <- fit$QEp
        }
      }

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
        out <- dplyr::arrange(out, .data$CHR, .data$POS)
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
        out <- dplyr::arrange(out, .data$CHR, .data$POS)
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
      out_name <- if (nzchar(out_prefix)) {
        paste0(out_prefix, "_", requested_scope, "_", feat, out_suffix)
      } else {
        paste0(requested_scope, "_", feat, out_suffix)
      }
      out_path <- file.path(out_dir, out_name)
      write.table(out,
        file = out_path, sep = "\t",
        quote = FALSE, row.names = FALSE, col.names = TRUE
      )
    }
  }

  res <- res[!vapply(res, is.null, logical(1))]
  # output directory echoed at the end of processing
  message("metaSummary: output directory = ", if (is.null(out_dir)) "NULL" else out_dir)
  return(invisible(res))
}
