#!/usr/bin/env Rscript
# step4_visualization.R
#
# Modes (professor requirement):
#   1) no --pheno and no --snp  -> one big combined plot (significant points across any pheno & snp)
#   2) --pheno only             -> Manhattan for that pheno
#   3) --snp only               -> Forest across phenos for that SNP (by default draw ALL phenos that contain the SNP; optional sig-only)
#   4) --pheno and --snp        -> Forest for that (pheno, snp) across studies
#
# Data format assumed for each meta file (tab-delimited):
#   SNP CHR POS stderr pval pval.het study1_est study1_stderr study2_est ...

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(qqman)
})

# ----------------------------- utilities -----------------------------

msg <- function(...) cat(sprintf(...), "\n")

safe_fread <- function(path, sep = "\t") {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- data.table::fread(path, sep = sep, data.table = FALSE)
  req <- c("SNP", "CHR", "POS")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop("Missing required columns in ", basename(path), ": ", paste(miss, collapse = ", "))
  df
}

discover_meta_files <- function(metaDir, methodPrefix = NULL, pattern = "step3_meta_.*\\.txt$") {
  files <- list.files(metaDir, pattern = pattern, full.names = TRUE)
  if (!is.null(methodPrefix) && nzchar(methodPrefix)) {
    files <- files[grepl(paste0("^", methodPrefix, "_"), basename(files))]
  }
  if (length(files) == 0) stop("No step3 meta files found in: ", metaDir)

  pheno <- sub(".*step3_meta_", "", basename(files))
  pheno <- sub("\\.txt$", "", pheno)

  data.frame(pheno = pheno, file = files, stringsAsFactors = FALSE) |>
    dplyr::arrange(pheno)
}

detect_studies <- function(df) {
  est_cols <- grep("^study[0-9]+_est$", names(df), value = TRUE)
  if (length(est_cols) == 0) stop("Cannot find any columns like study1_est, study2_est, ...")
  ids <- sub("_est$", "", est_cols)
  stderr_cols <- paste0(ids, "_stderr")
  miss <- setdiff(stderr_cols, names(df))
  if (length(miss) > 0) stop("Missing stderr columns: ", paste(miss, collapse = ", "))
  list(ids = ids, est_cols = est_cols, stderr_cols = stderr_cols)
}

pick_sig_field <- function(df) {
  # always use pval
  if ("pval" %in% names(df)) return("pval")
  stop("No pval column found.")
}

is_sig_row <- function(row, pCut = 1e-5) {
  # only use pval for significance
  if ("pval" %in% names(row) && !is.na(row$pval)) return(row$pval <= pCut)
  FALSE
}

parse_xlim <- function(xlim_str) {
  if (is.null(xlim_str) || is.na(xlim_str) || !nzchar(xlim_str)) return(NULL)
  parts <- stringr::str_split(xlim_str, ",", simplify = TRUE)
  if (ncol(parts) != 2) return(NULL)
  out <- suppressWarnings(as.numeric(parts))
  if (any(is.na(out))) return(NULL)
  out
}

sanitize_filename <- function(x) {
  x <- gsub("[:/\\\\]", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-\\.]+", "_", x)
  x
}

# ----------------------------- Manhattan helpers -----------------------------

add_bp_cum <- function(df) {
  df <- df |>
    dplyr::mutate(
      CHR = as.integer(gsub("^chr", "", as.character(.data$CHR), ignore.case = TRUE)),
      POS = as.numeric(.data$POS)
    ) |>
    dplyr::filter(!is.na(.data$CHR), !is.na(.data$POS))

  chr_max <- df |>
    dplyr::group_by(.data$CHR) |>
    dplyr::summarise(chr_len = max(.data$POS, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$CHR) |>
    dplyr::mutate(offset = dplyr::lag(cumsum(.data$chr_len), default = 0))

  df2 <- df |>
    dplyr::left_join(chr_max, by = "CHR") |>
    dplyr::mutate(BPcum = .data$POS + .data$offset)

  axis_df <- chr_max |>
    dplyr::mutate(center = .data$offset + .data$chr_len / 2)

  list(df = df2, axis_df = axis_df)
}


plot_manhattan <- function(df, outFile, title = NULL,
                           pCut = 1e-5,
                           onlySig = FALSE,
                           width = 11, height = 4.5, dpi = 300) {

  if (!requireNamespace("qqman", quietly = TRUE)) {
    stop("Package 'qqman' is required.")
  }

  # ---- hard-coded header ----
  required_cols <- c("SNP", "CHR", "POS", "pval")
  miss <- setdiff(required_cols, names(df))
  if (length(miss) > 0) {
    stop("Missing column(s): ", paste(miss, collapse = ", "))
  }

  # ---- prepare ----
  df <- df |>
    dplyr::filter(!is.na(.data$pval), .data$pval > 0)

  if (onlySig) {
    df <- dplyr::filter(df, .data$pval <= pCut)
  }

  if (nrow(df) == 0)
    stop("No points to plot.")

  man <- df |>
    dplyr::transmute(
      CHR = as.integer(.data$CHR),
      BP  = as.integer(.data$POS),
      P   = as.numeric(.data$pval),
      SNP = as.character(.data$SNP)
    )

  msg("Saving: %s", outFile)

  grDevices::jpeg(outFile, width = width, height = height,
                  units = "in", res = dpi)

  qqman::manhattan(
    man,
    chr = "CHR",
    bp  = "BP",
    p   = "P",
    snp = "SNP",
    logp = TRUE,
    col = c("grey20", "grey70"),
    cex = 0.55,
    main = title %||% "",
    suggestiveline = pCut,
    genomewideline = 5e-8
  )

  grDevices::dev.off()

  invisible(man)
}

# ----------------------------- Forest plot helpers -----------------------------

to_long_study <- function(df_rows, studies, y_col, ciMult = 1.96, studyLabels = NULL) {
  # df_rows: data.frame with multiple y categories (pheno or SNP)
  long <- lapply(studies$ids, function(st) {
    tibble::tibble(
        y = df_rows[[y_col]],
        Study = st,
        est = df_rows[[paste0(st, "_est")]],
        se  = df_rows[[paste0(st, "_stderr")]]
    )
    }) |> dplyr::bind_rows()


  long <- long |> dplyr::mutate(
    lower = est - ciMult * se,
    upper = est + ciMult * se
  )

  # apply study labels if provided
  if (!is.null(studyLabels) && length(studyLabels) == length(studies$ids)) {
    long$Study <- factor(long$Study, levels = studies$ids, labels = studyLabels)
  } else {
    long$Study <- factor(long$Study, levels = studies$ids)
  }

  long
}

forest_plot <- function(df_long, y_levels, outFile, title = NULL, xlab = "Effect",
                        xlim_num = NULL, het_y = NULL, width = 10, height = 7, dpi = 300,
                        show_meta = FALSE, df_meta = NULL) {
  df_long <- df_long |> dplyr::mutate(y = factor(y, levels = y_levels))
  study_lvls <- levels(df_long$Study)

  # palette (up to 8 nice colors; recycle if more)
  base_cols <- c("#F8766D","#CD9600","#7CAE00","#00BE67","#00BFC4","#00A9FF","#C77CFF","#FF61CC")
  pal <- setNames(rep(base_cols, length.out = length(study_lvls)), study_lvls)

  g <- ggplot2::ggplot()

  # heterogeneity highlight background (optional)
  if (!is.null(het_y) && length(het_y) > 0) {
    df_area <- tibble::tibble(y = factor(het_y, levels = y_levels))
    g <- g + ggplot2::geom_rect(
      data = df_area,
      ggplot2::aes(
        xmin = as.numeric(y) - 0.5,
        xmax = as.numeric(y) + 0.5,
        ymin = -Inf, ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "yellow", alpha = 0.18
    )
  }

  g <- g +
    ggplot2::geom_errorbar(
      data = df_long,
      ggplot2::aes(x = y, ymin = lower, ymax = upper, group = Study, color = Study),
      position = ggplot2::position_dodge(width = 0.65),
      width = 0, linewidth = 0.55
    ) +
    ggplot2::geom_point(
      data = df_long,
      ggplot2::aes(x = y, y = est, color = Study),
      position = ggplot2::position_dodge(width = 0.65),
      size = 2.6
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7, color = "#990000") +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      x = NULL,
      y = xlab,
      title = title %||% "",
      color = "Study"
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.9),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(10, 12, 10, 20),
      legend.position = "right"
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(reverse = TRUE))

  # optional meta overlay
  if (show_meta && !is.null(df_meta) && nrow(df_meta) > 0) {
    df_meta <- df_meta |> dplyr::mutate(y = factor(y, levels = y_levels))
    g <- g +
      ggplot2::geom_errorbar(
        data = df_meta,
        ggplot2::aes(x = y, ymin = lower, ymax = upper),
        inherit.aes = FALSE,
        width = 0, linewidth = 0.8, color = "black"
      ) +
      ggplot2::geom_point(
        data = df_meta,
        ggplot2::aes(x = y, y = est),
        inherit.aes = FALSE,
        shape = 18, size = 3.2, color = "black"
      )
  }

  if (!is.null(xlim_num)) {
    g <- g + ggplot2::ylim(xlim_num[1], xlim_num[2])
  }

  msg("Saving: %s", outFile)
  ggplot2::ggsave(outFile, g, width = width, height = height, dpi = dpi)
  invisible(g)
}

`%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

# ----------------------------- Mode implementations -----------------------------

# Mode A: one big combined plot (significant points across any pheno & snp)
#' @export
mode_big_combined <- function(metaIndex, outFile, pCut = 1e-5,
                              maxPoints = 200000, sep = "\t",
                              width = 12, height = 5, dpi = 300) {
  all_hits <- list()

  for (i in seq_len(nrow(metaIndex))) {
    ph <- metaIndex$pheno[i]
    f  <- metaIndex$file[i]
    df <- safe_fread(f, sep = sep)

    # use pval only (qval deprecated)
    if (!("pval" %in% names(df))) next
    df <- df |> dplyr::mutate(pheno = ph) |> dplyr::filter(!is.na(.data$pval), .data$pval > 0, .data$pval <= pCut)

    if (nrow(df) == 0) next

    # keep pval so plot_manhattan can pick sigField
    keep <- df |> dplyr::select(dplyr::any_of(c("SNP", "CHR", "POS", "pval")), .data$pheno)
    all_hits[[length(all_hits) + 1]] <- keep
  }

  if (length(all_hits) == 0) stop("No significant hits found across any pheno (try looser pCut).")

  big <- dplyr::bind_rows(all_hits)

  # downsample by pval
  if (nrow(big) > maxPoints) {
    if ("pval" %in% names(big)) {
      big <- dplyr::arrange(big, .data$pval) |> utils::head(maxPoints)
    } else {
      big <- utils::head(big, maxPoints)
    }
  }

  plot_manhattan(
    df = big,
    outFile = outFile,
    title = "Combined significant hits across phenotypes",
    pCut = pCut,
    onlySig = FALSE,         # already filtered to sig
    width = width,
    height = height,
    dpi = dpi
  )
}


# Mode B: Manhattan for a given pheno
#' @export
mode_pheno_manhattan <- function(metaIndex, phenoName, outFile, pCut,
                                 sep = "\t", onlySig = FALSE,
                                 width = 12, height = 4.5, dpi = 300) {
  row <- metaIndex |> dplyr::filter(.data$pheno == phenoName)
  if (nrow(row) == 0) stop("Cannot find meta file for pheno: ", phenoName)

  df <- safe_fread(row$file[1], sep = sep)

  plot_manhattan(
    df = df,
    outFile = outFile,
    title = paste0("Manhattan: ", phenoName),
    pCut = pCut,
    onlySig = onlySig,
    width = width, height = height, dpi = dpi
  )
}


# Mode C: SNP fixed; forest across phenos
# sigOnlyPheno is OPTIONAL (user request): default FALSE => draw ALL phenos that contain this SNP.
#' @export
mode_snp_forest_across_phenos <- function(metaIndex, snp, outFile,
                                         pCut = 1e-5,
                                         sigOnlyPheno = FALSE,
                                         ciMult = 1.96,
                                         studyLabels = NULL,
                                         sep = "\t",
                                         xlim_str = NA_character_,
                                         width = 10, height = 7, dpi = 300,
                                         show_meta = TRUE) {
  rows <- list()

  for (i in seq_len(nrow(metaIndex))) {
    ph <- metaIndex$pheno[i]
    f  <- metaIndex$file[i]
    df <- safe_fread(f, sep = sep)

    if (!("SNP" %in% names(df))) next
    r <- df[df$SNP == snp, , drop = FALSE]
    if (nrow(r) == 0) next

    # optionally filter by significance across phenos
    if (sigOnlyPheno) {
      if (!is_sig_row(r[1, , drop = FALSE], pCut = pCut)) next
    }

    r$pheno <- ph
    rows[[length(rows) + 1]] <- r
  }

  if (length(rows) == 0) {
    if (sigOnlyPheno) {
      stop("No phenotypes found (with this SNP) passing significance cutoff. Try --sigOnlyPheno FALSE or loosen pCut.")
    } else {
      stop("No phenotypes contain this SNP in meta files: ", snp)
    }
  }

  df_all <- dplyr::bind_rows(rows)

  # sort phenos by significance (use pval)
  sigField <- pick_sig_field(df_all)
  df_all <- df_all |> dplyr::mutate(score = .data[[sigField]])
  df_all <- df_all |> dplyr::arrange(score)

  studies <- detect_studies(df_all)

  df_long <- to_long_study(df_all, studies, y_col = "pheno", ciMult = ciMult, studyLabels = studyLabels)

  # meta overall
  df_meta <- NULL
  if (show_meta && all(c("est", "stderr") %in% names(df_all))) {
    df_meta <- df_all |>
      dplyr::transmute(
        y = pheno,
        est = est,
        se  = stderr,
        lower = est - ciMult * stderr,
        upper = est + ciMult * stderr
      )
  }

  # highlight phenos with heterogeneity (use pval.het)
  het_y <- character(0)
  if ("pval.het" %in% names(df_all)) {
    het_y <- df_all$pheno[!is.na(df_all$pval.het)]
  }

  y_levels <- df_all$pheno

  xlim_num <- parse_xlim(xlim_str)

  forest_plot(
    df_long = df_long,
    y_levels = y_levels,
    outFile = outFile,
    title = paste0("SNP forest across phenotypes: ", snp),
    xlab = "Effect (study-specific; meta in black)",
    xlim_num = xlim_num,
    het_y = het_y,
    width = width, height = height, dpi = dpi,
    show_meta = show_meta,
    df_meta = df_meta
  )
}

# Mode D: pheno + snp fixed; forest across studies
#' @export
mode_pheno_snp_forest <- function(metaIndex, pheno, snp, outFile,
                                  ciMult = 1.96,
                                  studyLabels = NULL,
                                  sep = "\t",
                                  xlim_str = NA_character_,
                                  width = 9, height = 4, dpi = 300,
                                  show_meta = TRUE) {
  row <- metaIndex |> dplyr::filter(.data$pheno == pheno)
  if (nrow(row) == 0) stop("Cannot find meta file for pheno: ", pheno)

  df <- safe_fread(row$file[1], sep = sep)
  r <- df[df$SNP == snp, , drop = FALSE]
  if (nrow(r) == 0) stop("SNP not found in this pheno meta file: ", snp)

  studies <- detect_studies(r)

  r$y <- pheno
  df_long <- to_long_study(r, studies, y_col = "y", ciMult = ciMult, studyLabels = studyLabels)

  df_meta <- NULL
  if (show_meta && all(c("est", "stderr") %in% names(r))) {
    df_meta <- r |>
      dplyr::transmute(
        y = y,
        est = est,
        se  = stderr,
        lower = est - ciMult * stderr,
        upper = est + ciMult * stderr
      )
  }

  het_y <- character(0)
  if ("pval.het" %in% names(r)) {
    if (!is.na(r$pval.het[1])) het_y <- pheno
  }

  xlim_num <- parse_xlim(xlim_str)

  forest_plot(
    df_long = df_long,
    y_levels = c(pheno),
    outFile = outFile,
    title = paste0("Forest: ", pheno, " @ ", snp),
    xlab = "Effect (study-specific; meta in black)",
    xlim_num = xlim_num,
    het_y = het_y,
    width = width, height = height, dpi = dpi,
    show_meta = show_meta,
    df_meta = df_meta
  )
}