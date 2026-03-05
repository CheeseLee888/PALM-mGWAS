#!/usr/bin/env Rscript
# step4_visualization.R
#
# Modes:
#   1) no --pheno and no --snp  -> one big combined plot (significant points across any pheno & snp)
#   2) --pheno only             -> Manhattan for that pheno
#   3) --snp only               -> Forest across phenos for that SNP (by default draw ALL phenos that contain the SNP; optional sig-only)
#   4) --pheno and --snp        -> Forest for that (pheno, snp) across studies
#
# Data format assumed for each meta file (tab-delimited):
#   SNP CHR POS stderr pval pval.het study1_est study1_stderr study2_est ...

# ----------------------------- utilities -----------------------------

#' Print a formatted message with newline
#' @export
msg <- function(...) cat(sprintf(...), "\n")

# Safely read a tabular file with required columns
safe_fread <- function(path, sep = "\t") {
  if (!file.exists(path)) stop("File not found: ", path)
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  df <- data.table::fread(path, sep = sep, data.table = FALSE)
  req <- c("SNP", "CHR", "POS")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop("Missing required columns in ", basename(path), ": ", paste(miss, collapse = ", "))
  df
}

# Normalize meta/step2 columns so downstream plotting works for single-study inputs
#
# If meta_* columns are absent but step2 columns exist, they are duplicated to
# meta_* names. Study columns are added when missing (single study only).
# @param df Data frame returned by `safe_fread()`.
# @param study_id Character scalar used when adding study*_est/stderr columns.
# @return Data frame with meta_* columns present.
standardize_meta_df <- function(df, study_id = "study1") {
  has_meta <- c("meta_pval", "meta_est", "meta_stderr") %in% names(df)
  if (!all(has_meta) && all(c("pval", "est", "stderr") %in% names(df))) {
    df$meta_pval <- df$pval
    df$meta_est <- df$est
    df$meta_stderr <- df$stderr
  }

  has_study_cols <- any(grepl("^study[0-9]+_est$", names(df)))
  if (!has_study_cols) {
    # ensure the synthetic study id matches detect_studies() expectation
    sid <- if (grepl("^study[0-9]+$", study_id)) study_id else "study1"
    fallback <- function(primary, backup) {
      if (!is.null(primary)) return(primary)
      backup
    }
    df[[paste0(sid, "_est")]] <- fallback(df$est, df$meta_est)
    df[[paste0(sid, "_stderr")]] <- fallback(df$stderr, df$meta_stderr)
  }

  df
}

# Read result file and harmonize columns (works for step3 meta or single-study step2)
read_result_file <- function(path, sep = "\t", study_id = NULL) {
  df <- safe_fread(path, sep = sep)

  if (is.null(study_id) || !nzchar(study_id)) {
    study_id <- basename(dirname(path))
    if (!nzchar(study_id)) study_id <- "study1"
  }

  standardize_meta_df(df, study_id = study_id)
}

#' Discover meta-analysis result files
#'
#' Utility to list meta files and extract phenotype names from filenames.
#'
#' @param metaDir Directory containing meta files.
#' @param pattern Regex pattern for filenames; default matches `step3_meta_*.txt`.
#'
#' @return A data frame with columns `pheno` and `file`.
#' @export
discover_meta_files <- function(metaDir,
                                pattern = "step3_meta_.*\\.txt$") {
  files <- list.files(metaDir,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop(
      "No result files found in: ", metaDir,
      "\nLooked for pattern: ", pattern,
      "\nFiles present: ",
      paste(list.files(metaDir), collapse = ", ")
    )
  }

  # Derive phenotype name from the first `.*` in the pattern (e.g., step3_meta_.*\\.txt$ or step2_allchr_.*\\.txt$)
  base_names <- basename(files)
  capture_pattern <- sub("\\.\\*", "(.*)", pattern)
  capture_pattern <- sub("\\*", "(.*)", capture_pattern, fixed = FALSE)
  matches <- regexec(capture_pattern, base_names)
  extracted <- regmatches(base_names, matches)
  pheno <- vapply(extracted, function(x) {
    if (length(x) >= 2) x[2] else NA_character_
  }, character(1))

  # Fallback: strip extension when capture fails
  if (any(is.na(pheno))) {
    pheno[is.na(pheno)] <- tools::file_path_sans_ext(base_names[is.na(pheno)])
  }

  data.frame(
    pheno = pheno,
    file = files,
    mode = "pattern",
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(pheno)
}

# Detect study columns (est/stderr) in a meta result data frame
detect_studies <- function(df) {
  est_cols <- grep("^study[0-9]+_est$", names(df), value = TRUE)
  if (length(est_cols) == 0) stop("Cannot find any columns like study1_est, study2_est, ...")
  ids <- sub("_est$", "", est_cols)
  stderr_cols <- paste0(ids, "_stderr")
  miss <- setdiff(stderr_cols, names(df))
  if (length(miss) > 0) stop("Missing stderr columns: ", paste(miss, collapse = ", "))
  list(ids = ids, est_cols = est_cols, stderr_cols = stderr_cols)
}

# Test if a row is significant at pCut
is_sig_row <- function(row, pCut = 1e-5) {
  if ("meta_pval" %in% names(row) && !is.na(row$meta_pval)) return(row$meta_pval <= pCut)
  FALSE
}

# Parse xlim string like "min,max" into numeric vector
parse_xlim <- function(xlim_str) {
  if (is.null(xlim_str) || is.na(xlim_str) || !nzchar(xlim_str)) {
    return(NULL)
  }
  parts <- stringr::str_split(xlim_str, ",", simplify = TRUE)
  if (ncol(parts) != 2) {
    return(NULL)
  }
  out <- suppressWarnings(as.numeric(parts))
  if (any(is.na(out))) {
    return(NULL)
  }
  out
}

#' Sanitize a string for safe filenames
#'
#' Replaces path separators and other disallowed characters with underscores.
#'
#' @param x Character vector to sanitize.
#' @return Sanitized character vector.
#' @export
sanitize_filename <- function(x) {
  x <- gsub("[:/\\\\]", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-\\.]+", "_", x)
  x
}

# Append a suffix before the filename extension
#
# @param path Path to modify.
# @param suffix Suffix string to insert (e.g., "_qq").
# @return Modified path with suffix inserted before extension.
with_suffix <- function(path, suffix) {
  ext <- tools::file_ext(path)
  base <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", path) else path
  ext_part <- if (nzchar(ext)) paste0(".", ext) else ""
  paste0(base, suffix, ext_part)
}

# Write top-N rows ordered by p-value
#
# @param df Data frame containing column `pval`.
# @param outFile Output path.
# @param n Number of rows to keep.
# @param mode Optional string; when \"step2\" keep step2-style columns.
# @return Invisibly returns the top-N data frame.
write_top_n <- function(df, outFile, n = 10, mode = NULL) {
  if (!("meta_pval" %in% names(df))) stop("Missing meta_pval column for top list.")
  pcol <- "meta_pval"
  top <- df |>
    dplyr::filter(!is.na(.data[[pcol]])) |>
    dplyr::arrange(.data[[pcol]]) |>
    utils::head(n)

  if (identical(mode, "step2")) {
    # match step2 column style
    # ensure pval column exists (use meta_pval fallback)
    if (!("pval" %in% names(top))) top$pval <- top[[pcol]]
    keep <- intersect(c("SNP", "CHR", "POS", "est", "stderr", "pval"), names(top))
    top <- top |> dplyr::select(dplyr::all_of(keep))
  }

  msg("Saving: %s", outFile)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(top, outFile, sep = "\t")
  } else {
    utils::write.table(top, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  invisible(top)
}

# ----------------------------- Manhattan helpers -----------------------------

# Add cumulative base-pair position for Manhattan plotting
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


# Draw a Manhattan plot and save to JPEG
plot_manhattan <- function(df, outFile, title = NULL,
                           pCut = 1e-5,
                           onlySig = FALSE,
                           width = NA_real_, height = NA_real_, dpi = 300) {
  if (!requireNamespace("qqman", quietly = TRUE)) {
    stop("Package 'qqman' is required.")
  }

  # ---- hard-coded header ----
  required_cols <- c("SNP", "CHR", "POS")
  miss <- setdiff(required_cols, names(df))
  if (length(miss) > 0) {
    stop("Missing column(s): ", paste(miss, collapse = ", "))
  }
  pcol <- "meta_pval"
  if (!(pcol %in% names(df))) stop("No p-value column (meta_pval) for Manhattan plot.")

  # ---- prepare ----
  df <- df |>
    dplyr::filter(!is.na(.data[[pcol]]), .data[[pcol]] > 0)

  if (onlySig) {
    df <- dplyr::filter(df, .data[[pcol]] <= pCut)
  }

  if (nrow(df) == 0) {
    stop("No points to plot.")
  }

  man <- df |>
    dplyr::transmute(
      CHR = as.integer(.data$CHR),
      BP  = as.integer(.data$POS),
      P   = as.numeric(.data[[pcol]]),
      SNP = as.character(.data$SNP)
    )

  auto_width <- if (is.na(width)) {
    max(11, length(unique(man$CHR)) * 0.45)
  } else {
    width
  }
  auto_height <- if (is.na(height)) 5 else height

  msg("Saving: %s", outFile)

  grDevices::jpeg(outFile,
    width = auto_width, height = auto_height,
    units = "in", res = dpi
  )

  qqman::manhattan(
    man,
    chr = "CHR",
    bp = "BP",
    p = "P",
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

# Draw a QQ plot and save to file
#
# @param df Data frame containing a p-value column named `pval`.
# @param outFile Output path (png/jpg/etc).
# @param title Optional plot title.
# @param width,height,dpi Device parameters.
plot_qq <- function(df, outFile, title = NULL,
                    width = NA_real_, height = NA_real_, dpi = 300) {
  if (!requireNamespace("qqman", quietly = TRUE)) {
    stop("Package 'qqman' is required.")
  }

  pcol <- "meta_pval"
  if (!(pcol %in% names(df))) stop("Missing p-value column (meta_pval).")

  df <- df |>
    dplyr::filter(!is.na(.data[[pcol]]), .data[[pcol]] > 0)

  if (nrow(df) == 0) stop("No valid p-values to plot.")

  msg("Saving: %s", outFile)

  auto_width <- if (is.na(width)) 6 else width
  auto_height <- if (is.na(height)) 6 else height

  grDevices::jpeg(outFile,
    width = auto_width, height = auto_height,
    units = "in", res = dpi
  )

  qqman::qq(df[[pcol]], main = title %||% "QQ Plot")

  grDevices::dev.off()

  invisible(NULL)
}

# ----------------------------- Forest plot helpers ----------------------------- # nolint: line_length_linter.

# Convert wide study columns to long format for forest plotting
to_long_study <- function(df_rows, studies, y_col, ciMult = 1.96, studyLabels = NULL) {
  # df_rows: data.frame with multiple y categories (pheno or SNP)
  long <- lapply(studies$ids, function(st) {
    tibble::tibble(
      y = df_rows[[y_col]],
      Study = st,
      est = df_rows[[paste0(st, "_est")]],
      se = df_rows[[paste0(st, "_stderr")]]
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

# Make a forest plot and save to file
forest_plot <- function(df_long, y_levels, outFile, title = NULL, xlab = "Effect",
                        xlim_num = NULL, het_y = NULL, width = NA_real_, height = NA_real_, dpi = 300,
                        show_meta = FALSE, df_meta = NULL, show_het = TRUE) {
  df_long <- df_long |> dplyr::mutate(y = factor(y, levels = y_levels))
  study_lvls <- levels(df_long$Study)

  # combine studies + meta as plotting elements so ordering is top→bottom
  meta_ok <- show_meta && length(study_lvls) > 1 && !is.null(df_meta) && nrow(df_meta) > 0
  if (show_meta && length(study_lvls) <= 1) {
    msg("Only one study present; skipping meta overlay.")
  }

  study_lvls_all <- c(study_lvls, if (meta_ok) "Meta" else NULL)

  df_plot <- df_long
  if (meta_ok) {
    df_meta <- df_meta |> dplyr::mutate(y = factor(y, levels = y_levels), Study = "Meta")
    df_plot <- dplyr::bind_rows(df_plot, df_meta)
  }

  df_plot <- df_plot |> dplyr::mutate(Study = factor(Study, levels = study_lvls_all))

  # palette (up to 8 nice colors; recycle if more) + black for meta
  base_cols <- c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")
  pal <- setNames(rep(base_cols, length.out = length(study_lvls_all)), study_lvls_all)
  if (meta_ok) pal["Meta"] <- "#000000"

  # shapes: circles for studies, diamond for meta
  shape_vals <- rep(16, length(study_lvls_all))
  if (meta_ok) shape_vals[length(shape_vals)] <- 18

  # compute vertical offsets so elements draw in order from top to bottom
  spacing <- 0.30
  offsets <- setNames(seq(0, by = -spacing, length.out = length(study_lvls_all)), study_lvls_all)
  df_plot <- df_plot |> dplyr::mutate(y_num = as.numeric(y) + offsets[as.character(Study)])

  g <- ggplot2::ggplot()

  # heterogeneity highlight background (optional)
  if (show_het && !is.null(het_y) && length(het_y) > 0) {
    df_area <- tibble::tibble(y = factor(het_y, levels = y_levels))
    df_area <- df_area |> dplyr::mutate(y_num = as.numeric(y))
    g <- g + ggplot2::geom_rect(
      data = df_area,
      ggplot2::aes(
        ymin = y_num - 0.5,
        ymax = y_num + 0.5,
        xmin = -Inf, xmax = Inf
      ),
      inherit.aes = FALSE,
      fill = "yellow", alpha = 0.18
    )
  }

  g <- g +
    ggplot2::geom_errorbarh(
      data = df_plot,
      ggplot2::aes(y = y_num, xmin = lower, xmax = upper, color = Study),
      height = 0,
      linewidth = 0.55
    ) +
    ggplot2::geom_point(
      data = df_plot,
      ggplot2::aes(y = y_num, x = est, color = Study, shape = Study),
      size = 2.6
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, color = "#990000") +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::scale_shape_manual(values = shape_vals) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(y_levels),
      labels = y_levels,
      expand = ggplot2::expansion(mult = c(0.04, 0.08))
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      x = xlab,
      y = NULL,
      title = title %||% "",
      color = "Study",
      shape = "Study"
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.9),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(10, 12, 10, 20),
      legend.position = "right"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(reverse = TRUE),
      shape = ggplot2::guide_legend(reverse = TRUE)
    )

  if (!is.null(xlim_num)) {
    g <- g + ggplot2::xlim(xlim_num[1], xlim_num[2])
  }

  # auto width/height when not provided: scale height by number of rows, width by studies
  auto_width <- if (is.na(width)) {
    max(7, length(study_lvls_all) * 1.1 + 3)
  } else {
    width
  }
  auto_height <- if (is.na(height)) {
    max(6, length(y_levels) * 0.60 + 2 + if (meta_ok) 0.6 else 0)
  } else {
    height
  }

  msg("Saving: %s", outFile)
  ggplot2::ggsave(outFile, g, width = auto_width, height = auto_height, dpi = dpi)
  invisible(g)
}

# Infix helper: return first non-empty string
`%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

# ----------------------------- Mode implementations -----------------------------

# Mode A: one big combined plot (best pheno per SNP)
#' Create a combined Manhattan plot using, for each SNP, the phenotype with the
#' smallest p-value. All points are shown (after best-per-SNP selection); no
#' significance filtering is applied.
#'
#' @param metaIndex Data frame with columns `pheno` and `file` pointing to meta
#'   result files (as produced by step3).
#' @param outFile Path to the output JPEG.
#' @param sep Field separator for meta files.
#' @param width,height,dpi Graphics device parameters passed to `ggsave`.
#' @param printCut Optional cutoff; SNP/phenotype pairs with best p below this
#'   value are printed to the console and written to `<outFile>_printCut.txt`.
#'
#' @return Invisibly returns the data frame passed to `qqman::manhattan()`.
#' @export
mode_big_combined <- function(metaIndex, outFile,
                              sep = "\t",
                              width = NA_real_, height = NA_real_, dpi = 300,
                              printCut = 1e-8) {
  all_hits <- list()

  for (i in seq_len(nrow(metaIndex))) {
    ph <- metaIndex$pheno[i]
    f <- metaIndex$file[i]
    df <- read_result_file(f, sep = sep, study_id = ph)

    if (!("meta_pval" %in% names(df))) {
      stop("Missing meta_pval column in meta file: ", f)
    }
    pcol <- "meta_pval"

    keep <- df |>
      dplyr::filter(!is.na(.data[[pcol]]), .data[[pcol]] > 0) |>
      dplyr::mutate(pheno = ph) |>
      dplyr::select(
        dplyr::any_of(c("SNP", "CHR", "POS", pcol)),
        .data$pheno
      )

    if (nrow(keep) == 0) next

    all_hits[[length(all_hits) + 1]] <- keep
  }

  if (length(all_hits) == 0) stop("No valid hits found across any pheno (check pval column).")

  big <- dplyr::bind_rows(all_hits)

  # pick the phenotype with the smallest p for each SNP
  best <- big |>
    dplyr::group_by(.data$SNP) |>
    dplyr::arrange(.data$meta_pval, .data$pheno) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()

  # optionally report pairs that beat the user threshold
  if (!is.na(printCut)) {
    hits_to_print <- best |>
      dplyr::filter(!is.na(.data$meta_pval), .data$meta_pval < printCut) |>
      dplyr::arrange(.data$meta_pval)

    # For single-study (step2-only), output columns: SNP CHR POS est stderr pval pheno
    if ("mode" %in% names(metaIndex) && all(metaIndex$mode == "step2")) {
      # ensure pval present
      if (!("pval" %in% names(hits_to_print))) {
        hits_to_print$pval <- hits_to_print$meta_pval
      }
      keep_cols <- c("SNP", "CHR", "POS", "est", "stderr", "pval", "pheno")
      keep_cols <- intersect(keep_cols, names(hits_to_print))
      hits_to_print <- hits_to_print |> dplyr::select(dplyr::all_of(keep_cols))
    }

    if (nrow(hits_to_print) > 0) {
      msg("SNP/pheno pairs with p < %g (best per SNP):", printCut)
      apply(hits_to_print, 1, function(r) {
        p_out <- if ("meta_pval" %in% names(r)) as.numeric(r[["meta_pval"]]) else as.numeric(r[["pval"]])
        msg("  %s\t%s\t%.3e", r[["SNP"]], r[["pheno"]], p_out)
        NULL
      })

      # save alongside plot, force .txt (strip any existing extension)
      base_no_ext <- tools::file_path_sans_ext(outFile)
      list_out <- paste0(base_no_ext, "_printCut.txt")
      msg("Saving list: %s", list_out)
      if (requireNamespace("data.table", quietly = TRUE)) {
        data.table::fwrite(hits_to_print, list_out, sep = "\t")
      } else {
        utils::write.table(hits_to_print, list_out, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    } else {
      msg("No SNP/pheno pairs found with p < %g (best per SNP).", printCut)
    }
  }

  plot_manhattan(
    df = best,
    outFile = outFile,
    title = "Best phenotype per SNP (minimum p across phenotypes)",
    onlySig = FALSE,
    width = width,
    height = height,
    dpi = dpi
  )
}


# Mode B: Manhattan for a given pheno
#' Manhattan plot for a single phenotype
#'
#' @inheritParams mode_big_combined
#' @param pCut P-value cutoff used for Manhattan suggestive line and optional
#'   filtering.
#' @param phenoName Phenotype name matching a row in `metaIndex$pheno`.
#' @param onlySig If TRUE, keep only points passing `pCut`.
#' @param qqOutFile Optional output path for QQ plot.
#' @param topOutFile Optional output path for top-N hits (ordered by p-value).
#' @param top_n How many hits to keep if `topOutFile` is provided.
#'
#' @export
mode_pheno_manhattan <- function(metaIndex, phenoName, outFile, pCut,
                                 sep = "\t", onlySig = FALSE,
                                 width = NA_real_, height = NA_real_, dpi = 300,
                                 qqOutFile = NULL, qq_width = NA_real_, qq_height = NA_real_,
                                 topOutFile = NULL, top_n = 10) {
  row <- metaIndex |> dplyr::filter(.data$pheno == phenoName)
  if (nrow(row) == 0) stop("Cannot find meta file for pheno: ", phenoName)

  df <- read_result_file(row$file[1], sep = sep, study_id = phenoName)
  plot_manhattan(
    df = df,
    outFile = outFile,
    title = paste0("Manhattan: ", phenoName),
    pCut = pCut,
    onlySig = onlySig,
    width = width, height = height, dpi = dpi
  )

  if (!is.null(qqOutFile)) {
    plot_qq(
      df = df,
      outFile = qqOutFile,
      title = paste0("QQ: ", phenoName),
      width = qq_width,
      height = qq_height,
      dpi = dpi
    )
  }

  # if (!is.null(topOutFile)) {
  #   write_top_n(
  #     df = df,
  #     outFile = topOutFile,
  #     n = top_n,
  #     mode = row$mode[1] %||% NULL
  #   )
  # }
}


# Mode C: SNP fixed; forest across phenos
# sigOnlyPheno is OPTIONAL (user request): default FALSE => draw ALL phenos that contain this SNP.
#' Forest plot for a SNP across all phenotypes
#'
#' @inheritParams mode_big_combined
#' @param pCut P-value cutoff used to decide significance across phenotypes
#'   when `sigOnlyPheno` is TRUE.
#' @param snp SNP ID to plot.
#' @param sigOnlyPheno If TRUE, keep only phenotypes where this SNP passes
#'   `pCut`.
#' @param ciMult Multiplier for confidence interval width.
#' @param studyLabels Optional labels replacing study IDs in the legend.
#' @param xlim_str Optional comma-separated numeric limits for the x-axis.
#' @param show_meta If TRUE, overlay meta-effect in black.
#'
#' @export
mode_snp_forest_across_phenos <- function(metaIndex, snp, outFile,
                                          pCut = 1e-5,
                                          sigOnlyPheno = FALSE,
                                          ciMult = 1.96,
                                          studyLabels = NULL,
                                          sep = "\t",
                                          xlim_str = NA_character_,
                                          width = NA_real_, height = NA_real_, dpi = 300,
                                          show_meta = TRUE,
                                          show_het = TRUE) {
  rows <- list()

  for (i in seq_len(nrow(metaIndex))) {
    ph <- metaIndex$pheno[i]
    f <- metaIndex$file[i]
    df <- read_result_file(f, sep = sep, study_id = ph)

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
  if (!("meta_pval" %in% names(df_all))) {
    stop("Missing meta_pval column after binding phenotypes.")
  }
  df_all <- df_all |>
    dplyr::mutate(score = .data$meta_pval) |>
    dplyr::arrange(score)

  studies <- detect_studies(df_all)

  df_long <- to_long_study(df_all, studies, y_col = "pheno", ciMult = ciMult, studyLabels = studyLabels)

  # meta overall
  df_meta <- NULL
  est_col <- "meta_est"
  se_col <- "meta_stderr"
  has_meta_cols <- all(c(est_col, se_col) %in% names(df_all))
  show_meta_flag <- show_meta && has_meta_cols
  if (show_meta && !has_meta_cols) {
    msg("meta_est/meta_stderr not found; skipping meta overlay for SNP forest.")
  }
  if (show_meta_flag) {
    df_meta <- df_all |>
      dplyr::transmute(
        y = pheno,
        est = .data[[est_col]],
        se = .data[[se_col]],
        lower = est - ciMult * .data[[se_col]],
        upper = est + ciMult * .data[[se_col]]
      )
  }

  # highlight phenos with heterogeneity (use pval.het)
  het_y <- character(0)
  het_col <- if ("meta_pval.het" %in% names(df_all)) "meta_pval.het" else NULL
  if (!is.null(het_col)) {
    het_y <- df_all$pheno[!is.na(df_all[[het_col]])]
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
    show_meta = show_meta_flag,
    df_meta = df_meta,
    show_het = show_het
  )
}

# Mode D: pheno + snp fixed; forest across studies
#' Forest plot for a given phenotype/SNP across studies
#'
#' @inheritParams mode_snp_forest_across_phenos
#' @param pheno Phenotype name to plot.
#'
#' @export
mode_pheno_snp_forest <- function(metaIndex, pheno, snp, outFile,
                                  ciMult = 1.96,
                                  studyLabels = NULL,
                                  sep = "\t",
                                  xlim_str = NA_character_,
                                  width = NA_real_, height = NA_real_, dpi = 300,
                                  show_meta = TRUE,
                                  show_het = TRUE) {
  row <- metaIndex |> dplyr::filter(.data$pheno == pheno)
  if (nrow(row) == 0) stop("Cannot find meta file for pheno: ", pheno)

  df <- read_result_file(row$file[1], sep = sep, study_id = pheno)
  r <- df[df$SNP == snp, , drop = FALSE]
  if (nrow(r) == 0) stop("SNP not found in this pheno meta file: ", snp)

  studies <- detect_studies(r)

  r$y <- pheno
  df_long <- to_long_study(r, studies, y_col = "y", ciMult = ciMult, studyLabels = studyLabels)

  df_meta <- NULL
  est_col <- "meta_est"
  se_col <- "meta_stderr"
  has_meta_cols <- all(c(est_col, se_col) %in% names(r))
  show_meta_flag <- show_meta && has_meta_cols
  if (show_meta && !has_meta_cols) {
    msg("meta_est/meta_stderr not found; skipping meta overlay for pheno/SNP forest.")
  }
  if (show_meta_flag) {
    df_meta <- r |>
      dplyr::transmute(
        y = y,
        est = .data[[est_col]],
        se = .data[[se_col]],
        lower = est - ciMult * .data[[se_col]],
        upper = est + ciMult * .data[[se_col]]
      )
  }

  het_y <- character(0)
  het_col <- if ("meta_pval.het" %in% names(r)) "meta_pval.het" else NULL
  if (!is.null(het_col) && !is.na(r[[het_col]][1])) het_y <- pheno

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
    show_meta = show_meta_flag,
    df_meta = df_meta,
    show_het = show_het
  )
}
