#!/usr/bin/env Rscript
# Step4 visualization helpers used by extdata/step4_visualization.R
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

  est_cols <- setdiff(grep("_est$", names(df), value = TRUE), "meta_est")
  stderr_cols <- setdiff(grep("_stderr$", names(df), value = TRUE), "meta_stderr")
  has_study_cols <- length(est_cols) > 0 && length(stderr_cols) > 0
  if (!has_study_cols) {
    sid <- "single_study"
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
  message("discover_meta_files: scanning ", metaDir, " with pattern ", pattern)
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

  base_names <- basename(files)
  extract_by_fixed_affixes <- function(names_vec, pattern_str) {
    wildcard_pos <- regexpr("\\.\\*", pattern_str)
    if (wildcard_pos[1] < 0) {
      wildcard_pos <- regexpr("\\*", pattern_str)
    }
    if (wildcard_pos[1] < 0) {
      return(rep(NA_character_, length(names_vec)))
    }

    prefix_pat <- substr(pattern_str, 1, wildcard_pos[1] - 1)
    suffix_pat <- substr(
      pattern_str,
      wildcard_pos[1] + attr(wildcard_pos, "match.length"),
      nchar(pattern_str)
    )

    unescape_literal <- function(x) {
      x <- sub("\\$$", "", x)
      x <- gsub("^\\^", "", x)
      x <- gsub("\\\\\\.", ".", x)
      x <- gsub("\\\\\\+", "+", x)
      x <- gsub("\\\\\\(", "(", x)
      x <- gsub("\\\\\\)", ")", x)
      x <- gsub("\\\\\\[", "[", x)
      x <- gsub("\\\\\\]", "]", x)
      x <- gsub("\\\\\\{", "{", x)
      x <- gsub("\\\\\\}", "}", x)
      x <- gsub("\\\\\\|", "|", x)
      x <- gsub("\\\\\\?", "?", x)
      x <- gsub("\\\\\\^", "^", x)
      x <- gsub("\\\\\\$", "$", x)
      x <- gsub("\\\\\\\\", "\\\\", x)
      x
    }

    prefix <- unescape_literal(prefix_pat)
    suffix <- unescape_literal(suffix_pat)

    bad_prefix <- nzchar(prefix) & !startsWith(names_vec, prefix)
    bad_suffix <- nzchar(suffix) & !endsWith(names_vec, suffix)
    out <- substring(names_vec, first = nchar(prefix) + 1)
    if (nzchar(suffix)) {
      suffix_n <- nchar(suffix)
      out <- ifelse(
        nchar(out) >= suffix_n,
        substr(out, 1, nchar(out) - suffix_n),
        ""
      )
    }
    out[bad_prefix | bad_suffix] <- NA_character_
    out
  }

  pheno <- extract_by_fixed_affixes(base_names, pattern)
  if (any(is.na(pheno) | !nzchar(pheno))) {
    # Fallback to regex capture when fixed prefix/suffix stripping is not enough.
    capture_pattern <- sub("\\.\\*", "(.*)", pattern)
    capture_pattern <- sub("\\*", "(.*)", capture_pattern, fixed = FALSE)
    matches <- regexec(capture_pattern, base_names)
    extracted <- regmatches(base_names, matches)
    pheno_regex <- vapply(extracted, function(x) {
      if (length(x) >= 2) x[2] else NA_character_
    }, character(1))

    replace_idx <- is.na(pheno) | !nzchar(pheno)
    pheno[replace_idx] <- pheno_regex[replace_idx]
  }

  if (any(is.na(pheno) | !nzchar(pheno))) {
    pheno[is.na(pheno) | !nzchar(pheno)] <- tools::file_path_sans_ext(base_names[is.na(pheno) | !nzchar(pheno)])
  }

  dup_pheno <- unique(pheno[duplicated(pheno)])
  if (length(dup_pheno) > 0) {
    dup_rows <- data.frame(pheno = pheno, file = files, stringsAsFactors = FALSE) |>
      dplyr::filter(.data$pheno %in% dup_pheno) |>
      dplyr::arrange(.data$pheno, .data$file)
    stop(
      "discover_meta_files: duplicated phenotype names extracted from filenames: ",
      paste(dup_pheno, collapse = ", "),
      "\nExamples:\n",
      paste(sprintf("  %s -> %s", dup_rows$pheno, basename(dup_rows$file)), collapse = "\n")
    )
  }

  data.frame(
    pheno = pheno,
    file = files,
    mode = "pattern",
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(.data$pheno) |>
    (\(x) {
      message("discover_meta_files: matched ", nrow(x), " file(s).")
      x
    })()
}

# Detect study columns (est/stderr) in a meta result data frame
detect_studies <- function(df) {
  est_cols <- grep("_est$", names(df), value = TRUE)
  est_cols <- setdiff(est_cols, "meta_est")
  if (length(est_cols) == 0) {
    stop("Cannot find any per-study columns ending in '_est'.")
  }
  ids <- sub("_est$", "", est_cols)
  stderr_cols <- paste0(ids, "_stderr")
  miss <- setdiff(stderr_cols, names(df))
  if (length(miss) > 0) stop("Missing stderr columns: ", paste(miss, collapse = ", "))
  list(ids = ids, est_cols = est_cols, stderr_cols = stderr_cols)
}

# Test if a row is significant at pCut
is_sig_row <- function(row, pCut = 5e-8) {
  if (is.na(pCut)) return(FALSE)
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

# Choose a headless-safe graphics device function for ggsave/base output.
#
# @param outFile Output path whose extension controls the device choice.
# @return A device function suitable for `ggsave(device=...)`.
ggsave_device_for_path <- function(outFile) {
  ext <- tolower(tools::file_ext(outFile))

  if (ext %in% c("png", "")) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      return(ragg::agg_png)
    }
    if (capabilities("cairo")) {
      return(function(filename, width, height, units = "in", res = 300, ...) {
        grDevices::png(
          filename = filename,
          width = width,
          height = height,
          units = units,
          res = res,
          type = "cairo",
          ...
        )
      })
    }
    return(grDevices::png)
  }

  if (ext %in% c("jpg", "jpeg")) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      return(ragg::agg_jpeg)
    }
    if (capabilities("cairo")) {
      return(function(filename, width, height, units = "in", res = 300, ...) {
        grDevices::jpeg(
          filename = filename,
          width = width,
          height = height,
          units = units,
          res = res,
          type = "cairo",
          ...
        )
      })
    }
    return(grDevices::jpeg)
  }

  if (ext == "pdf") {
    return(grDevices::pdf)
  }

  stop("Unsupported output extension: ", outFile)
}

# Open a graphics device based on the output file extension.
#
# @param outFile Output path.
# @param width,height,dpi Device size in inches and resolution.
open_plot_device <- function(outFile, width, height, dpi = 300) {
  ext <- tolower(tools::file_ext(outFile))
  out_dir <- dirname(outFile)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (ext %in% c("png", "")) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(outFile, width = width, height = height, units = "in", res = dpi)
      return(invisible("agg_png"))
    }
    if (capabilities("cairo")) {
      grDevices::png(outFile, width = width, height = height, units = "in", res = dpi, type = "cairo")
      return(invisible("png_cairo"))
    }
    grDevices::png(outFile, width = width, height = height, units = "in", res = dpi)
    return(invisible("png"))
  }

  if (ext %in% c("jpg", "jpeg")) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_jpeg(outFile, width = width, height = height, units = "in", res = dpi)
      return(invisible("agg_jpeg"))
    }
    if (capabilities("cairo")) {
      grDevices::jpeg(outFile, width = width, height = height, units = "in", res = dpi, type = "cairo")
      return(invisible("jpeg_cairo"))
    }
    grDevices::jpeg(outFile, width = width, height = height, units = "in", res = dpi)
    return(invisible("jpeg"))
  }

  if (ext == "pdf") {
    grDevices::pdf(outFile, width = width, height = height, onefile = FALSE)
    return(invisible("pdf"))
  }

  stop("Unsupported output extension for forest plot: ", outFile)
}

# Format numeric values for forest-plot annotations.
fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

# Format p-values for concise display.
fmt_pval <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 1e-4) return(format(x, digits = 2, scientific = TRUE))
  formatC(x, digits = 4, format = "f")
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
                           pCut = 5e-8,
                           width = NA_real_, height = NA_real_, dpi = 300,
                           plotMinP = NA_real_) {
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

  man <- man |>
    dplyr::mutate(logp = -log10(.data$P))

  if (!is.na(plotMinP)) {
    if (!is.numeric(plotMinP) || length(plotMinP) != 1 || plotMinP <= 0 || plotMinP >= 1) {
      stop("plotMinP must be a single number in (0, 1) or NA.")
    }
    plot_min_p_logp <- -log10(plotMinP)
    cap_bump <- max(0.2, plot_min_p_logp * 0.02)
    man <- man |>
      dplyr::mutate(PLOT_P = dplyr::if_else(.data$P < plotMinP, plot_min_p_logp + cap_bump, .data$logp))
  } else {
    cap_bump <- 0
    man <- man |>
      dplyr::mutate(PLOT_P = .data$logp)
  }

  auto_width <- if (is.na(width)) {
    max(11, length(unique(man$CHR)) * 0.45)
  } else {
    width
  }
  auto_height <- if (is.na(height)) 5 else height

  msg("Saving: %s", outFile)
  open_plot_device(outFile,
    width = auto_width, height = auto_height,
    dpi = dpi
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  qqman::manhattan(
    man,
    chr = "CHR",
    bp = "BP",
    p = "PLOT_P",
    snp = "SNP",
    logp = FALSE,
    col = c("grey20", "grey70"),
    cex = 0.55,
    main = title %||% "",
    ylab = expression(-log[10](italic(P))),
    ylim = c(0, max(man$PLOT_P, na.rm = TRUE) + max(0.15, cap_bump * 0.5)),
    suggestiveline = FALSE,
    genomewideline = FALSE
  )

  graphics::abline(h = -log10(5e-8), lty = 2, lwd = 1, col = "red")

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

  open_plot_device(outFile,
    width = auto_width, height = auto_height,
    dpi = dpi
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  qqman::qq(df[[pcol]], main = title %||% "QQ Plot")

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
  spacing <- 0.28
  offsets <- setNames(
    seq(from = -((length(study_lvls_all) - 1) / 2) * spacing,
        by = spacing,
        length.out = length(study_lvls_all)),
    study_lvls_all
  )
  offset_min <- min(offsets)
  offset_max <- max(offsets)
  df_plot <- df_plot |> dplyr::mutate(y_num = as.numeric(y) + offsets[as.character(Study)])

  g <- ggplot2::ggplot()

  # heterogeneity highlight background (optional)
  if (show_het && !is.null(het_y) && length(het_y) > 0) {
    df_area <- tibble::tibble(y = factor(het_y, levels = y_levels))
    df_area <- df_area |> dplyr::mutate(
      y_center = as.numeric(y),
      raw_ymin = y_center + offset_min - spacing / 2,
      raw_ymax = y_center + offset_max + spacing / 2
    )
    # clamp band to avoid overlap between adjacent rows
    band_half <- 0.48
    df_area <- df_area |> dplyr::mutate(
      ymin = pmax(y_center - band_half, raw_ymin),
      ymax = pmin(y_center + band_half, raw_ymax)
    )
    g <- g + ggplot2::geom_rect(
      data = df_area,
      ggplot2::aes(
        ymin = ymin,
        ymax = ymax,
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

  # auto width/height when not provided: keep plots readable without exceeding
  # ggsave's default 50-inch guard when the phenotype count is large.
  auto_width <- if (is.na(width)) {
    max(7, length(study_lvls_all) * 1.1 + 3)
  } else {
    width
  }
  auto_height <- if (is.na(height)) {
    per_row <- dplyr::case_when(
      length(y_levels) <= 40 ~ 0.60,
      length(y_levels) <= 80 ~ 0.45,
      TRUE ~ 0.32
    )
    min(48, max(6, length(y_levels) * per_row + 2 + if (meta_ok) 0.6 else 0))
  } else {
    height
  }

  msg("Saving: %s", outFile)
  ggplot2::ggsave(
    outFile,
    g,
    width = auto_width,
    height = auto_height,
    dpi = dpi,
    limitsize = FALSE,
    device = ggsave_device_for_path(outFile)
  )
  invisible(g)
}

# Draw a single phenotype/SNP forest plot with study information columns.
#
# Uses `metafor::forest.default()` so the figure can carry study labels,
# estimates, standard errors, inverse-variance weights, and textual CI output.
forest_plot_single_pheno <- function(r, pheno, snp, outFile,
                                     ciMult = 1.96,
                                     studyLabels = NULL,
                                     xlim_num = NULL,
                                     width = NA_real_, height = NA_real_, dpi = 300,
                                     show_meta = TRUE,
                                     show_het = TRUE) {
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required for forest plots with study details.")
  }

  studies <- detect_studies(r)
  ids <- studies$ids
  yi <- as.numeric(r[1, paste0(ids, "_est"), drop = TRUE])
  sei <- as.numeric(r[1, paste0(ids, "_stderr"), drop = TRUE])
  keep <- !is.na(yi) & !is.na(sei) & is.finite(yi) & is.finite(sei) & sei > 0

  if (!any(keep)) {
    stop("No valid study estimate/stderr pairs found for SNP ", snp, " in phenotype ", pheno, ".")
  }

  ids <- ids[keep]
  yi <- yi[keep]
  sei <- sei[keep]
  ci_lb <- yi - ciMult * sei
  ci_ub <- yi + ciMult * sei
  weights <- (1 / (sei^2))
  weights_pct <- weights / sum(weights) * 100

  if (!is.null(studyLabels) && length(studyLabels) == length(studies$ids)) {
    study_names <- studyLabels[keep]
  } else {
    study_names <- ids
  }
  study_names[study_names == "single_study"] <- "single study"

  meta_ok <- show_meta &&
    length(ids) > 1 &&
    all(c("meta_est", "meta_stderr") %in% names(r)) &&
    !is.na(r$meta_est[1]) &&
    !is.na(r$meta_stderr[1]) &&
    is.finite(r$meta_est[1]) &&
    is.finite(r$meta_stderr[1]) &&
    r$meta_stderr[1] > 0

  meta_est <- if (meta_ok) as.numeric(r$meta_est[1]) else NA_real_
  meta_se <- if (meta_ok) as.numeric(r$meta_stderr[1]) else NA_real_
  meta_lb <- if (meta_ok) meta_est - ciMult * meta_se else NA_real_
  meta_ub <- if (meta_ok) meta_est + ciMult * meta_se else NA_real_
  meta_p <- if ("meta_pval" %in% names(r)) as.numeric(r$meta_pval[1]) else NA_real_
  het_p <- if ("meta_pval.het" %in% names(r)) as.numeric(r$meta_pval.het[1]) else NA_real_

  axis_rng <- range(c(ci_lb, ci_ub, if (meta_ok) c(meta_lb, meta_ub)), na.rm = TRUE)
  if (!all(is.finite(axis_rng))) axis_rng <- c(-1, 1)
  if (diff(axis_rng) == 0) axis_rng <- axis_rng + c(-0.5, 0.5)

  if (is.null(xlim_num)) {
    pad <- max(0.15 * diff(axis_rng), 0.2)
    alim <- c(axis_rng[1] - pad, axis_rng[2] + pad)
  } else {
    alim <- xlim_num
  }

  span <- diff(alim)
  plot_xmin <- alim[1] - 2.45 * span
  plot_xmax <- alim[2] + 1.65 * span
  est_x <- alim[1] - 1.55 * span
  se_x <- alim[1] - 0.98 * span
  wt_x <- alim[1] - 0.42 * span
  ci_x <- alim[2] + 0.38 * span

  k <- length(yi)
  rows <- seq(from = k + 1, to = 2, by = -1)
  header_y <- k + 2.4
  meta_row <- 1

  auto_width <- if (is.na(width)) 10.2 else width
  auto_height <- if (is.na(height)) max(4.8, 2.7 + 0.48 * k + if (meta_ok) 0.4 else 0) else height

  tbl_est <- fmt_num(yi)
  tbl_se <- fmt_num(sei)
  tbl_wt <- paste0(fmt_num(weights_pct, digits = 1), "%")
  tbl_ci <- paste0(fmt_num(yi), " [", fmt_num(ci_lb), ", ", fmt_num(ci_ub), "]")

  summary_note <- paste0(
    "SNP: ", snp,
    " | CHR: ", as.character(r$CHR[1]),
    " | POS: ", as.character(r$POS[1]),
    if (!is.na(meta_p)) paste0(" | meta P: ", fmt_pval(meta_p)) else "",
    if (!is.na(het_p)) paste0(" | het P: ", fmt_pval(het_p)) else ""
  )

  msg("Saving: %s", outFile)
  open_plot_device(outFile, width = auto_width, height = auto_height, dpi = dpi)
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(mar = c(4.8, 4.2, 4.4, 2), xpd = NA)

  metafor::forest.default(
    x = yi,
    sei = sei,
    slab = study_names,
    rows = rows,
    xlim = c(plot_xmin, plot_xmax),
    alim = alim,
    refline = NA,
    xlab = "Effect",
    annotate = FALSE,
    header = FALSE,
    pch = 15,
    psize = 0.9 + 1.6 * sqrt(weights_pct / max(weights_pct)),
    efac = c(0, 1),
    cex = 0.9
  )

  if (meta_ok) {
    meta_col <- if (show_het && !is.na(het_p)) "#8B0000" else "#000000"
    metafor::addpoly.default(
      x = meta_est,
      sei = meta_se,
      row = meta_row,
      mlab = "Meta",
      annotate = FALSE,
      col = meta_col,
      border = meta_col
    )
    graphics::text(ci_x, meta_row, paste0(fmt_num(meta_est), " [", fmt_num(meta_lb), ", ", fmt_num(meta_ub), "]"), pos = 4, cex = 0.88, xpd = NA)
  }

  graphics::text(plot_xmin, header_y, "Study", pos = 4, font = 2, cex = 0.92, xpd = NA)
  graphics::text(est_x, header_y, "Effect", font = 2, cex = 0.92, xpd = NA)
  graphics::text(se_x, header_y, "SE", font = 2, cex = 0.92, xpd = NA)
  graphics::text(wt_x, header_y, "Weight", font = 2, cex = 0.92, xpd = NA)
  graphics::text(ci_x, header_y, "Effect [95% CI]", pos = 4, font = 2, cex = 0.92, xpd = NA)

  graphics::text(est_x, rows, tbl_est, cex = 0.88, xpd = NA)
  graphics::text(se_x, rows, tbl_se, cex = 0.88, xpd = NA)
  graphics::text(wt_x, rows, tbl_wt, cex = 0.88, xpd = NA)
  graphics::text(ci_x, rows, tbl_ci, pos = 4, cex = 0.88, xpd = NA)

  graphics::title(main = paste0("Forest: ", pheno, " @ ", snp))
  graphics::mtext(summary_note, side = 3, line = 0.6, cex = 0.78)
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
#' @param pCut Optional cutoff. When using the combined mode, SNP/phenotype
#'   pairs with best p below this value are printed to the console and written
#'   to `<outFile>_pCut.txt`.
#' @param plotMinP Optional Manhattan plotting threshold for p-value
#'   compression. Points with `P < plotMinP` are drawn slightly above
#'   `-log10(plotMinP)` instead of stretching the full y-axis.
#'
#' @return Invisibly returns the data frame passed to `qqman::manhattan()`.
#' @export
mode_big_combined <- function(metaIndex, outFile,
                              sep = "\t",
                              width = NA_real_, height = NA_real_, dpi = 300,
                              pCut = 5e-8,
                              plotMinP = NA_real_) {
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
      dplyr::mutate(feature = ph) |>
      dplyr::select(
        dplyr::any_of(c("SNP", "CHR", "POS", "meta_est", "meta_stderr", pcol, "est", "stderr")),
        .data$feature
      )

    if (nrow(keep) == 0) next

    all_hits[[length(all_hits) + 1]] <- keep
  }

  if (length(all_hits) == 0) stop("No valid hits found across any pheno (check pval column).")

  big <- dplyr::bind_rows(all_hits)

  # pick the phenotype with the smallest p for each SNP
  best <- big |>
    dplyr::group_by(.data$SNP) |>
    dplyr::arrange(.data$meta_pval, .data$feature) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()

  # optionally report pairs that beat the user threshold
  if (!is.na(pCut)) {
    hits_to_print <- best |>
      dplyr::filter(!is.na(.data$meta_pval), .data$meta_pval < pCut) |>
      dplyr::arrange(.data$meta_pval)

    if ("meta_est" %in% names(hits_to_print)) hits_to_print$est <- hits_to_print$meta_est
    if ("meta_stderr" %in% names(hits_to_print)) hits_to_print$stderr <- hits_to_print$meta_stderr
    if ("meta_pval" %in% names(hits_to_print)) hits_to_print$pval <- hits_to_print$meta_pval
    keep_cols <- c("SNP", "est", "stderr", "pval", "feature")
    keep_cols <- intersect(keep_cols, names(hits_to_print))
    hits_to_print <- hits_to_print |> dplyr::select(dplyr::all_of(keep_cols))

    if (nrow(hits_to_print) > 0) {
      msg("SNP/pheno pairs with p < %g (best per SNP):", pCut)
      apply(hits_to_print, 1, function(r) {
        p_out <- if ("meta_pval" %in% names(r)) as.numeric(r[["meta_pval"]]) else as.numeric(r[["pval"]])
        msg("  %s\t%s\t%.3e", r[["SNP"]], r[["feature"]], p_out)
        NULL
      })

      # save alongside plot, force .txt (strip any existing extension)
      base_no_ext <- tools::file_path_sans_ext(outFile)
      list_out <- paste0(base_no_ext, "_pCut.txt")
      msg("Saving list: %s", list_out)
      if (requireNamespace("data.table", quietly = TRUE)) {
        data.table::fwrite(hits_to_print, list_out, sep = "\t")
      } else {
        utils::write.table(hits_to_print, list_out, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    } else {
      msg("No SNP/pheno pairs found with p < %g (best per SNP).", pCut)
    }
  }

  plot_manhattan(
    df = best,
    outFile = outFile,
    title = "Best phenotype per SNP (minimum p across phenotypes)",
    width = width,
    height = height,
    dpi = dpi,
    plotMinP = plotMinP
  )
}


# Mode B: Manhattan for a given pheno
#' Manhattan plot for a single phenotype
#'
#' @inheritParams mode_big_combined
#' @param phenoName Phenotype name matching a row in `metaIndex$pheno`.
#' @param qqOutFile Optional output path for QQ plot.
#' @param topOutFile Optional output path for top-N hits (ordered by p-value).
#' @param top_n How many hits to keep if `topOutFile` is provided.
#'
#' @export
mode_pheno_manhattan <- function(metaIndex, phenoName, outFile,
                                 sep = "\t",
                                 width = NA_real_, height = NA_real_, dpi = 300,
                                 qqOutFile = NULL, qq_width = NA_real_, qq_height = NA_real_,
                                 topOutFile = NULL, top_n = 10,
                                 plotMinP = NA_real_) {
  row <- metaIndex |> dplyr::filter(.data$pheno == .env$phenoName)
  if (nrow(row) == 0) stop("Cannot find meta file for pheno: ", phenoName)

  df <- read_result_file(row$file[1], sep = sep, study_id = phenoName)
  plot_manhattan(
    df = df,
    outFile = outFile,
    title = paste0("Manhattan: ", phenoName),
    pCut = NA_real_,
    width = width, height = height, dpi = dpi,
    plotMinP = plotMinP
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
#' Forest plots for a SNP across all phenotypes, one file per phenotype
#'
#' @inheritParams mode_big_combined
#' @param pCut P-value cutoff used to filter phenotypes for this SNP. Use `NA`
#'   to disable filtering and keep all phenotypes containing the SNP.
#' @param snp SNP ID to plot.
#' @param ciMult Multiplier for confidence interval width.
#' @param studyLabels Optional labels replacing study IDs in the legend.
#' @param xlim_str Optional comma-separated numeric limits for the x-axis.
#' @param show_meta If TRUE, overlay meta-effect in black.
#'
#' @export
mode_snp_forest_across_phenos <- function(metaIndex, snp, outFile,
                                          pCut = 5e-8,
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
    if (!is.na(pCut)) {
      if (!is_sig_row(r[1, , drop = FALSE], pCut = pCut)) next
    }

    r$pheno <- ph
    rows[[length(rows) + 1]] <- r
  }

  if (length(rows) == 0) {
    if (!is.na(pCut)) {
      stop("No phenotypes found (with this SNP) passing significance cutoff. Try pCut=NA or loosen pCut.")
    } else {
      stop("No phenotypes contain this SNP in meta files: ", snp)
    }
  }

  df_all <- dplyr::bind_rows(rows)
  msg("Mode C: retained %d phenotype(s) for SNP %s.", nrow(df_all), snp)

  # sort phenos by significance (use pval)
  if (!("meta_pval" %in% names(df_all))) {
    stop("Missing meta_pval column after binding phenotypes.")
  }
  df_all <- df_all |>
    dplyr::mutate(score = .data$meta_pval) |>
    dplyr::arrange(.data$score)

  xlim_num <- parse_xlim(xlim_str)
  out_files <- character(nrow(df_all))

  for (i in seq_len(nrow(df_all))) {
    ph <- df_all$pheno[i]
    msg("Mode C: plotting %d/%d phenotype(s): %s", i, nrow(df_all), ph)
    ph_out <- with_suffix(outFile, paste0("_", sanitize_filename(ph)))
    forest_plot_single_pheno(
      r = df_all[i, , drop = FALSE],
      pheno = ph,
      snp = snp,
      outFile = ph_out,
      ciMult = ciMult,
      studyLabels = studyLabels,
      xlim_num = xlim_num,
      width = width, height = height, dpi = dpi,
      show_meta = show_meta,
      show_het = show_het
    )
    out_files[i] <- ph_out
  }

  invisible(data.frame(
    pheno = df_all$pheno,
    outFile = out_files,
    stringsAsFactors = FALSE
  ))
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
  row <- metaIndex |> dplyr::filter(.data$pheno == .env$pheno)
  if (nrow(row) == 0) stop("Cannot find meta file for pheno: ", pheno)

  df <- read_result_file(row$file[1], sep = sep, study_id = pheno)
  r <- df[df$SNP == snp, , drop = FALSE]
  if (nrow(r) == 0) stop("SNP not found in this pheno meta file: ", snp)

  xlim_num <- parse_xlim(xlim_str)
  forest_plot_single_pheno(
    r = r,
    pheno = pheno,
    snp = snp,
    outFile = outFile,
    ciMult = ciMult,
    studyLabels = studyLabels,
    xlim_num = xlim_num,
    width = width, height = height, dpi = dpi,
    show_meta = show_meta,
    show_het = show_het
  )
}
