suppressPackageStartupMessages({
    library(PALMmGWAS)
    library(optparse)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
  local_impl <- file.path(dirname(script_path), "..", "R", "reporting.R")
  if (file.exists(local_impl)) {
    source(local_impl)
  }
}

default_r_plot_file <- "Rplots.pdf"
old_device_option <- getOption("device")
options(device = function(...) {
  grDevices::pdf(
    file = tempfile(pattern = "Rplots_", tmpdir = tempdir(), fileext = ".pdf"),
    ...
  )
})
on.exit(options(device = old_device_option), add = TRUE)

raw_args <- commandArgs(trailingOnly = TRUE)

arg_supplied <- function(flag, args = raw_args) {
  bare <- paste0("--", flag)
  prefix <- paste0(bare, "=")
  any(args == bare | startsWith(args, prefix))
}

option_list <- list(
  make_option(c("--metaDir"), type = "character", default = "",
              help = "Directory containing step3 meta files [default %default]"),
  make_option(c("--plotPrefix"), type = "character", default = "",
              help = "Prefix for plot outputs, including directory and optional filename prefix [default %default]"),
  make_option(c("--feature"), type = "character", default = NA,
              help = "Feature name (suffix in step3/step2 filename)"),
  make_option(c("--snp"), type = "character", default = NA,
              help = "SNP ID, e.g. chr1:123:A:G (must match SNP column exactly)"),
  make_option(c("--pCut"), type = "character", default = "5e-8",
              help = "When neither --feature nor --snp is given, print SNP/feature pairs whose best p across features is below this cutoff; when only --snp is given, filter features by this cutoff before writing one forest plot per feature. Use NA to disable filtering/printing [default %default]"),
  make_option(c("--showMeta"), type = "logical", default = TRUE,
              help = "Only valid when --snp is given; overlay meta est/stderr in black if available [default TRUE]"),
  make_option(c("--showHet"), type = "logical", default = TRUE,
              help = "Only valid when --snp is given; emphasize heterogeneity information in SNP-based forest plots when available [default TRUE]"),
  make_option(c("--width"), type = "character", default = NA_character_,
              help = "Plot width inches; NA lets the script auto-size"),
  make_option(c("--height"), type = "character", default = NA_character_,
              help = "Plot height inches; NA lets the script auto-size"),
  make_option(c("--plotMinP"), type = "character", default = "NA",
              help = "Optional Manhattan plotting threshold for p-value compression; points with P < plotMinP are compressed near the threshold and colored red instead of stretching the full y-axis. Use NA to disable compression")
)

opt <- parse_args(OptionParser(option_list = option_list))

cleanup_default_rplots <- function() {
  if (!file.exists(default_r_plot_file)) return(invisible(FALSE))
  unlink(default_r_plot_file)
  invisible(TRUE)
}

cleanup_default_rplots()
on.exit(cleanup_default_rplots(), add = TRUE)

# coerce width/height strings (including "NA"/"null"/empty) to numeric or NA_real_
parse_dim <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  if (is.na(x)) return(NA_real_)
  if (is.character(x)) {
    up <- toupper(trimws(x))
    if (up %in% c("", "NA", "NULL")) return(NA_real_)
  }
  val <- suppressWarnings(as.numeric(x))
  if (is.na(val)) return(NA_real_)
  val
}

parse_pcut <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  if (is.na(x)) return(NA_real_)
  up <- toupper(trimws(as.character(x)))
  if (up %in% c("", "NA", "NULL")) return(NA_real_)
  val <- suppressWarnings(as.numeric(x))
  if (is.na(val)) {
    stop("--pCut must be numeric or NA.")
  }
  val
}

parse_probability <- function(x, flag_name) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  if (is.na(x)) return(NA_real_)
  up <- toupper(trimws(as.character(x)))
  if (up %in% c("", "NA", "NULL")) return(NA_real_)
  val <- suppressWarnings(as.numeric(x))
  if (is.na(val) || val <= 0 || val >= 1) {
    stop(flag_name, " must be a number in (0, 1) or NA.")
  }
  val
}

width_in  <- parse_dim(opt$width)
height_in <- parse_dim(opt$height)
p_cut <- parse_pcut(opt$pCut)
plot_min_p <- parse_probability(opt$plotMinP, "--plotMinP")

metaDir <- opt$metaDir
plotPrefix <- opt$plotPrefix
plotPrefix <- sub("_+$", "", plotPrefix)

prefixed_out <- function(suffix) {
  if (!nzchar(plotPrefix)) return(suffix)
  if (grepl("[/\\\\]$", plotPrefix)) return(paste0(plotPrefix, suffix))
  paste0(plotPrefix, "_", suffix)
}

metaIndex <- discover_meta_files(metaDir)

feature <- if (!is.na(opt$feature)) opt$feature else NULL
snp   <- if (!is.na(opt$snp)) opt$snp else NULL

if (arg_supplied("pCut") && !((is.null(feature) && is.null(snp)) || (is.null(feature) && !is.null(snp)))) {
  stop("--pCut can only be used when neither --feature nor --snp is specified, or when only --snp is specified.")
}

if (arg_supplied("showMeta") && is.null(snp)) {
  stop("--showMeta requires --snp.")
}

if (arg_supplied("showHet") && is.null(snp)) {
  stop("--showHet requires --snp.")
}

# ----------------------------- dispatch -----------------------------

msg("MetaDir: %s", metaDir)
msg("PlotPrefix: %s", plotPrefix)
msg("Found %d files.", nrow(metaIndex))
msg("Resolved pCut: %s", if (is.na(p_cut)) "NA" else format(p_cut, scientific = TRUE))
msg("Resolved plotMinP: %s", if (is.na(plot_min_p)) "NA" else format(plot_min_p, scientific = TRUE))
msg("Resolved width x height: %s x %s", if (is.na(width_in)) "auto" else as.character(width_in), if (is.na(height_in)) "auto" else as.character(height_in))

if (is.null(feature) && is.null(snp)) {
  # Mode A
  # Big combined plot: show best phenotype per SNP, no significance filtering.
  outFile <- prefixed_out("combined_hits.png")
  dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
  msg("Output file/base: %s", outFile)
  msg("Reporting mode: combined Manhattan across phenotypes.")
  msg("Mode A behavior: pCut %s", if (is.na(p_cut)) "disabled" else paste0("enabled at ", format(p_cut, scientific = TRUE)))
  mode_big_combined(
    metaIndex = metaIndex,
    outFile = outFile,
    sep = "\t",
    width = width_in, height = height_in, dpi = 300,
    pCut = p_cut,
    plotMinP = plot_min_p
  )

} else if (!is.null(feature) && is.null(snp)) {
  # Mode B
  outFile <- prefixed_out(paste0("manhattan_", sanitize_filename(feature), ".png"))
  qq_out <- prefixed_out(paste0("qq_", sanitize_filename(feature), ".png"))
  dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(qq_out), recursive = TRUE, showWarnings = FALSE)
  msg("Output file/base: %s", outFile)
  msg("QQ output file: %s", qq_out)
  msg("Reporting mode: Manhattan and qq for feature %s.", feature)
  msg("Mode B behavior: pCut is ignored in this mode.")
  # keep auxiliary outputs aligned with main outFile
  base_no_ext <- sub("\\.[^.]+$", "", outFile)
  top_out <- paste0(base_no_ext, "_top10.txt")
  mode_pheno_manhattan(
    metaIndex = metaIndex,
    phenoName = feature,
    outFile = outFile,
    sep = "\t",
    width = width_in, height = height_in, dpi = 300,
    qqOutFile = qq_out,
    topOutFile = top_out,
    top_n = 10,
    plotMinP = plot_min_p
  )

} else if (is.null(feature) && !is.null(snp)) {
  # Mode C
  outFile <- prefixed_out(paste0("forest_", sanitize_filename(snp), ".png"))
  dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
  msg("Output file/base: %s", outFile)
  msg("Reporting mode: per-phenotype forest plots for SNP %s.", snp)
  msg("Mode C behavior: pCut %s; showMeta=%s; showHet=%s; one file is generated for each retained phenotype.",
      if (is.na(p_cut)) "disabled" else paste0("enabled at ", format(p_cut, scientific = TRUE)),
      opt$showMeta, opt$showHet)
  if (!is.na(plot_min_p)) {
    msg("Mode C behavior: plotMinP is ignored in this mode.")
  }
  mode_snp_forest_across_phenos(
    metaIndex = metaIndex,
    snp = snp,
    outFile = outFile,
    pCut = p_cut,
    sep = "\t",
    width = width_in, height = height_in, dpi = 300,
    show_meta = opt$showMeta,
    show_het = opt$showHet
  )

} else {
  # Mode D
  outFile <- prefixed_out(paste0("forest_", sanitize_filename(feature), "_", sanitize_filename(snp), ".png"))
  dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
  msg("Output file/base: %s", outFile)
  msg("Reporting mode: forest for feature %s and SNP %s.", feature, snp)
  msg("Mode D behavior: showMeta=%s; showHet=%s; pCut is ignored in this mode.", opt$showMeta, opt$showHet)
  if (!is.na(plot_min_p)) {
    msg("Mode D behavior: plotMinP is ignored in this mode.")
  }
  mode_pheno_snp_forest(
    metaIndex = metaIndex,
    pheno = feature,
    snp = snp,
    outFile = outFile,
    sep = "\t",
    width = width_in, height = height_in, dpi = 300,
    show_meta = opt$showMeta,
    show_het = opt$showHet
  )
}

msg("Done.")
