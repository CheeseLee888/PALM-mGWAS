suppressPackageStartupMessages({
    library(PALMGWAS)
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
  make_option(c("--inputPrefix"), type = "character", default = "",
              help = "Directory containing Step3 meta files or single-study Step2 files [default %default]"),
  make_option(c("--outputPrefix"), type = "character", default = "",
              help = "Prefix for plot outputs, including directory and optional filename prefix [default %default]"),
  make_option(c("--feature"), type = "character", default = NA,
              help = "Feature name (suffix in step3/step2 filename)"),
  make_option(c("--snp"), type = "character", default = NA,
              help = "SNP ID, e.g. chr1:123:A:G (must match SNP column exactly)"),
  make_option(c("--pCut"), type = "character", default = "5e-8",
              help = "P-value cutoff for Manhattan reference lines, Manhattan/QQ point highlighting, combined hit output, and SNP-only forest filtering [default %default]"),
  make_option(c("--width"), type = "character", default = NA_character_,
              help = "Plot width inches; NA lets the script auto-size"),
  make_option(c("--height"), type = "character", default = NA_character_,
              help = "Plot height inches; NA lets the script auto-size"),
  make_option(c("--plotMinP"), type = "character", default = "NA",
              help = "Optional Manhattan and QQ plotting threshold for p-value compression; must be no larger than --pCut. Use NA to disable compression")
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
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    stop("--pCut must be a number in (0, 1).")
  }
  up <- toupper(trimws(as.character(x)))
  if (up %in% c("", "NA", "NULL")) {
    stop("--pCut must be a number in (0, 1).")
  }
  val <- suppressWarnings(as.numeric(x))
  if (is.na(val) || val <= 0 || val >= 1) {
    stop("--pCut must be a number in (0, 1).")
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
if (!is.na(plot_min_p) && plot_min_p > p_cut) {
  stop("--plotMinP should not be larger than --pCut; otherwise, the Manhattan plot cannot show the --pCut reference line.")
}

inputPrefix <- opt$inputPrefix
outputPrefix <- opt$outputPrefix
outputPrefix <- sub("_+$", "", outputPrefix)

prefixed_out <- function(suffix) {
  if (!nzchar(outputPrefix)) return(suffix)
  if (grepl("[/\\\\]$", outputPrefix)) return(paste0(outputPrefix, suffix))
  paste0(outputPrefix, "_", suffix)
}

metaIndex <- discover_meta_files(inputPrefix)

feature <- if (!is.na(opt$feature)) opt$feature else NULL
snp   <- if (!is.na(opt$snp)) opt$snp else NULL

user_supplied_pcut <- arg_supplied("pCut")
if (user_supplied_pcut && !is.null(feature) && !is.null(snp)) {
  stop("--pCut is ignored when both --feature and --snp are specified.")
}

# ----------------------------- dispatch -----------------------------

msg("InputPrefix: %s", inputPrefix)
msg("OutputPrefix: %s", outputPrefix)
msg("Found %d feature(s) across %d result file(s).", nrow(metaIndex), sum(lengths(metaIndex$files)))
msg("Resolved pCut: %s", format(p_cut, scientific = TRUE))
msg("Resolved plotMinP: %s", if (is.na(plot_min_p)) "NA" else format(plot_min_p, scientific = TRUE))
msg("Resolved width x height: %s x %s", if (is.na(width_in)) "auto" else as.character(width_in), if (is.na(height_in)) "auto" else as.character(height_in))

if (is.null(feature) && is.null(snp)) {
  # Mode A
  # Big combined plot: show best phenotype per SNP, no significance filtering.
  outFile <- prefixed_out("manhattan_combined.png")
  dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
  msg("Output file/base: %s", outFile)
  msg("Reporting mode: combined Manhattan across phenotypes.")
  msg("Mode A behavior: pCut reference/highlighting enabled at %s", format(p_cut, scientific = TRUE))
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
  msg("Mode B behavior: pCut reference/highlighting enabled at %s", format(p_cut, scientific = TRUE))
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
    pCut = p_cut,
    plotMinP = plot_min_p
  )

} else if (is.null(feature) && !is.null(snp)) {
  # Mode C
  outFile <- prefixed_out(paste0("forest_", sanitize_filename(snp), ".png"))
  dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
  msg("Output file/base: %s", outFile)
  msg("Reporting mode: per-phenotype forest plots for SNP %s.", snp)
  msg("Mode C behavior: pCut %s; one file is generated for each retained phenotype.",
      paste0("enabled at ", format(p_cut, scientific = TRUE)))
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
    show_meta = TRUE,
    show_het = TRUE
  )

} else {
  # Mode D
  outFile <- prefixed_out(paste0("forest_", sanitize_filename(snp), "_", sanitize_filename(feature), ".png"))
  dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
  msg("Output file/base: %s", outFile)
  msg("Reporting mode: forest for feature %s and SNP %s.", feature, snp)
  if (user_supplied_pcut) {
    msg("Mode D behavior: pCut is ignored in this mode.")
  }
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
    show_meta = TRUE,
    show_het = TRUE
  )
}

msg("Done.")
