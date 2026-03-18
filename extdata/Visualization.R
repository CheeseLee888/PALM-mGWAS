suppressPackageStartupMessages({
    library(PALMmGWAS)
    library(optparse)
})

raw_args <- commandArgs(trailingOnly = TRUE)

arg_supplied <- function(flag, args = raw_args) {
  bare <- paste0("--", flag)
  prefix <- paste0(bare, "=")
  any(args == bare | startsWith(args, prefix))
}

option_list <- list(
  make_option(c("--metaDir"), type = "character", default = "",
              help = "Directory containing step3 meta files [default %default]"),
  make_option(c("--plotDir"), type = "character", default = "",
              help = "Directory to save plots [default %default]"),
  make_option(c("--pattern"), type = "character", default = "",
              help = "Regex pattern for meta filenames [default %default]"),
  make_option(c("--pheno"), type = "character", default = NA,
              help = "Phenotype name (suffix in step3 meta filename)"),
  make_option(c("--snp"), type = "character", default = NA,
              help = "SNP ID, e.g. chr1:123:A:G (must match SNP column exactly)"),
  make_option(c("--PlotOutputFile"), type = "character", default = NA,
              help = "Output file name (png/pdf). If not set, auto-named."),
  make_option(c("--pCut"), type = "character", default = "1e-5",
              help = "When neither --pheno nor --snp is given, print SNP/pheno pairs whose best p across phenos is below this cutoff; when only --snp is given, filter phenotypes by this cutoff. Use NA to disable filtering/printing [default %default]"),
  make_option(c("--showMeta"), type = "logical", default = TRUE,
              help = "Only valid when --snp is given; overlay meta est/stderr in black if available [default TRUE]"),
  make_option(c("--showHet"), type = "logical", default = TRUE,
              help = "Only valid when --snp is given; highlight heterogeneity rows in yellow [default TRUE]"),
  make_option(c("--width"), type = "character", default = NA_character_,
              help = "Plot width inches; NA lets the script auto-size"),
  make_option(c("--height"), type = "character", default = NA_character_,
              help = "Plot height inches; NA lets the script auto-size"),
  make_option(c("--manhattanCap"), type = "character", default = NA_character_,
              help = "Optional Manhattan y-axis cap on the -log10(P) scale; a dashed line is drawn at the cap and points above it are shown slightly above that line")
)

opt <- parse_args(OptionParser(option_list = option_list))

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

parse_positive_numeric <- function(x, flag_name) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  if (is.na(x)) return(NA_real_)
  up <- toupper(trimws(as.character(x)))
  if (up %in% c("", "NA", "NULL")) return(NA_real_)
  val <- suppressWarnings(as.numeric(x))
  if (is.na(val) || val <= 0) {
    stop(flag_name, " must be a positive number or NA.")
  }
  val
}

width_in  <- parse_dim(opt$width)
height_in <- parse_dim(opt$height)
p_cut <- parse_pcut(opt$pCut)
manhattan_cap <- parse_positive_numeric(opt$manhattanCap, "--manhattanCap")

metaDir <- opt$metaDir
plotDir <- opt$plotDir
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

metaIndex <- discover_meta_files(metaDir, pattern = opt$pattern)

pheno <- if (!is.na(opt$pheno)) opt$pheno else NULL
snp   <- if (!is.na(opt$snp)) opt$snp else NULL

if (arg_supplied("pCut") && !((is.null(pheno) && is.null(snp)) || (is.null(pheno) && !is.null(snp)))) {
  stop("--pCut can only be used when neither --pheno nor --snp is specified, or when only --snp is specified.")
}

if (arg_supplied("showMeta") && is.null(snp)) {
  stop("--showMeta requires --snp.")
}

if (arg_supplied("showHet") && is.null(snp)) {
  stop("--showHet requires --snp.")
}

# auto output filename
auto_out <- function() {
  base <- ""
  if (is.null(pheno) && is.null(snp)) {
    return(file.path(plotDir, paste0(base, "combined_hits.png")))
  }
  if (!is.null(pheno) && is.null(snp)) {
    return(file.path(plotDir, paste0(base, "manhattan_", sanitize_filename(pheno), ".png")))
  }
  if (is.null(pheno) && !is.null(snp)) {
    tag <- sanitize_filename(snp)
    return(file.path(plotDir, paste0(base, "forest_", tag, ".png")))
  }
  # both
  return(file.path(plotDir, paste0(base, "forest_", sanitize_filename(pheno), "_", sanitize_filename(snp), ".png")))
}

outFile <- if (!is.na(opt$PlotOutputFile)) opt$PlotOutputFile else auto_out()

# ----------------------------- dispatch -----------------------------

msg("MetaDir: %s", metaDir)
msg("PlotDir: %s", plotDir)
msg("Pattern: %s", opt$pattern)
msg("Found %d files.", nrow(metaIndex))
msg("Resolved pCut: %s", if (is.na(p_cut)) "NA" else format(p_cut, scientific = TRUE))
msg("Resolved manhattanCap: %s", if (is.na(manhattan_cap)) "NA" else as.character(manhattan_cap))
msg("Resolved width x height: %s x %s", if (is.na(width_in)) "auto" else as.character(width_in), if (is.na(height_in)) "auto" else as.character(height_in))
msg("Output file: %s", outFile)

if (is.null(pheno) && is.null(snp)) {
  # Mode A
  # Big combined plot: show best phenotype per SNP, no significance filtering.
  msg("Visualization mode: combined Manhattan across phenotypes.")
  msg("Mode A behavior: pCut %s", if (is.na(p_cut)) "disabled" else paste0("enabled at ", format(p_cut, scientific = TRUE)))
  mode_big_combined(
    metaIndex = metaIndex,
    outFile = outFile,
    sep = "\t",
    width = width_in, height = height_in, dpi = 300,
    pCut = p_cut,
    manhattanCap = manhattan_cap
  )

} else if (!is.null(pheno) && is.null(snp)) {
  # Mode B
  msg("Visualization mode: Manhattan and qq for phenotype %s.", pheno)
  msg("Mode B behavior: pCut is ignored in this mode.")
  # keep auxiliary outputs aligned with main outFile
  base_no_ext <- sub("\\.[^.]+$", "", outFile)
  qq_out <- paste0(base_no_ext, "_qq.png")
  top_out <- paste0(base_no_ext, "_top10.txt")
  mode_pheno_manhattan(
    metaIndex = metaIndex,
    phenoName = pheno,
    outFile = outFile,
    sep = "\t",
    width = width_in, height = height_in, dpi = 300,
    qqOutFile = qq_out,
    topOutFile = top_out,
    top_n = 10,
    manhattanCap = manhattan_cap
  )

} else if (is.null(pheno) && !is.null(snp)) {
  # Mode C
  msg("Visualization mode: forest across phenotypes for SNP %s.", snp)
  msg("Mode C behavior: pCut %s; showMeta=%s; showHet=%s",
      if (is.na(p_cut)) "disabled" else paste0("enabled at ", format(p_cut, scientific = TRUE)),
      opt$showMeta, opt$showHet)
  if (!is.na(manhattan_cap)) {
    msg("Mode C behavior: manhattanCap is ignored in this mode.")
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
  msg("Visualization mode: forest for phenotype %s and SNP %s.", pheno, snp)
  msg("Mode D behavior: showMeta=%s; showHet=%s; pCut is ignored in this mode.", opt$showMeta, opt$showHet)
  if (!is.na(manhattan_cap)) {
    msg("Mode D behavior: manhattanCap is ignored in this mode.")
  }
  mode_pheno_snp_forest(
    metaIndex = metaIndex,
    pheno = pheno,
    snp = snp,
    outFile = outFile,
    sep = "\t",
    width = width_in, height = height_in, dpi = 300,
    show_meta = opt$showMeta,
    show_het = opt$showHet
  )
}

msg("Done.")
