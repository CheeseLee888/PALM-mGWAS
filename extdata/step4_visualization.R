suppressPackageStartupMessages({
    library(PALMmGWAS)
    library(stringr)
    library(optparse)
})

option_list <- list(
  make_option(c("--metaDir"), type = "character", default = "output/meta",
              help = "Directory containing step3 meta files [default %default]"),
  make_option(c("--plotDir"), type = "character", default = "output/plot",
              help = "Directory to save plots [default %default]"),
  make_option(c("--methodPrefix"), type = "character", default = NA,
              help = "Optional method prefix filter, e.g. palm1 or palm2"),
  make_option(c("--pheno"), type = "character", default = NA,
              help = "Phenotype name (suffix in step3 meta filename)"),
  make_option(c("--snp"), type = "character", default = NA,
              help = "SNP ID, e.g. chr1:123:A:G (must match SNP column exactly)"),
  make_option(c("--out"), type = "character", default = NA,
              help = "Output file name (png/pdf). If not set, auto-named."),
  make_option(c("--sep"), type = "character", default = "\t",
              help = "Input delimiter [default tab]"),

  make_option(c("--prefer"), type = "character", default = "qval",
              help = "Prefer 'qval' or 'pval' for significance & y-axis in manhattan [default %default]"),
  make_option(c("--qCut"), type = "double", default = 0.05,
              help = "q-value cutoff [default %default]"),
  make_option(c("--pCut"), type = "double", default = 1e-2,
              help = "p-value cutoff (used if qval missing) [default %default]"),

  make_option(c("--sigOnlyPheno"), action = "store_true", default = FALSE,
              help = "When --snp only: only plot phenos passing significance cutoff (default: plot ALL phenos containing SNP)."),
  make_option(c("--onlySigBig"), action = "store_true", default = TRUE,
              help = "When no pheno/snp: big plot draws only significant hits (default TRUE)."),
  make_option(c("--maxPoints"), type = "integer", default = 200000,
              help = "Max points to plot (downsample by significance) [default %default]"),

  make_option(c("--ciMult"), type = "double", default = 1.96,
              help = "CI multiplier (1.96 ~ 95%% CI) [default %default]"),
  make_option(c("--hetQvalCut"), type = "double", default = 0.1,
              help = "Highlight heterogeneity if qval.het <= this [default %default]"),

  make_option(c("--studyLabels"), type = "character", default = NA,
              help = "Optional comma-separated labels for studies, e.g. 'FR-CRC,DE-CRC,UKB' (must match number of study*_est columns)."),
  make_option(c("--showMeta"), action = "store_true", default = TRUE,
              help = "Overlay meta est/stderr in black if available [default TRUE]"),

  make_option(c("--xlim"), type = "character", default = NA,
              help = "Effect axis limits for forest plots: 'min,max' (e.g. '-1.5,1.5')"),
  make_option(c("--width"), type = "double", default = 10,
              help = "Plot width inches [default %default]"),
  make_option(c("--height"), type = "double", default = 6,
              help = "Plot height inches [default %default]"),
  make_option(c("--dpi"), type = "integer", default = 300,
              help = "DPI for png [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

prefer <- match.arg(tolower(opt$prefer), choices = c("qval", "pval"))

metaDir <- opt$metaDir
plotDir <- opt$plotDir
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

methodPrefix <- if (!is.na(opt$methodPrefix)) opt$methodPrefix else NULL
metaIndex <- discover_meta_files(metaDir, methodPrefix = methodPrefix)

pheno <- if (!is.na(opt$pheno)) opt$pheno else NULL
snp   <- if (!is.na(opt$snp)) opt$snp else NULL

studyLabels <- NULL
if (!is.na(opt$studyLabels)) {
  studyLabels <- str_split(opt$studyLabels, ",", simplify = TRUE) |> as.character()
  studyLabels <- trimws(studyLabels)
  if (length(studyLabels) == 0) studyLabels <- NULL
}

# auto output filename
auto_out <- function() {
  base <- if (!is.null(methodPrefix)) methodPrefix else "palm"
  if (is.null(pheno) && is.null(snp)) {
    return(file.path(plotDir, paste0(base, "_step4_combined_hits.png")))
  }
  if (!is.null(pheno) && is.null(snp)) {
    return(file.path(plotDir, paste0(base, "_step4_manhattan_", sanitize_filename(pheno), ".png")))
  }
  if (is.null(pheno) && !is.null(snp)) {
    tag <- sanitize_filename(snp)
    if (opt$sigOnlyPheno) {
      return(file.path(plotDir, paste0(base, "_step4_forest_snp_", tag, "_sigPheno.png")))
    } else {
      return(file.path(plotDir, paste0(base, "_step4_forest_snp_", tag, "_allPheno.png")))
    }
  }
  # both
  return(file.path(plotDir, paste0(base, "_step4_forest_", sanitize_filename(pheno), "_", sanitize_filename(snp), ".png")))
}

outFile <- if (!is.na(opt$out)) opt$out else auto_out()

# ----------------------------- dispatch -----------------------------

msg("MetaDir: %s", metaDir)
msg("PlotDir: %s", plotDir)
msg("Found %d meta files.", nrow(metaIndex))

if (is.null(pheno) && is.null(snp)) {
  # Mode A
  # Big combined plot: by definition only significant hits; controlled by --onlySigBig (default TRUE)
  # If user sets onlySigBig=FALSE, it could explode; we keep the flag but still cap by maxPoints.
  if (opt$onlySigBig) {
    mode_big_combined(
      metaIndex = metaIndex,
      outFile = outFile,
      qCut = opt$qCut,
      pCut = opt$pCut,
      prefer = prefer,
      maxPoints = opt$maxPoints,
      sep = opt$sep,
      width = opt$width, height = opt$height, dpi = opt$dpi
    )
  } else {
    # If not onlySigBig, we still read all but then it will be huge.
    # Safer behavior: still only plot top maxPoints by significance from all files.
    msg("Warning: --onlySigBig FALSE can be huge; will still cap to --maxPoints by best significance.")
    mode_big_combined(
      metaIndex = metaIndex,
      outFile = outFile,
      qCut = 1,                 # effectively keep all for qval
      pCut = 1,                 # keep all for pval
      prefer = prefer,
      maxPoints = opt$maxPoints,
      sep = opt$sep,
      width = opt$width, height = opt$height, dpi = opt$dpi
    )
  }

} else if (!is.null(pheno) && is.null(snp)) {
  # Mode B
  mode_pheno_manhattan(
    metaIndex = metaIndex,
    phenoName = pheno,
    outFile = outFile,
    qCut = opt$qCut,
    pCut = opt$pCut,
    prefer = prefer,
    sep = opt$sep,
    onlySig = FALSE,           # Manhattan normally shows all; you can change if you want
    maxPoints = opt$maxPoints,
    width = opt$width, height = opt$height, dpi = opt$dpi
  )

} else if (is.null(pheno) && !is.null(snp)) {
  # Mode C
  mode_snp_forest_across_phenos(
    metaIndex = metaIndex,
    snp = snp,
    outFile = outFile,
    qCut = opt$qCut,
    pCut = opt$pCut,
    prefer = prefer,
    sigOnlyPheno = opt$sigOnlyPheno,     # user requested optional; default FALSE => draw all phenos
    ciMult = opt$ciMult,
    hetQvalCut = opt$hetQvalCut,
    studyLabels = studyLabels,
    sep = opt$sep,
    xlim_str = opt$xlim,
    width = opt$width, height = opt$height, dpi = opt$dpi,
    show_meta = opt$showMeta
  )

} else {
  # Mode D
  mode_pheno_snp_forest(
    metaIndex = metaIndex,
    pheno = pheno,
    snp = snp,
    outFile = outFile,
    ciMult = opt$ciMult,
    hetQvalCut = opt$hetQvalCut,
    studyLabels = studyLabels,
    sep = opt$sep,
    xlim_str = opt$xlim,
    width = opt$width, height = opt$height, dpi = opt$dpi,
    show_meta = opt$showMeta
  )
}

msg("Done.")
