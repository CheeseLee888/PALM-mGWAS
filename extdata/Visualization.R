suppressPackageStartupMessages({
    library(PALMmGWAS)
    library(optparse)
})

option_list <- list(
  make_option(c("--metaDir"), type = "character", default = "",
              help = "Directory containing step3 meta files [default %default]"),
  make_option(c("--plotDir"), type = "character", default = "",
              help = "Directory to save plots [default %default]"),
  make_option(c("--pheno"), type = "character", default = NA,
              help = "Phenotype name (suffix in step3 meta filename)"),
  make_option(c("--snp"), type = "character", default = NA,
              help = "SNP ID, e.g. chr1:123:A:G (must match SNP column exactly)"),
  make_option(c("--out"), type = "character", default = NA,
              help = "Output file name (png/pdf). If not set, auto-named."),
  make_option(c("--pCut"), type = "double", default = 1e-5,
              help = "p-value cutoff (used if qval missing) [default %default]"),

  make_option(c("--sigOnlyPheno"), type = "logical", default = FALSE,
              help = "When --snp only: only plot phenos passing significance cutoff (default: plot ALL phenos containing SNP)."),
  make_option(c("--printCut"), type = "double", default = 1e-8,
              help = "When no pheno/snp: print SNP/pheno pairs whose best p across phenos is below this cutoff [default %default]."),
  make_option(c("--showMeta"), type = "logical", default = TRUE,
              help = "Overlay meta est/stderr in black if available [default TRUE]"),
  make_option(c("--showHet"), type = "logical", default = TRUE,
              help = "Highlight heterogeneity rows in yellow [default TRUE]"),
  make_option(c("--width"), type = "double", default = 10,
              help = "Plot width inches [default %default]"),
  make_option(c("--height"), type = "double", default = 6,
              help = "Plot height inches [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))


metaDir <- opt$metaDir
plotDir <- opt$plotDir
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

metaIndex <- discover_meta_files(metaDir)

pheno <- if (!is.na(opt$pheno)) opt$pheno else NULL
snp   <- if (!is.na(opt$snp)) opt$snp else NULL

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
    if (opt$sigOnlyPheno) {
      return(file.path(plotDir, paste0(base, "forest_snp_", tag, "_sigPheno.png")))
    } else {
      return(file.path(plotDir, paste0(base, "forest_snp_", tag, "_allPheno.png")))
    }
  }
  # both
  return(file.path(plotDir, paste0(base, "forest_", sanitize_filename(pheno), "_", sanitize_filename(snp), ".png")))
}

outFile <- if (!is.na(opt$out)) opt$out else auto_out()

# ----------------------------- dispatch -----------------------------

msg("MetaDir: %s", metaDir)
msg("PlotDir: %s", plotDir)
msg("Found %d meta files.", nrow(metaIndex))

if (is.null(pheno) && is.null(snp)) {
  # Mode A
  # Big combined plot: show best phenotype per SNP, no significance filtering.
  mode_big_combined(
    metaIndex = metaIndex,
    outFile = outFile,
    sep = "\t",
    width = opt$width, height = opt$height, dpi = 300,
    printCut = opt$printCut
  )

} else if (!is.null(pheno) && is.null(snp)) {
  # Mode B
  # keep auxiliary outputs aligned with main outFile
  base_no_ext <- sub("\\.[^.]+$", "", outFile)
  qq_out <- paste0(base_no_ext, "_qq.png")
  top_out <- paste0(base_no_ext, "_top10.txt")
  mode_pheno_manhattan(
    metaIndex = metaIndex,
    phenoName = pheno,
    outFile = outFile,
    pCut = opt$pCut,
    sep = "\t",
    onlySig = FALSE,           # Manhattan normally shows all; you can change if you want
    width = opt$width, height = opt$height, dpi = 300,
    qqOutFile = qq_out,
    topOutFile = top_out,
    top_n = 10
  )

} else if (is.null(pheno) && !is.null(snp)) {
  # Mode C
  mode_snp_forest_across_phenos(
    metaIndex = metaIndex,
    snp = snp,
    outFile = outFile,
    pCut = opt$pCut,
    sigOnlyPheno = opt$sigOnlyPheno,     # user requested optional; default FALSE => draw all phenos
    sep = "\t",
    width = opt$width, height = opt$height, dpi = 300,
    show_meta = opt$showMeta,
    show_het = opt$showHet
  )

} else {
  # Mode D
  mode_pheno_snp_forest(
    metaIndex = metaIndex,
    pheno = pheno,
    snp = snp,
    outFile = outFile,
    sep = "\t",
    width = opt$width, height = opt$height, dpi = 300,
    show_meta = opt$showMeta,
    show_het = opt$showHet
  )
}

msg("Done.")
