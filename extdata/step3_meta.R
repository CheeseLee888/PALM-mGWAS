#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

option_list <- list(
  make_option("--studyDirList", type="character", default="",
              help="TSV: each line 'studyID<TAB>dir'"),
  make_option("--inPrefix", type="character", default="",
              help="Input step2 prefix"),
  make_option("--outDir", type="character", default="",
              help="Output directory"),
  make_option("--outPrefix", type="character", default="",
              help="Output file prefix"),
  make_option("--keepHet", type="character", default="TRUE",
              help="TRUE/FALSE"),
  make_option("--metaMethod", type="character", default="EE",
              help="metafor::rma method"),
  make_option("--featuresFile", type="character", default="NULL",
              help="Optional: feature names, one per line; NULL=infer")
)

opt <- parse_args(OptionParser(option_list = option_list))

# normalize "NULL"
if (toupper(opt$featuresFile) == "NULL" || !nzchar(opt$featuresFile)) opt$featuresFile <- NULL
opt$keepHet <- toupper(opt$keepHet) %in% c("TRUE","T","1","YES","Y")

if (!nzchar(opt$studyDirList) || !file.exists(opt$studyDirList))
  stop("Missing/invalid --studyDirList")
if (!nzchar(opt$outDir))
  stop("Missing --outDir")

# read studyDirList (studyID \t dir)
sd <- read.table(opt$studyDirList, header = FALSE, sep = "", stringsAsFactors = FALSE)
if (ncol(sd) < 2) stop("studyDirList must have >=2 columns: studyID and dir")
study_dirs <- setNames(as.character(sd[[2]]), as.character(sd[[1]]))


study_dirs <- setNames(as.character(sd[[2]]), as.character(sd[[1]]))

# optional features
features <- NULL
if (!is.null(opt$featuresFile)) {
  if (!file.exists(opt$featuresFile)) stop("Cannot find --featuresFile")
  features <- readLines(opt$featuresFile, warn=FALSE)
  features <- features[nzchar(features)]
  if (length(features) == 0) features <- NULL
}

dir.create(opt$outDir, recursive=TRUE, showWarnings=FALSE)

metaSummary(
  study_dirs = study_dirs,
  in_prefix  = opt$inPrefix,
  features   = features,
  out_dir    = opt$outDir,
  out_prefix = opt$outPrefix,
  keep_het   = opt$keepHet,
  meta.method = opt$metaMethod
)
