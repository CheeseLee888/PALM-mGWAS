#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

option_list <- list(
  make_option("--studyDirFile", type="character", default="",
              help="txt: each line 'studyID<TAB>dir'"),
  make_option("--pattern", type="character", default="",
              help="Regex pattern for input step2 filenames [default %default]"),
  make_option("--features", type="character", default="",
              help="Optional comma-separated feature name(s) to meta-analyze [default all matched features]"),
  make_option("--metaDir", type="character", default="",
              help="Output directory for meta files"),
  make_option("--metaPrefix", type="character", default="",
              help="Output file prefix for meta files")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!nzchar(opt$studyDirFile) || !file.exists(opt$studyDirFile))
  stop("Missing/invalid --studyDirFile")
if (!nzchar(opt$metaDir))
  stop("Missing --metaDir")

# read studyDirFile (studyID \t dir)
sd <- read.table(opt$studyDirFile, header = FALSE, sep = "", stringsAsFactors = FALSE)
if (ncol(sd) < 2) stop("studyDirFile must have >=2 columns: studyID and dir")
study_dirs <- setNames(as.character(sd[[2]]), as.character(sd[[1]]))

dir.create(opt$metaDir, recursive=TRUE, showWarnings=FALSE)

feature_subset <- NULL
if (nzchar(opt$features)) {
  feature_subset <- strsplit(opt$features, ",", fixed = TRUE)[[1]]
  feature_subset <- trimws(feature_subset)
  feature_subset <- feature_subset[nzchar(feature_subset)]
  if (length(feature_subset) == 0L) {
    feature_subset <- NULL
  }
}

metaSummary(
  study_dirs = study_dirs,
  pattern    = opt$pattern,
  features   = feature_subset,
  out_dir    = opt$metaDir,
  out_prefix = opt$metaPrefix,
  keep_het   = TRUE
)
