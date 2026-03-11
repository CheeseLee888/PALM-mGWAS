#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

option_list <- list(
  make_option("--studyFile", type="character", default="",
              help="TSV: each line 'studyID<TAB>dir'"),
  make_option("--pattern", type="character", default="",
              help="Regex pattern for input step2 filenames [default %default]"),
  make_option("--metaDir", type="character", default="",
              help="Output directory for meta files"),
  make_option("--metaPrefix", type="character", default="",
              help="Output file prefix for meta files")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!nzchar(opt$studyFile) || !file.exists(opt$studyFile))
  stop("Missing/invalid --studyFile")
if (!nzchar(opt$metaDir))
  stop("Missing --metaDir")

# read studyFile (studyID \t dir)
sd <- read.table(opt$studyFile, header = FALSE, sep = "", stringsAsFactors = FALSE)
if (ncol(sd) < 2) stop("studyFile must have >=2 columns: studyID and dir")
study_dirs <- setNames(as.character(sd[[2]]), as.character(sd[[1]]))

dir.create(opt$metaDir, recursive=TRUE, showWarnings=FALSE)

metaSummary(
  study_dirs = study_dirs,
  pattern    = opt$pattern,
  features   = NULL,
  out_dir    = opt$metaDir,
  out_prefix = opt$metaPrefix,
  keep_het   = TRUE
)
