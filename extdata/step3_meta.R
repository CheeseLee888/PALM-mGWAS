#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
  local_impl <- file.path(dirname(script_path), "..", "R", "metaSummary.R")
  if (file.exists(local_impl)) {
    source(local_impl)
  }
}

option_list <- list(
  make_option("--studyDirFile", type="character", default="",
              help="txt: each line 'studyID<TAB>dir'"),
  make_option("--inputPrefix", type="character", default="",
              help="Shared Step2 base prefix [default %default]"),
  make_option("--chrom", type="character", default="NULL",
              help="Step2 scope: NULL for allchr, or 1..22 [default %default]"),
  make_option("--featureColList", type="character", default="NULL",
              help="Optional comma-separated feature name(s) to meta-analyze [default %default]"),
  make_option("--metaPrefix", type="character", default="",
              help="Full output prefix for meta files, e.g. example/output/meta/step3_meta")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!nzchar(opt$studyDirFile) || !file.exists(opt$studyDirFile))
  stop("Missing/invalid --studyDirFile")
if (!nzchar(opt$inputPrefix))
  stop("Missing --inputPrefix")
if (!nzchar(opt$metaPrefix))
  stop("Missing --metaPrefix")

# read studyDirFile (studyID \t dir)
sd <- read.table(opt$studyDirFile, header = FALSE, sep = "", stringsAsFactors = FALSE)
if (ncol(sd) < 2) stop("studyDirFile must have >=2 columns: studyID and dir")
study_dirs <- setNames(as.character(sd[[2]]), as.character(sd[[1]]))

meta_out_dir <- dirname(opt$metaPrefix)
meta_out_prefix <- sub("_+$", "", basename(opt$metaPrefix))
if (!meta_out_dir %in% c("", ".")) {
  dir.create(meta_out_dir, recursive = TRUE, showWarnings = FALSE)
}

feature_subset <- NULL
feature_flag <- trimws(opt$featureColList)
if (nzchar(feature_flag) && toupper(feature_flag) != "NULL") {
  feature_subset <- strsplit(feature_flag, ",", fixed = TRUE)[[1]]
  feature_subset <- trimws(feature_subset)
  feature_subset <- feature_subset[nzchar(feature_subset)]
  if (length(feature_subset) == 0L) {
    feature_subset <- NULL
  }
}
chrom_flag <- trimws(opt$chrom)
if (!nzchar(chrom_flag) || toupper(chrom_flag) == "NULL") {
  opt$chrom <- NULL
}

metaSummary(
  study_dirs = study_dirs,
  inputPrefix = opt$inputPrefix,
  chrom = opt$chrom,
  featureColList = feature_subset,
  out_dir    = meta_out_dir,
  out_prefix = meta_out_prefix,
  keep_het   = TRUE
)
