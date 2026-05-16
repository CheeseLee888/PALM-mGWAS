#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMGWAS)
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
  make_option("--inputPrefixFile", type="character", default="",
              help="txt: each line 'studyID<TAB>step2Prefix'"),
  make_option("--chrom", type="character", default="NULL",
              help="Step2 scope: NULL for allchr, or 1..22 [default %default]"),
  make_option("--featureList", type="character", default="NULL",
              help="Optional comma-separated feature name(s) to meta-analyze [default %default]"),
  make_option("--meta.method", type="character", default="EE",
              help="Meta-analysis method passed to metafor::rma.uni() [default %default]"),
  make_option("--outputPrefix", type="character", default="",
              help="Full output prefix for meta files, e.g. example/output/meta/step3_meta")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!nzchar(opt$outputPrefix))
  stop("Missing --outputPrefix")

if (!nzchar(opt$inputPrefixFile) || !file.exists(opt$inputPrefixFile)) {
  stop("Missing/invalid --inputPrefixFile")
}
sd <- read.table(opt$inputPrefixFile, header = FALSE, sep = "", stringsAsFactors = FALSE)
if (ncol(sd) < 2) stop("inputPrefixFile must have >=2 columns: studyID and Step2 prefix")
input_prefixes <- setNames(as.character(sd[[2]]), as.character(sd[[1]]))
study_dirs <- setNames(dirname(input_prefixes), names(input_prefixes))

meta_out_dir <- dirname(opt$outputPrefix)
meta_out_prefix <- sub("_+$", "", basename(opt$outputPrefix))
if (!meta_out_dir %in% c("", ".")) {
  dir.create(meta_out_dir, recursive = TRUE, showWarnings = FALSE)
}

feature_subset <- NULL
feature_flag <- trimws(opt$featureList)
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
  inputPrefix = input_prefixes,
  chrom = opt$chrom,
  featureList = feature_subset,
  out_dir    = meta_out_dir,
  out_prefix = meta_out_prefix,
  keep_het   = TRUE,
  meta.method = opt$meta.method
)
