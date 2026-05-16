#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMGWAS)
  library(optparse)
})

option_list <- list(
    make_option("--genoFile",
        type = "character", default = "",
        help = ""
    ),
    make_option("--nullModelPrefix",
        type = "character", default = "",
        help = ""
    ),
    make_option("--outputPrefix",
        type = "character", default = "",
        help = ""
    ),
    make_option("--chrom",
        type = "character", default = "",
        help = ""
    ),
    make_option("--featureList",
        type = "character", default = "",
        help = ""
    ),
    make_option("--minMAF",
        type = "double", default = 0.05,
        help = ""
    ),
    make_option("--minMAC",
        type = "integer", default = 5,
        help = ""
    ),
    make_option("--maxMissing",
        type = "double", default = 0.15,
        help = ""
    ),
    make_option("--impute_method",
        type = "character", default = "best_guess",
        help = ""
    ),
    make_option("--snpInfoFile",
        type = "character", default = "NULL",
        help = "")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$snpInfoFile) || !nzchar(opt$snpInfoFile) || toupper(opt$snpInfoFile) == "NULL") {
  opt$snpInfoFile <- NULL
}
if (is.null(opt$featureList) || !nzchar(opt$featureList) || toupper(opt$featureList) == "NULL") {
  opt$featureList <- NULL
}

message("step2.1: PALM summary started.")
message("step2.1: output prefix = ", opt$outputPrefix)
message("step2.1: chromosome = ", if (is.null(opt$chrom) || !nzchar(opt$chrom)) "NULL" else opt$chrom)
message(
  "step2.1: featureList = ",
  if (is.null(opt$featureList)) {
    "NULL (all modeled features)"
  } else {
    opt$featureList
  }
)

getSummary(
  genoFile = opt$genoFile,
  nullModelPrefix = opt$nullModelPrefix,
  outputPrefix = opt$outputPrefix,
  chrom = opt$chrom,
  featureList = opt$featureList,
  minMAF = opt$minMAF,
  minMAC = opt$minMAC,
  maxMissing = opt$maxMissing,
  impute_method = opt$impute_method,
  snpInfoFile = opt$snpInfoFile
)
