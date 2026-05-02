#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

option_list <- list(
    make_option("--genoFile",
        type = "character", default = "",
        help = ""
    ),
    make_option("--NULLObjPrefix",
        type = "character", default = "",
        help = ""
    ),
    make_option("--SummaryPrefix",
        type = "character", default = "",
        help = ""
    ),
    make_option("--chrom",
        type = "character", default = "",
        help = ""
    ),
    make_option("--featureColList",
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
    make_option("--SnpInfoFile",
        type = "character", default = "NULL",
        help = ""
    ),
    make_option("--useCluster",
        type = "logical", default = FALSE,
        help = ""),
    make_option("--clusterFile",
        type = "character", default = "NULL",
        help = "")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$SnpInfoFile) || !nzchar(opt$SnpInfoFile) || toupper(opt$SnpInfoFile) == "NULL") {
  opt$SnpInfoFile <- NULL
}
if (is.null(opt$featureColList) || !nzchar(opt$featureColList) || toupper(opt$featureColList) == "NULL") {
  opt$featureColList <- NULL
}
if (is.null(opt$clusterFile) || !nzchar(opt$clusterFile) || toupper(opt$clusterFile) == "NULL") {
  opt$clusterFile <- NULL
}

message("step2.1: PALM summary started.")
message("step2.1: summary prefix = ", opt$SummaryPrefix)
message("step2.1: chromosome = ", if (is.null(opt$chrom) || !nzchar(opt$chrom)) "NULL" else opt$chrom)
message(
  "step2.1: featureColList = ",
  if (is.null(opt$featureColList)) {
    "NULL (all modeled features)"
  } else {
    opt$featureColList
  }
)

getSummary(
  genoFile = opt$genoFile,
  NULLObjPrefix = opt$NULLObjPrefix,
  SummaryPrefix = opt$SummaryPrefix,
  chrom = opt$chrom,
  featureColList = opt$featureColList,
  minMAF = opt$minMAF,
  minMAC = opt$minMAC,
  maxMissing = opt$maxMissing,
  impute_method = opt$impute_method,
  SnpInfoFile = opt$SnpInfoFile,
  useCluster = opt$useCluster,
  clusterFile = opt$clusterFile
)
