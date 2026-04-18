#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

option_list <- list(
    make_option("--abdFile",
        type = "character", default = "",
        help = ""
    ),
    make_option("--covFile",
        type = "character", default = "",
        help = ""
    ),
    make_option("--covarColList",
        type = "character", default = "",
        help = ""
    ),
    make_option("--depthCol",
        type = "character", default = "",
        help = ""
    ),
    make_option("--prev.filter",
        type = "double", default = 0.1,
        help = ""
    ),
    make_option("--NULLObjPrefix",
        type = "character", default = "",
        help = ""
    ),
    make_option("--FeatureInfoFile",
        type = "character", default = "NULL",
        help = ""
    ),
    make_option("--FeatureNameListFile",
        type = "character", default = "NULL",
        help = ""
    )
)

opt <- parse_args(OptionParser(option_list = option_list))
# normalize covFile from optparse (character) to R NULL
if (is.null(opt$covFile) || !nzchar(opt$covFile) || toupper(opt$covFile) == "NULL") {
  opt$covFile <- NULL
}
if (is.null(opt$covarColList) || !nzchar(opt$covarColList) || toupper(opt$covarColList) == "NULL") {
  opt$covarColList <- NULL
}
if (is.null(opt$depthCol) || !nzchar(opt$depthCol) || toupper(opt$depthCol) == "NULL") {
  opt$depthCol <- NULL
}
if (is.null(opt$FeatureInfoFile) || !nzchar(opt$FeatureInfoFile) || toupper(opt$FeatureInfoFile) == "NULL") {
  opt$FeatureInfoFile <- NULL
}
if (is.null(opt$FeatureNameListFile) || !nzchar(opt$FeatureNameListFile) || toupper(opt$FeatureNameListFile) == "NULL") {
  opt$FeatureNameListFile <- NULL
}

fitNULL(
  abdFile = opt$abdFile,
  covFile = opt$covFile,
  covarColList = opt$covarColList,
  depthCol = opt$depthCol,
  prev.filter = opt$prev.filter,
  FeatureInfoFile = opt$FeatureInfoFile,
  FeatureNameListFile = opt$FeatureNameListFile,
  NULLObjPrefix = opt$NULLObjPrefix
)
