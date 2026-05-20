#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMGWAS)
  library(optparse)
})

option_list <- list(
    make_option("--abdFile",
        type = "character", default = "",
        help = ""
    ),
    make_option("--covFile",
        type = "character", default = NULL,
        help = ""
    ),
    make_option("--covarColList",
        type = "character", default = NULL,
        help = ""
    ),
    make_option("--depthCol",
        type = "character", default = NULL,
        help = ""
    ),
    make_option("--clusterCol",
        type = "character", default = NULL,
        help = ""
    ),
    make_option("--prev.filter",
        type = "double", default = 0.1,
        help = ""
    ),
    make_option("--nullModelPrefix",
        type = "character", default = "",
        help = ""
    ),
    make_option("--featureInfoFile",
        type = "character", default = NULL,
        help = ""
    )
)

opt <- parse_args(OptionParser(option_list = option_list))
# normalize omitted, empty-string, and literal NULL command-line values to R NULL
if (is.null(opt$covFile) || !nzchar(opt$covFile) || toupper(opt$covFile) == "NULL") {
  opt$covFile <- NULL
}
if (is.null(opt$covarColList) || !nzchar(opt$covarColList) || toupper(opt$covarColList) == "NULL") {
  opt$covarColList <- NULL
}
if (is.null(opt$depthCol) || !nzchar(opt$depthCol) || toupper(opt$depthCol) == "NULL") {
  opt$depthCol <- NULL
}
if (is.null(opt$clusterCol) || !nzchar(opt$clusterCol) || toupper(opt$clusterCol) == "NULL") {
  opt$clusterCol <- NULL
}
if (is.null(opt$featureInfoFile) || !nzchar(opt$featureInfoFile) || toupper(opt$featureInfoFile) == "NULL") {
  opt$featureInfoFile <- NULL
}

fitNULL(
  abdFile = opt$abdFile,
  covFile = opt$covFile,
  covarColList = opt$covarColList,
  depthCol = opt$depthCol,
  clusterCol = opt$clusterCol,
  prev.filter = opt$prev.filter,
  featureInfoFile = opt$featureInfoFile,
  nullModelPrefix = opt$nullModelPrefix
)
