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
    make_option("--depthFile",
        type = "character", default = "",
        help = ""
    ),
    make_option("--depth.filter",
        type = "double", default = 0,
        help = ""
    ),
    make_option("--prev.filter",
        type = "double", default = 0.1,
        help = ""
    ),
    make_option("--NULLmodelFile",
        type = "character", default = "",
        help = ""
    )
)

opt <- parse_args(OptionParser(option_list = option_list))
# normalize covFile from optparse (character) to R NULL
if (is.null(opt$covFile) || !nzchar(opt$covFile) || toupper(opt$covFile) == "NULL") {
  opt$covFile <- NULL
}
if (is.null(opt$depthFile) || !nzchar(opt$depthFile) || toupper(opt$depthFile) == "NULL") {
  opt$depthFile <- NULL
}

fitNULL(
  abdFile = opt$abdFile,
  covFile = opt$covFile,
  depthFile = opt$depthFile,
  depth.filter = opt$depth.filter,
  prev.filter = opt$prev.filter,
  NULLmodelFile = opt$NULLmodelFile
)
