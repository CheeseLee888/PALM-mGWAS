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
    make_option("--outputPrefix",
        type = "character", default = "",
        help = ""
    )
)

opt <- parse_args(OptionParser(option_list = option_list))
# normalize covFile from optparse (character) to R NULL
if (is.null(opt$covFile) || !nzchar(opt$covFile) || toupper(opt$covFile) == "NULL") {
  opt$covFile <- NULL
}

fitNULL(
  abdFile = opt$abdFile,
  covFile = opt$covFile,
  outputPrefix = opt$outputPrefix
)
