#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

option_list <- list(
    make_option("--genoPrefix",
        type = "character", default = "",
        help = ""
    ),
    make_option("--NULLmodelFile",
        type = "character", default = "",
        help = ""
    ),
    make_option("--PALMOutputFile",
        type = "character", default = "",
        help = ""
    ),
    make_option("--chrom",
        type = "character", default = "",
        help = ""
    ),
    make_option("--correct",
        type = "character", default = "NULL",
        help = ""
    ),
    make_option("--useCluster",
        type = "logical", default = TRUE,
        help = "")
)

opt <- parse_args(OptionParser(option_list = option_list))
# normalize correct from optparse (character) to R NULL
if (is.null(opt$correct) || !nzchar(opt$correct) || toupper(opt$correct) == "NULL") {
  opt$correct <- NULL
}

getSummary(
  genoPrefix = opt$genoPrefix,
  NULLmodelFile = opt$NULLmodelFile,
  PALMOutputFile = opt$PALMOutputFile,
  chrom = opt$chrom,
  correct = opt$correct,
  useCluster = opt$useCluster
)
