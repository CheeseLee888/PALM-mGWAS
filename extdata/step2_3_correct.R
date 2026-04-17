#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

option_list <- list(
  make_option("--inputPrefix",
    type = "character", default = "",
    help = ""
  ),
  make_option("--outputPrefix",
    type = "character", default = "NULL",
    help = ""
  ),
  make_option("--NULLmodelFile",
    type = "character", default = "NULL",
    help = ""
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$outputPrefix) || !nzchar(opt$outputPrefix) || toupper(opt$outputPrefix) == "NULL") {
  opt$outputPrefix <- opt$inputPrefix
}
if (is.null(opt$NULLmodelFile) || !nzchar(opt$NULLmodelFile) || toupper(opt$NULLmodelFile) == "NULL") {
  opt$NULLmodelFile <- NULL
}

ptm <- proc.time()
message("step2.3: correction started.")
message("step2.3: input prefix = ", opt$inputPrefix)
message("step2.3: output prefix = ", opt$outputPrefix)
message("step2.3: correction mode = median")

correctSummary(
  inputPrefix = opt$inputPrefix,
  outputPrefix = opt$outputPrefix,
  NULLmodelFile = opt$NULLmodelFile
)

elapsed <- proc.time() - ptm
message(
  sprintf(
    "step2.3: correction finished. user=%.2fs system=%.2fs elapsed=%.2fs",
    unname(elapsed["user.self"]),
    unname(elapsed["sys.self"]),
    unname(elapsed["elapsed"])
  )
)
