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
    make_option("--outputSnpFile",
        type = "character", default = "NULL",
        help = ""
    ),
    make_option("--correct",
        type = "character", default = "NULL",
        help = ""
    ),
    make_option("--useCluster",
        type = "logical", default = FALSE,
        help = "")
)

opt <- parse_args(OptionParser(option_list = option_list))
# normalize correct from optparse (character) to R NULL
if (is.null(opt$correct) || !nzchar(opt$correct) || toupper(opt$correct) == "NULL") {
  opt$correct <- NULL
}
if (is.null(opt$outputSnpFile) || !nzchar(opt$outputSnpFile) || toupper(opt$outputSnpFile) == "NULL") {
  opt$outputSnpFile <- NULL
}
if (is.null(opt$featureColList) || !nzchar(opt$featureColList) || toupper(opt$featureColList) == "NULL") {
  opt$featureColList <- NULL
}

ptm <- proc.time()
message("step2: PALM summary started.")
message("step2: output prefix = ", opt$PALMOutputFile)
message("step2: chromosome = ", if (is.null(opt$chrom) || !nzchar(opt$chrom)) "NULL" else opt$chrom)
message(
  "step2: featureColList = ",
  if (is.null(opt$featureColList)) {
    "NULL (all modeled features)"
  } else {
    opt$featureColList
  }
)

getSummary(
  genoFile = opt$genoFile,
  NULLmodelFile = opt$NULLmodelFile,
  PALMOutputFile = opt$PALMOutputFile,
  chrom = opt$chrom,
  featureColList = opt$featureColList,
  minMAF = opt$minMAF,
  minMAC = opt$minMAC,
  outputSnpFile = opt$outputSnpFile,
  correct = opt$correct,
  useCluster = opt$useCluster
)

elapsed <- proc.time() - ptm
message(
  sprintf(
    "step2: PALM summary finished. user=%.2fs system=%.2fs elapsed=%.2fs",
    unname(elapsed["user.self"]),
    unname(elapsed["sys.self"]),
    unname(elapsed["elapsed"])
  )
)
