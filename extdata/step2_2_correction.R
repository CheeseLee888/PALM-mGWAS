#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
  local_impl <- file.path(dirname(script_path), "..", "R", "correctSummary.R")
  if (file.exists(local_impl)) {
    source(local_impl)
  }
}

option_list <- list(
  make_option("--inputPrefix",
    type = "character", default = "",
    help = ""
  ),
  make_option("--chrom",
    type = "character", default = "NULL",
    help = ""
  ),
  make_option("--overwriteOutput",
    type = "character", default = "TRUE",
    help = ""
  ),
  make_option("--correct",
    type = "character", default = "median",
    help = ""
  ),
  make_option("--NULLmodelFile",
    type = "character", default = "NULL",
    help = ""
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
chrom_flag <- trimws(opt$chrom)
if (!nzchar(chrom_flag) || toupper(chrom_flag) == "NULL") {
  opt$chrom <- NULL
}
overwrite_flag <- toupper(trimws(opt$overwriteOutput))
if (!overwrite_flag %in% c("TRUE", "FALSE")) {
  stop("--overwriteOutput must be TRUE or FALSE.")
}
opt$overwriteOutput <- identical(overwrite_flag, "TRUE")
opt$correct <- tolower(trimws(opt$correct))
if (!opt$correct %in% c("median", "tune")) {
  stop("--correct must be median or tune.")
}
if (is.null(opt$NULLmodelFile) || !nzchar(opt$NULLmodelFile) || toupper(opt$NULLmodelFile) == "NULL") {
  opt$NULLmodelFile <- NULL
}

ptm <- proc.time()
message("step2.2: correction started.")
message("step2.2: input prefix = ", opt$inputPrefix)
message("step2.2: chromosome scope = ", if (is.null(opt$chrom)) "NULL" else opt$chrom)
message("step2.2: overwrite output = ", opt$overwriteOutput)
message("step2.2: correction mode = ", opt$correct)
message("step2.2: NULL model file = ", if (is.null(opt$NULLmodelFile)) "NULL" else opt$NULLmodelFile)

correctSummary(
  inputPrefix = opt$inputPrefix,
  chrom = opt$chrom,
  overwriteOutput = opt$overwriteOutput,
  correct = opt$correct,
  NULLmodelFile = opt$NULLmodelFile
)

elapsed <- proc.time() - ptm
message(
  sprintf(
    "step2.2: correction finished. user=%.2fs system=%.2fs elapsed=%.2fs",
    unname(elapsed["user.self"]),
    unname(elapsed["sys.self"]),
    unname(elapsed["elapsed"])
  )
)
