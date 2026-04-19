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

ptm <- proc.time()
message("step2.3: correction started.")
message("step2.3: input prefix = ", opt$inputPrefix)
message("step2.3: chromosome scope = ", if (is.null(opt$chrom)) "NULL" else opt$chrom)
message("step2.3: overwrite output = ", opt$overwriteOutput)
message("step2.3: correction mode = median")

correctSummary(
  inputPrefix = opt$inputPrefix,
  chrom = opt$chrom,
  overwriteOutput = opt$overwriteOutput
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
