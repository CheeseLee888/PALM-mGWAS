#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
  local_impl <- file.path(dirname(script_path), "..", "R", "mergeSummary.R")
  if (file.exists(local_impl)) {
    source(local_impl)
  }
}


option_list <- list(
  make_option("--inputPrefix",
    type = "character", default = "",
    help = ""
  ),
  make_option("--outputPrefix",
    type = "character", default = "NULL",
    help = ""
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$outputPrefix) || !nzchar(opt$outputPrefix) || toupper(opt$outputPrefix) == "NULL") {
  opt$outputPrefix <- NULL
}

message("step2.2: merge started.")
message("step2.2: input prefix = ", opt$inputPrefix)
message("step2.2: output prefix = ", if (is.null(opt$outputPrefix)) "NULL (auto)" else opt$outputPrefix)

mergeSummary(
  inputPrefix = opt$inputPrefix,
  outputPrefix = opt$outputPrefix
)

message("step2.2: merge finished.")
