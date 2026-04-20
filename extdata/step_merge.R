#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
  local_impl <- file.path(dirname(script_path), "..", "R", "merge.R")
  if (file.exists(local_impl)) {
    source(local_impl)
  }
}


option_list <- list(
  make_option("--inputPrefix",
    type = "character", default = "",
    help = ""
  ),
  make_option("--featureColList",
    type = "character", default = "",
    help = ""
  ),
  make_option("--outputPrefix",
    type = "character", default = "NULL",
    help = ""
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$featureColList) || !nzchar(opt$featureColList) || toupper(opt$featureColList) == "NULL") {
  opt$featureColList <- NULL
}
if (is.null(opt$outputPrefix) || !nzchar(opt$outputPrefix) || toupper(opt$outputPrefix) == "NULL") {
  opt$outputPrefix <- NULL
}

message("merge: started.")
message("merge: input prefix = ", opt$inputPrefix)
message(
  "merge: featureColList = ",
  if (is.null(opt$featureColList)) {
    "NULL (all features)"
  } else {
    opt$featureColList
  }
)
message("merge: output prefix = ", if (is.null(opt$outputPrefix)) "NULL (auto)" else opt$outputPrefix)

mergeResults(
  inputPrefix = opt$inputPrefix,
  featureColList = opt$featureColList,
  outputPrefix = opt$outputPrefix
)

message("merge: finished.")
