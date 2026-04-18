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
  make_option("--overwriteOutput",
    type = "character", default = "TRUE",
    help = ""
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
overwrite_flag <- toupper(trimws(opt$overwriteOutput))
if (!overwrite_flag %in% c("TRUE", "FALSE")) {
  stop("--overwriteOutput must be TRUE or FALSE.")
}
opt$overwriteOutput <- identical(overwrite_flag, "TRUE")

ptm <- proc.time()
message("step2.3: correction started.")
message("step2.3: input prefix = ", opt$inputPrefix)
message("step2.3: overwrite output = ", opt$overwriteOutput)
message("step2.3: correction mode = median")

correctSummary(
  inputPrefix = opt$inputPrefix,
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
