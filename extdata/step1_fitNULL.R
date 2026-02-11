#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

options(stringsAsFactors = F)

## load R libraries

library(SAIGEQTL)

require(optparse) # install.packages("optparse")

print(sessionInfo())

## set list of cmd line arguments
option_list <- list(
  make_option("--grmFile",
    type = "character", default = "",
    help = "Path to save the GRM RDS file"
  ),
  make_option("--abdFile",
    type = "character", default = "",
    help = "Required. Path to the abundance file."
  ),
  make_option("--covFile",
    type = "character", default = "",
    help = "Required. Path to the covariate file."
  ), 
  make_option("--covarColList",
    type = "character", default = "ALL",
    help = "List of covariates (comma separated)"
  ),
  make_option("--outputPrefix",
    type = "character", default = "",
    help = "Required. Path and prefix of the output files [default='']"
  ),
  make_option("--isCovariateOffset",
    type = "logical", default = TRUE,
    help = "Optional. Whether to estimate fixed effect coeffciets. [default, 'TRUE']"
  ),
  make_option("--useGRMtoFitNULL", type = "logical", default = TRUE, help = ""),
  make_option("--batch_idx",
    type = "integer", default = 1,
    help = "Index of the batch to process (starting from 1) [default=1]"
  )
)


## list of options
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

## covariates dealing
if (is.null(opt$covarColList) || length(opt$covarColList) == 0 || toupper(opt$covarColList)=="NULL") {
  cat("No covariates specified\n")
  covars <- character(0)
} else if (toupper(opt$covarColList)=="ALL") {
  cat("Using all covariates from covFile header\n")
  cov <- read_table_with_id(opt$covFile)
  cov_header <- colnames(cov)
  ## Remove ID columns
  drop_cols <- c("IID")
  covars <- setdiff(cov_header, drop_cols)
  if (!length(covars)) {
    stop("No covariate columns found in covFile after removing ID and offset columns. Please set '--covarColList= NULL' and try again.\n")
  }
} else {
  cat("Using user-specified covariates from --covarColList\n")
  covars <- strsplit(opt$covarColList, ",")[[1]]
  covars <- trimws(covars)
  covars <- covars[nzchar(covars)]
}

# convertoNumeric <- function(x, stringOutput) {
#   y <- tryCatch(expr = as.numeric(x), warning = function(w) {
#     return(NULL)
#   })
#   if (is.null(y)) {
#     stop(stringOutput, " is not numeric\n")
#   } else {
#     cat(stringOutput, " is ", y, "\n")
#   }
#   return(y)
# }

if (basename(toupper(opt$grmFile)) == "NULL") {
  opt$useGRMtoFitNULL <- FALSE
}

print(opt)

set.seed(1)
fitNULLGLMM_multiV(
  grmFile = opt$grmFile,
  abdFile = opt$abdFile,
  covFile = opt$covFile,
  covarColList = covars,
  outputPrefix = opt$outputPrefix,
  isCovariateOffset = opt$isCovariateOffset,
  useGRMtoFitNULL = opt$useGRMtoFitNULL,
  batch_idx = as.integer(opt$batch_idx)
)
