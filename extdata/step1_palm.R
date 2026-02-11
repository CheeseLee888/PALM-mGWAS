#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(PALM)
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

read_firstcol_as_rownames <- function(file) {
  df <- data.table::fread(
    file = file,
    data.table   = FALSE,
    check.names  = FALSE
  )
  id <- as.character(df[[1]])
  rownames(df) <- id
  df[[1]] <- NULL
  return(df)
}

abd <- read_firstcol_as_rownames(opt$abdFile)
abd <- as.matrix(abd)

if(is.null(opt$covFile)) {
  # depth is calculated as the row sums of rel.abd
  cat("Fitting PALM null model without covariates.\n")
  modglmm <- palm.null.model(
    rel.abd = abd,
    prev.filter = 0
    )
}else {
  cov <- read_firstcol_as_rownames(opt$covFile)
  modglmm <- palm.null.model(
    rel.abd = abd,
    covariate.adjust = cov,
    prev.filter = 0
    )
}

save(modglmm, file = paste0(opt$outputPrefix, ".rda"))
cat("Done. PALM null model saved to", paste0(opt$outputPrefix, ".rda"), "\n")

