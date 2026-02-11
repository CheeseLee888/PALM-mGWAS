#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

options(stringsAsFactors = F)

## load R libraries

library(SAIGEQTL)

require(optparse) # install.packages("optparse")

print(sessionInfo())

## command line options
option_list <- list(
  make_option(
    "--genoFile",
    type    = "character", default = "",
    help    = "Required. Genotype input: PLINK prefix (no extension) or a .gds file path."
  ),
  make_option(
    "--grmFile",
    type    = "character", default = "",
    help    = "Required. Output RDS file for GRM (list with K and sample.id)."
  )
)

## list of options
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

if (!nzchar(opt$genoFile)) stop("Please specify --genoFile (PLINK prefix or .gds).")
if (!nzchar(opt$grmFile))  stop("Please specify --grmFile (output RDS path).")


generateGRM(
  genoFile = opt$genoFile,
  grmFile  = opt$grmFile
)
