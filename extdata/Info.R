#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

option_list <- list(
  make_option("--genoFile", type = "character", default = NULL,
              help = "PLINK prefix (without .bed/.bim/.fam)"),
  make_option("--abdFile", type = "character", default = NULL,
              help = "Path to the microbiome abundance data file"),
  make_option("--outputSnpFile", type = "character", default = NULL,
              help = "Output path for SNP info"),
  make_option("--outputFeatureFile", type = "character", default = NULL,
              help = "Output path for feature info"),
  make_option("--outputSeqDepthFile", type = "character", default = NULL,
              help = "Output path for sequencing depth info")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$genoFile) || is.null(opt$abdFile) ||
    is.null(opt$outputSnpFile) || is.null(opt$outputFeatureFile) || is.null(opt$outputSeqDepthFile)) {
  stop("All arguments (--genoFile, --abdFile, --outputSnpFile, --outputFeatureFile, --outputSeqDepthFile) are required.")
}

run_info(
  geno_file = opt$genoFile,
  abd_file = opt$abdFile,
  output_snp = opt$outputSnpFile,
  output_feature = opt$outputFeatureFile,
  output_seqdepth = opt$outputSeqDepthFile
)

cat("Finished generating info files:\n")
cat("  SNP info: ", opt$outputSnpFile, "\n", sep = "")
cat("  Feature info: ", opt$outputFeatureFile, "\n", sep = "")
cat("  SeqDepth info: ", opt$outputSeqDepthFile, "\n", sep = "")
