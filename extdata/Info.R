#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
})

option_list <- list(
  make_option("--genoFile", type = "character", default = NULL,
              help = "Genotype input: PLINK prefix, VCF(.vcf/.vcf.gz/.vcf.bgz), or BGEN(.bgen)"),
  make_option("--abdFile", type = "character", default = NULL,
              help = "Path to the microbiome abundance data file"),
  make_option("--vcfField", type = "character", default = "DS",
              help = "VCF FORMAT field to import for VCF input: DS or GT [default %default]"),
  make_option("--alleleOrder", type = "character", default = "NULL",
              help = "BGEN allele order: ref-first, ref-last, ref-unknown, or NULL [default %default]"),
  make_option("--keepTemp", type = "logical", default = FALSE,
              help = "Keep temporary converted PLINK files for VCF/BGEN input [default %default]"),
  make_option("--outputSnpFile", type = "character", default = NULL,
              help = "Output path for SNP info"),
  make_option("--outputFeatureFile", type = "character", default = NULL,
              help = "Output path for feature info"),
  make_option("--outputSeqDepthFile", type = "character", default = NULL,
              help = "Output path for sequencing depth info")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$genoFile) || !nzchar(opt$genoFile) || is.null(opt$abdFile) ||
    is.null(opt$outputSnpFile) || is.null(opt$outputFeatureFile) || is.null(opt$outputSeqDepthFile)) {
  stop("All arguments (--genoFile, --abdFile, --outputSnpFile, --outputFeatureFile, --outputSeqDepthFile) are required.")
}

if (is.null(opt$alleleOrder) || !nzchar(opt$alleleOrder) || toupper(opt$alleleOrder) == "NULL") {
  opt$alleleOrder <- NULL
}

run_info(
  genoFile = opt$genoFile,
  abd_file = opt$abdFile,
  output_snp = opt$outputSnpFile,
  output_feature = opt$outputFeatureFile,
  output_seqdepth = opt$outputSeqDepthFile,
  vcf_field = opt$vcfField,
  allele_order = opt$alleleOrder,
  keep_temp = opt$keepTemp
)

cat("Finished generating info files:\n")
cat("  SNP info: ", opt$outputSnpFile, "\n", sep = "")
cat("  Feature info: ", opt$outputFeatureFile, "\n", sep = "")
cat("  SeqDepth info: ", opt$outputSeqDepthFile, "\n", sep = "")
