#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMmGWAS)
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--abdFile",  type = "character",
              help = "Abundance table path"),
  make_option("--covFile",  type = "character",
              help = "Covariate table path"),
  make_option("--genoFile", type = "character",
              help = "Genotype input: PLINK prefix, VCF(.vcf/.vcf.gz/.vcf.bgz), or BGEN(.bgen)"),
  make_option("--vcfField", type = "character", default = "DS",
              help = "VCF FORMAT field to import for VCF input: DS or GT [default %default]"),
  make_option("--alleleOrder", type = "character", default = "NULL",
              help = "BGEN allele order: ref-first, ref-last, ref-unknown, or NULL [default %default]"),
  make_option("--keepTemp", type = "logical", default = FALSE,
              help = "Keep temporary converted PLINK files for VCF/BGEN input [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$abdFile) || is.null(opt$covFile) || is.null(opt$genoFile)) {
  stop("abdFile, covFile, and genoFile must all be provided.")
}

read_ids_firstcol <- function(path) {
  df <- fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(df) < 1) stop("File has no columns: ", path)
  if (anyNA(df)) stop("Missing values detected in df: ", path)
  ids <- as.character(df[[1]])
  if (anyDuplicated(ids)) stop("Duplicated IID in file: ", path)
  ids
}

read_fam_iid <- function(prefix) {
  fam <- paste0(prefix, ".fam")
  if (!file.exists(fam)) stop("Missing .fam file: ", fam)
  fam_df <- fread(fam, data.table = FALSE, header = FALSE)
  if (ncol(fam_df) < 2) stop("Invalid .fam (need >=2 cols): ", fam)
  ids <- as.character(fam_df[[2]])  # IID
  if (anyDuplicated(ids)) stop("Duplicated IID in .fam: ", fam)
  ids
}

reorder_file_to_ref <- function(path, ref_ids) {
  df <- fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(df) < 1) stop("File has no columns: ", path)

  ids <- as.character(df[[1]])
  ord <- match(ref_ids, ids)
  if (anyNA(ord)) {
    stop("Cannot reorder: some ref IIDs not found in file: ", path)
  }

  df2 <- df[ord, , drop = FALSE]
  # overwrite original
  fwrite(df2, file = path, sep = "\t", quote = FALSE, na = "NA", col.names = TRUE)
}

cat("Checking IIDs...\n")

abd_ids <- read_ids_firstcol(opt$abdFile)
cov_ids <- read_ids_firstcol(opt$covFile)
if (is.null(opt$alleleOrder) || !nzchar(opt$alleleOrder) || toupper(opt$alleleOrder) == "NULL") {
  opt$alleleOrder <- NULL
}

# 1) Strict set check between abundance/covariates
if (!setequal(abd_ids, cov_ids)) stop("IID set mismatch: abdFile vs covFile.")

# 2) Always use genotype sample order as the reference order
changed <- FALSE
geno_input <- PALMmGWAS:::prepare_plink_input(
  genoFile = opt$genoFile,
  vcfField = opt$vcfField,
  alleleOrder = opt$alleleOrder,
  keepTemp = opt$keepTemp,
  tempLabel = "check_tmp"
)
if (!opt$keepTemp && length(geno_input$cleanup) > 0L) {
  on.exit(unlink(geno_input$cleanup, force = TRUE), add = TRUE)
}
geno_ids <- read_fam_iid(geno_input$prefix)
if (!setequal(abd_ids, geno_ids)) stop("IID set mismatch: abdFile vs genoFile (.fam).")

if (!identical(abd_ids, geno_ids)) {
  cat("Order mismatch: reordering abdFile to match genotype IID order...\n")
  reorder_file_to_ref(opt$abdFile, geno_ids)
  changed <- TRUE
}

if (!identical(cov_ids, geno_ids)) {
  cat("Order mismatch: reordering covFile to match genotype IID order...\n")
  reorder_file_to_ref(opt$covFile, geno_ids)
  changed <- TRUE
}

abd_ids2 <- read_ids_firstcol(opt$abdFile)
cov_ids2 <- read_ids_firstcol(opt$covFile)

if (!identical(abd_ids2, geno_ids)) stop("After reorder, abdFile IID order still mismatched.")
if (!identical(cov_ids2, geno_ids)) stop("After reorder, covFile IID order still mismatched.")

if (changed) {
  cat("Done: abd/cov reordered and overwritten to match genotype IID order.\n")
} else {
  cat("ID check passed: abd, cov, and genotype are already identical and aligned.\n")
}
