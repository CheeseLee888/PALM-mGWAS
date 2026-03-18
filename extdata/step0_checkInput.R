#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--abdFile",  type = "character"),
  make_option("--covFile",  type = "character"),
  make_option("--genoFile", type = "character")
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

infer_geno_format <- function(path) {
  lower <- tolower(path)
  if (grepl("\\.(vcf|vcf\\.gz|vcf\\.bgz)$", lower)) {
    return("vcf")
  }
  if (grepl("\\.bgen$", lower)) {
    return("bgen")
  }
  "plink"
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
geno_format <- infer_geno_format(opt$genoFile)

# 1) Strict set check between abundance/covariates
if (!setequal(abd_ids, cov_ids)) stop("IID set mismatch: abdFile vs covFile.")

# 2) If PLINK input, use .fam as the reference order; otherwise keep abd order
changed <- FALSE

if (identical(geno_format, "plink")) {
  geno_ids <- read_fam_iid(opt$genoFile)
  if (!setequal(abd_ids, geno_ids)) stop("IID set mismatch: abdFile vs genoFile (.fam).")

  if (!identical(abd_ids, geno_ids)) {
    cat("Order mismatch: reordering abdFile to match .fam IID order...\n")
    reorder_file_to_ref(opt$abdFile, geno_ids)
    changed <- TRUE
  }

  if (!identical(cov_ids, geno_ids)) {
    cat("Order mismatch: reordering covFile to match .fam IID order...\n")
    reorder_file_to_ref(opt$covFile, geno_ids)
    changed <- TRUE
  }

  abd_ids2 <- read_ids_firstcol(opt$abdFile)
  cov_ids2 <- read_ids_firstcol(opt$covFile)

  if (!identical(abd_ids2, geno_ids)) stop("After reorder, abdFile IID order still mismatched.")
  if (!identical(cov_ids2, geno_ids)) stop("After reorder, covFile IID order still mismatched.")

  if (changed) {
    cat("Done: abd/cov reordered and overwritten to match geno (.fam) IID order.\n")
  } else {
    cat("ID check passed: abd, cov, geno (.fam) are already identical and aligned.\n")
  }
} else {
  if (!identical(abd_ids, cov_ids)) {
    cat("Order mismatch: reordering covFile to match abdFile IID order...\n")
    reorder_file_to_ref(opt$covFile, abd_ids)
    changed <- TRUE
  }

  cov_ids2 <- read_ids_firstcol(opt$covFile)
  if (!identical(abd_ids, cov_ids2)) stop("After reorder, covFile IID order still mismatched with abdFile.")

  if (changed) {
    cat("Done: covFile reordered to match abdFile IID order. Genotype order will be handled in step2 conversion.\n")
  } else {
    cat("ID check passed: abd and cov are already identical and aligned. Genotype order will be handled in step2 conversion.\n")
  }
}
