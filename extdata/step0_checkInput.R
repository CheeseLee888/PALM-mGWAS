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
  make_option("--abdAlignedFile", type = "character",
              help = "Output path for aligned abundance table"),
  make_option("--covAlignedFile", type = "character",
              help = "Output path for aligned covariate table"),
  make_option("--covarColList", type = "character", default = "NULL",
              help = "Optional comma-separated covariate columns used in Step1; if NULL, all non-ID covariate columns are used. Samples missing these columns are removed [default %default]"),
  make_option("--depthCol", type = "character", default = "NULL",
              help = "Optional covariate column name used as sequencing depth [default %default]"),
  make_option("--depth.filter", type = "double", default = 0,
              help = "Sample-level depth threshold; samples with depth <= threshold are removed before ID matching [default %default]"),
  make_option("--genoFile", type = "character",
              help = "Genotype input: PLINK prefix, VCF(.vcf/.vcf.gz/.vcf.bgz), or BGEN(.bgen)"),
  make_option("--vcfField", type = "character", default = "DS",
              help = "VCF FORMAT field to import for VCF input: DS or GT [default %default]"),
  make_option("--alleleOrder", type = "character", default = "NULL",
              help = "BGEN allele order: ref-first, ref-last, ref-unknown, or NULL [default %default]"),
  make_option("--keepTemp", type = "logical", default = FALSE,
              help = "Keep temporary converted PLINK files for VCF/BGEN input [default %default]"),
  make_option("--outputSeqDepthFile", type = "character", default = "NULL",
              help = "Optional output file for sequencing depth info used by Step0 filtering [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$abdFile) || is.null(opt$covFile) || is.null(opt$genoFile) ||
    is.null(opt$abdAlignedFile) || is.null(opt$covAlignedFile)) {
  stop("abdFile, covFile, genoFile, abdAlignedFile, and covAlignedFile must all be provided.")
}

normalize_col_list <- function(x, arg_name) {
  if (is.null(x) || !nzchar(x) || toupper(x) == "NULL") {
    return(NULL)
  }
  cols <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  cols <- cols[nzchar(cols)]
  if (!length(cols)) {
    stop("'", arg_name, "' did not contain any valid column names.")
  }
  cols
}

default_covariate_cols <- function(df) {
  cols <- colnames(df)
  if (length(cols) <= 1L) {
    return(character(0))
  }
  cols[-1]
}

read_abd_table <- function(path) {
  df <- fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(df) < 1) stop("File has no columns: ", path)
  if (anyNA(df)) stop("Missing values detected in abundance file: ", path)
  ids <- as.character(df[[1]])
  if (anyNA(ids) || any(!nzchar(ids))) stop("Missing/empty IID detected in abundance file: ", path)
  if (anyDuplicated(ids)) stop("Duplicated IID in file: ", path)
  df
}

read_cov_table <- function(path) {
  df <- fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(df) < 1) stop("File has no columns: ", path)
  ids <- as.character(df[[1]])
  if (anyNA(ids) || any(!nzchar(ids))) stop("Missing/empty IID detected in covariate file: ", path)
  if (anyDuplicated(ids)) stop("Duplicated IID in file: ", path)
  df
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

reorder_df_to_ref <- function(df, ref_ids, path_label) {
  ids <- as.character(df[[1]])
  ord <- match(ref_ids, ids)
  if (anyNA(ord)) {
    stop("Cannot reorder: some ref IIDs not found in file: ", path_label)
  }
  df[ord, , drop = FALSE]
}

cat("Checking IIDs...\n")

abd_df <- read_abd_table(opt$abdFile)
cov_df <- read_cov_table(opt$covFile)
cat(
  "Loaded input tables: abd=", nrow(abd_df),
  " sample(s), cov=", nrow(cov_df), " sample(s).\n",
  sep = ""
)
if (is.null(opt$alleleOrder) || !nzchar(opt$alleleOrder) || toupper(opt$alleleOrder) == "NULL") {
  opt$alleleOrder <- NULL
}
if (is.null(opt$outputSeqDepthFile) || !nzchar(opt$outputSeqDepthFile) || toupper(opt$outputSeqDepthFile) == "NULL") {
  opt$outputSeqDepthFile <- NULL
}
opt$covarColList <- normalize_col_list(opt$covarColList, "covarColList")
opt$depthCol <- normalize_col_list(opt$depthCol, "depthCol")
if (!is.null(opt$depthCol) && length(opt$depthCol) != 1L) {
  stop("'depthCol' must specify exactly one column name.")
}
if (!is.numeric(opt$depth.filter) || length(opt$depth.filter) != 1L || is.na(opt$depth.filter) || opt$depth.filter < 0) {
  stop("--depth.filter must be a single non-negative numeric value.")
}

if (is.null(opt$covarColList)) {
  opt$covarColList <- default_covariate_cols(cov_df)
  if (length(opt$covarColList) > 0L) {
    cat("covarColList not provided: defaulting to all covariate columns in covFile.\n")
  } else {
    cat("covarColList not provided: covFile has no non-ID covariate columns.\n")
  }
}

required_cov_cols <- unique(c(opt$covarColList, opt$depthCol))
missing_cov_cols <- setdiff(required_cov_cols, colnames(cov_df))
if (length(missing_cov_cols) > 0L) {
  stop("Required covariate column(s) not found in covFile: ", paste(missing_cov_cols, collapse = ", "))
}
if (length(required_cov_cols) > 0L) {
  cat("Step1-required covariate/depth columns: ", paste(required_cov_cols, collapse = ", "), "\n", sep = "")
  missing_counts <- colSums(is.na(cov_df[, required_cov_cols, drop = FALSE]))
  missing_summary <- paste(names(missing_counts), missing_counts, sep = "=")
  cat("Missing counts by required covariate/depth column: ", paste(missing_summary, collapse = ", "), "\n", sep = "")
} else {
  cat("No Step1-required covariate/depth columns requested for missingness filtering.\n")
}

filtered <- FALSE
abd_ids_all <- as.character(abd_df[[1]])
if (is.null(opt$depthCol)) {
  depth <- rowSums(as.matrix(abd_df[, -1, drop = FALSE]), na.rm = TRUE)
  names(depth) <- abd_ids_all
  cat("Depth filter source: row sums of abdFile.\n")
} else {
  depth <- suppressWarnings(as.numeric(cov_df[[opt$depthCol]]))
  if (any(is.na(depth) & !is.na(cov_df[[opt$depthCol]]))) {
    stop("Requested depth column in covFile cannot be safely converted to numeric: ", opt$depthCol)
  }
  names(depth) <- as.character(cov_df[[1]])
  cat("Depth filter source: covFile column '", opt$depthCol, "'.\n", sep = "")
}
cat(
  "Depth summary before filtering: min=", min(depth),
  ", median=", stats::median(depth),
  ", max=", max(depth), ".\n",
  sep = ""
)
if (!is.null(opt$outputSeqDepthFile)) {
  cat("Generating DepthInfo from the exact depth values used by Step0 sample filtering...\n")
  seqdepth_df <- PALMmGWAS:::seqdepth_info_from_values(names(depth), depth)
  dir.create(dirname(opt$outputSeqDepthFile), recursive = TRUE, showWarnings = FALSE)
  fwrite(seqdepth_df, file = opt$outputSeqDepthFile, sep = "\t", quote = FALSE, na = "NA", col.names = TRUE)
  cat(
    "DepthInfo finished: ", nrow(seqdepth_df),
    " sample(s) written to ", opt$outputSeqDepthFile, ".\n",
    sep = ""
  )
} else {
  cat("DepthInfo skipped: --outputSeqDepthFile is NULL.\n")
}

if (opt$depth.filter > 0) {
  keep_ids <- names(depth)[depth > opt$depth.filter]
  removed_n <- length(depth) - length(keep_ids)
  cat("Applying depth.filter=", opt$depth.filter, ": keeping ", length(keep_ids), " sample(s), removing ", removed_n, ".\n", sep = "")
  if (!length(keep_ids)) {
    stop("No samples remain after applying --depth.filter=", opt$depth.filter)
  }
  abd_df <- abd_df[as.character(abd_df[[1]]) %in% keep_ids, , drop = FALSE]
  cov_df <- cov_df[as.character(cov_df[[1]]) %in% keep_ids, , drop = FALSE]
  filtered <- filtered || removed_n > 0L
}

if (length(required_cov_cols) > 0L) {
  keep_cov_complete <- stats::complete.cases(cov_df[, required_cov_cols, drop = FALSE])
  removed_n <- sum(!keep_cov_complete)
  cat("Removing ", removed_n, " sample(s) with missing values in Step1-required covariates/depth columns.\n", sep = "")
  if (removed_n > 0L) {
    keep_ids <- as.character(cov_df[[1]])[keep_cov_complete]
    cov_df <- cov_df[keep_cov_complete, , drop = FALSE]
    abd_df <- abd_df[as.character(abd_df[[1]]) %in% keep_ids, , drop = FALSE]
    filtered <- TRUE
  }
}

abd_ids <- as.character(abd_df[[1]])
cov_ids <- as.character(cov_df[[1]])
cat(
  "After sample-level filtering: abd=", length(abd_ids),
  " sample(s), cov=", length(cov_ids), " sample(s).\n",
  sep = ""
)

# 1) Keep matched IDs across abundance and covariates after sample filtering
common_abd_cov_ids <- intersect(abd_ids, cov_ids)
if (!length(common_abd_cov_ids)) {
  stop("No matched IID remain across abdFile and covFile after filtering.")
}
if (length(common_abd_cov_ids) < length(abd_ids) || length(common_abd_cov_ids) < length(cov_ids)) {
  cat(
    "Keeping ", length(common_abd_cov_ids),
    " matched sample(s) shared by abdFile and covFile after filtering.\n",
    sep = ""
  )
  abd_df <- abd_df[abd_ids %in% common_abd_cov_ids, , drop = FALSE]
  cov_df <- cov_df[cov_ids %in% common_abd_cov_ids, , drop = FALSE]
  filtered <- TRUE
}
abd_ids <- as.character(abd_df[[1]])
cov_ids <- as.character(cov_df[[1]])
cat("Matched abd/cov sample count: ", length(abd_ids), ".\n", sep = "")

# 2) Use genotype sample order as the reference order; genotype may contain extra samples
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
cat("Genotype sample count: ", length(geno_ids), ".\n", sep = "")
missing_in_geno <- setdiff(abd_ids, geno_ids)
if (length(missing_in_geno) > 0) {
  stop("IID set mismatch: filtered abd/cov samples missing from genoFile (.fam): ",
       paste(utils::head(missing_in_geno, 5), collapse = ", "))
}

ref_ids <- geno_ids[geno_ids %in% abd_ids]
cat(
  "Retained samples present in genotype: ", length(ref_ids),
  ". Extra genotype-only samples: ", length(geno_ids) - length(ref_ids), ".\n",
  sep = ""
)

if (!identical(abd_ids, ref_ids)) {
  cat("Order mismatch: reordering abdFile to match filtered genotype IID order...\n")
  abd_df <- reorder_df_to_ref(abd_df, ref_ids, opt$abdFile)
  changed <- TRUE
}

if (!identical(cov_ids, ref_ids)) {
  cat("Order mismatch: reordering covFile to match filtered genotype IID order...\n")
  cov_df <- reorder_df_to_ref(cov_df, ref_ids, opt$covFile)
  changed <- TRUE
}

abd_ids2 <- as.character(abd_df[[1]])
cov_ids2 <- as.character(cov_df[[1]])

if (!identical(abd_ids2, ref_ids)) stop("After reorder, abdFile IID order still mismatched.")
if (!identical(cov_ids2, ref_ids)) stop("After reorder, covFile IID order still mismatched.")

dir.create(dirname(opt$abdAlignedFile), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$covAlignedFile), recursive = TRUE, showWarnings = FALSE)
fwrite(abd_df, file = opt$abdAlignedFile, sep = "\t", quote = FALSE, na = "NA", col.names = TRUE)
fwrite(cov_df, file = opt$covAlignedFile, sep = "\t", quote = FALSE, na = "NA", col.names = TRUE)
cat(
  "Final aligned sample count: abd=", nrow(abd_df),
  ", cov=", nrow(cov_df), ".\n",
  sep = ""
)

if (changed || filtered) {
  cat("Done: wrote filtered/aligned abd file to ", opt$abdAlignedFile, "\n", sep = "")
  cat("Done: wrote filtered/aligned cov file to ", opt$covAlignedFile, "\n", sep = "")
} else {
  cat("ID check passed: wrote aligned copies without modifying original input files.\n")
}
