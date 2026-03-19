#!/usr/bin/env Rscript

suppressPackageStartupMessages({
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

find_executable <- function(cmd) {
  path <- Sys.which(cmd)
  if (!nzchar(path)) {
    return("")
  }
  path
}

run_command_checked <- function(command, args) {
  cat("Running command:", command, paste(args, collapse = " "), "\n")
  status <- system2(command, args = args)
  if (!identical(status, 0L)) {
    stop("Command failed with exit status ", status, ": ", command)
  }
}

sanitize_local_name <- function(x) {
  x <- as.character(x)[1]
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_+", "", x)
  x <- sub("_+$", "", x)
  if (!nzchar(x)) {
    x <- "geno"
  }
  x
}

prepare_plink_input <- function(genoFile,
                                genoFormat,
                                vcfField,
                                alleleOrder,
                                keepTemp) {
  if (identical(genoFormat, "plink")) {
    bed <- paste0(genoFile, ".bed")
    bim <- paste0(genoFile, ".bim")
    fam <- paste0(genoFile, ".fam")
    for (f in c(bed, bim, fam)) {
      if (!file.exists(f)) stop("Missing PLINK file: ", f)
    }
    return(list(prefix = genoFile, cleanup = character(0)))
  }

  if (!file.exists(genoFile)) {
    stop("Genotype input file not found: ", genoFile)
  }

  temp_prefix <- file.path(
    tempdir(),
    paste0(sanitize_local_name(basename(genoFile)), "_check_tmp_", genoFormat, "_plink")
  )
  cleanup <- paste0(temp_prefix, c(".bed", ".bim", ".fam", ".log", ".nosex"))

  if (identical(genoFormat, "vcf")) {
    vcfField <- toupper(trimws(as.character(vcfField)[1]))
    if (!vcfField %in% c("DS", "GT")) {
      stop("`vcfField` must be either 'DS' or 'GT'.")
    }

    plink2_exec <- find_executable("plink2")
    if (nzchar(plink2_exec)) {
      args <- c("--vcf", genoFile)
      if (identical(vcfField, "DS")) {
        args <- c(args, "dosage=DS")
      }
      args <- c(args, "--double-id", "--keep-allele-order", "--make-bed", "--out", temp_prefix)
      run_command_checked(plink2_exec, args)
    } else {
      plink_exec <- find_executable("plink")
      if (!nzchar(plink_exec)) {
        stop("Neither 'plink2' nor 'plink' was found in PATH; cannot convert VCF input.")
      }
      args <- c("--vcf", genoFile)
      if (identical(vcfField, "DS")) {
        args <- c(args, "dosage=DS")
      }
      args <- c(args, "--double-id", "--keep-allele-order", "--make-bed", "--out", temp_prefix)
      run_command_checked(plink_exec, args)
    }
  } else if (identical(genoFormat, "bgen")) {
    plink2_exec <- find_executable("plink2")
    if (!nzchar(plink2_exec)) {
      stop("BGEN input requires plink2, but executable 'plink2' was not found in PATH.")
    }
    if (is.null(alleleOrder) || !nzchar(as.character(alleleOrder)[1]) || toupper(as.character(alleleOrder)[1]) == "NULL") {
      alleleOrder <- "ref-last"
    }
    alleleOrder <- trimws(as.character(alleleOrder)[1])
    if (!alleleOrder %in% c("ref-first", "ref-last", "ref-unknown")) {
      stop("`alleleOrder` must be 'ref-first', 'ref-last', or 'ref-unknown' for BGEN input.")
    }

    sample_candidates <- c(
      sub("\\.bgen$", ".sample", genoFile, ignore.case = TRUE),
      sub("\\.bgen$", ".samples", genoFile, ignore.case = TRUE),
      paste0(genoFile, ".sample"),
      paste0(genoFile, ".samples")
    )
    sample_candidates <- unique(sample_candidates[file.exists(sample_candidates)])
    sample_file <- character(0)
    if (length(sample_candidates) > 0L) {
      sample_file <- sample_candidates[1]
      cat("Detected BGEN sample file:", sample_file, "\n")
    } else {
      cat("No BGEN sample file detected; plink2 will rely on sample IDs embedded in the .bgen file.\n")
    }

    args <- c("--bgen", genoFile, alleleOrder)
    if (length(sample_file) == 1L && nzchar(sample_file)) {
      args <- c(args, "--sample", sample_file)
    }
    args <- c(args, "--make-bed", "--out", temp_prefix)
    run_command_checked(plink2_exec, args)
  } else {
    stop("Unsupported inferred genotype format: ", genoFormat)
  }

  for (f in paste0(temp_prefix, c(".bed", ".bim", ".fam"))) {
    if (!file.exists(f)) {
      stop("Conversion to temporary PLINK files failed; missing output: ", f)
    }
  }

  if (keepTemp) {
    cat("Keeping temporary converted PLINK files with prefix:", temp_prefix, "\n")
  }

  list(prefix = temp_prefix, cleanup = cleanup)
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
if (is.null(opt$alleleOrder) || !nzchar(opt$alleleOrder) || toupper(opt$alleleOrder) == "NULL") {
  opt$alleleOrder <- NULL
}

# 1) Strict set check between abundance/covariates
if (!setequal(abd_ids, cov_ids)) stop("IID set mismatch: abdFile vs covFile.")

# 2) Always use genotype sample order as the reference order
changed <- FALSE
geno_input <- prepare_plink_input(
  genoFile = opt$genoFile,
  genoFormat = geno_format,
  vcfField = opt$vcfField,
  alleleOrder = opt$alleleOrder,
  keepTemp = opt$keepTemp
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
