#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

# options(stringsAsFactors=F, scipen = 999)
options(stringsAsFactors = F)
library(SAIGEQTL)


BLASctl_installed <- require(RhpcBLASctl)
library(optparse)
library(data.table)
library(methods)
print(sessionInfo())

option_list <- list(
  make_option("--inFile",
    type = "character", default = "",
    help = "Path to geno file. Support many formats including bed/bim/fam, bgen, gds and txt. If using bed/bim/fam, only go with prefix; Otherwise full file name."
  ),
  make_option("--NULLmodelFile",
    type = "character", default = "",
    help = "Path to the input file containing the glmm model, which is output from previous step. Will be used by load()"
  ),
  make_option("--PALMOutputFile",
    type = "character", default = "",
    help = "Path to the output file containing assoc test results"
  ),
  make_option("--chrom",
    type = "character", default = "",
    help = paste(
      "Chromosome to analyze (e.g. 1, 2, X).",
      "Currently only supported for PLINK bed/bim/fam input.",
      "If empty, use all variants in inFile."
    )
  ),
  make_option("--phenoCol",
    type = "character", default = "",
    help = "Column name for the phenotype to be tested in the NULLmodelFile."
  ),
  make_option("--minMAF",
    type = "numeric", default = 0,
    help = "Minimum minor allele frequency of markers to test. By default 0."
  )
)


parser <- OptionParser(usage = "%prog [options]", option_list = option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)


# --------------------------
# Handle --chrom / Subset PLINK
# --------------------------
inFile_to_use <- opt$inFile  # Default: directly use the user-provided inFile

if (!is.null(opt$chrom) && toupper(opt$chrom) != "NULL") {
  message("Chromosome specified: ", opt$chrom)

  # Currently, --chrom is only implemented for PLINK bed/bim/fam
  bed_file <- paste0(opt$inFile, ".bed")
  bim_file <- paste0(opt$inFile, ".bim")
  fam_file <- paste0(opt$inFile, ".fam")

  if (!file.exists(bed_file) || !file.exists(bim_file) || !file.exists(fam_file)) {
    stop(
      "--chrom is currently only supported when --inFile is a PLINK prefix. ",
      "Expected files: ", bed_file, ", ", bim_file, ", ", fam_file
    )
  }

  # Temporary prefix: total prefix + .chrX_tmp
  tmp_prefix <- sprintf("%s.chr%s_tmp", opt$inFile, opt$chrom)

  plink_args <- c(
    "--bfile", opt$inFile,
    "--chr",  opt$chrom,
    "--make-bed",
    "--out", tmp_prefix
  )

  message("Subsetting PLINK by chromosome using command:")
  message("  plink ", paste(plink_args, collapse = " "))

  status <- system2("plink", args = plink_args)

  if (!identical(status, 0L)) {
    stop("PLINK failed with non-zero exit status: ", status)
  }

  # Use this sub-prefix to run SPAGMMATtest
  inFile_to_use <- tmp_prefix

}else{
  message("No chromosome specified; using all variants in inFile.")
}

# --------------------------
# Call SPAGMMATtest
# --------------------------
SPAGMMATtest(
  inFile = inFile_to_use,
  phenoCol = opt$phenoCol,
  NULLmodelFile = opt$NULLmodelFile,
  PALMOutputFile = opt$PALMOutputFile,
  minMAF = opt$minMAF
)

# --------------------------
# Cleanup temporary files
# --------------------------
if (!is.null(opt$chrom) && toupper(opt$chrom) != "NULL") {
  tmp_files <- Sys.glob(paste0(tmp_prefix, ".*"))
  if (length(tmp_files) > 0L) {
    message(
      "Cleaning up temporary PLINK files: ",
      paste(tmp_files, collapse = ", ")
    )
    try(file.remove(tmp_files), silent = TRUE)
  }
}
