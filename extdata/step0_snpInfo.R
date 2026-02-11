#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(snpStats)
})

# -----------------------------
# CLI
# -----------------------------
option_list <- list(
  make_option(c("--genoFile"), type="character", default=NULL,
              help="PLINK prefix (without .bed/.bim/.fam)"),
  make_option(c("--outputFile"), type="character", default=NULL,
              help="Output file path. [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$genoFile)) stop("--genoFile is required")

# -----------------------------
# Helpers
# -----------------------------
# Parse SNP string like "chr1:1119172:G:A"
parse_snp <- function(snp_ids) {
  parts <- strsplit(snp_ids, ":", fixed = TRUE)
  lens  <- vapply(parts, length, integer(1))
  if (any(lens < 4)) {
    bad <- snp_ids[lens < 4][1]
    stop("SNP id format error. Expect 'chr<CHR>:<POS>:<A1>:<A2>'. Bad SNP: ", bad)
  }

  chr_raw <- vapply(parts, `[[`, character(1), 1)
  pos_raw <- vapply(parts, `[[`, character(1), 2)
  a1      <- vapply(parts, `[[`, character(1), 4)
  a2      <- vapply(parts, `[[`, character(1), 3)

  chr <- gsub("^(chr|CHR)", "", chr_raw)
  pos <- suppressWarnings(as.integer(pos_raw))
  if (any(is.na(pos))) {
    bad <- snp_ids[is.na(pos)][1]
    stop("POS parse failed for SNP: ", bad)
  }

  data.frame(
    CHR = chr,
    SNP = snp_ids,
    POS = pos,
    A1  = a1,
    A2  = a2,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

# -----------------------------
# Read PLINK
# -----------------------------
bed <- paste0(opt$genoFile, ".bed")
bim <- paste0(opt$genoFile, ".bim")
fam <- paste0(opt$genoFile, ".fam")

message("Reading PLINK: ", opt$genoFile)
plink <- snpStats::read.plink(bed, bim, fam)

# Convert genotypes to numeric 0/1/2/NA (counts of allele.2)
G <- as(plink$genotypes, "numeric")
snps <- colnames(G)
if (is.null(snps)) stop("No SNP names found in PLINK genotypes.")

# -----------------------------
# Compute N and AF (A2)
# -----------------------------
N   <- colSums(!is.na(G))
MAC <- colSums(G, na.rm = TRUE)   # sum of A2 dosages
AF  <- rep(NA_real_, length(snps))
ok  <- N > 0
AF[ok] <- (MAC[ok] / N[ok]) / 2

# -----------------------------
# Build output
# -----------------------------
info <- parse_snp(snps)
info$N  <- as.integer(N)
info$AF <- as.numeric(AF)

out <- info[, c("CHR","SNP","POS","A1","A2","N","AF")]

message("Writing: ", opt$outputFile)
data.table::fwrite(out, file = opt$outputFile, sep = "\t", quote = FALSE, na = "NA")
