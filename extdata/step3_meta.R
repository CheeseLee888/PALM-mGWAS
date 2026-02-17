#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(PALM)
  library(data.table)
})

option_list <- list(
  make_option("--files", type="character", default="",
              help="Comma-separated .rda files from step2 (each contains object `res`)."),
  make_option("--outPrefix", type="character", default="",
              help="Output prefix. Will save <outPrefix>_objects.rda and <outPrefix>_<feature>.txt"),
  make_option("--p.adjust.method", type="character", default="fdr",
              help="p.adjust method used by PALM meta (default: fdr)"),
  make_option("--meta.method", type="character", default="EE",
              help="meta method for PALM (default: EE)")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (!nzchar(opt$files)) stop("Please provide --files (comma-separated .rda files).")
if (!nzchar(opt$outPrefix)) stop("Please provide --outPrefix.")

# Parse files vector
files <- trimws(unlist(strsplit(opt$files, ",")))
files <- files[nzchar(files)]
if (length(files) < 1) stop("No valid files parsed from --files.")
for (f in files) if (!file.exists(f)) stop("Missing file: ", f)

# -------- helper: load one step2 .rda and extract one-study summary --------
load_one_study <- function(rda_file) {
  e <- new.env(parent = emptyenv())
  load(rda_file, envir = e)

  # Expect object `res` (your step2 saved result) or fall back to `summary.stats`
  if (exists("res", envir = e, inherits = FALSE)) {
    rr <- get("res", envir = e)
  } else if (exists("summary.stats", envir = e, inherits = FALSE)) {
    rr <- get("summary.stats", envir = e)
  } else {
    stop("File ", rda_file, " does not contain `res` or `summary.stats`.")
  }

  # rr could be:
  # (A) list length 1: rr$Study$est, rr$Study$stderr, rr$Study$n
  # (B) already a single-study list: rr$est, rr$stderr, rr$n
  if (is.list(rr) && length(rr) == 1 && is.list(rr[[1]]) &&
      all(c("est","stderr","n") %in% names(rr[[1]]))) {
    one <- rr[[1]]
  } else if (is.list(rr) && all(c("est","stderr","n") %in% names(rr))) {
    one <- rr
  } else {
    stop("Unexpected structure in ", rda_file,
         ". Expect `res$Study$est/stderr/n` or `res$est/stderr/n`.")
  }

  # normalize to list(est=matrix, stderr=matrix, n=integer)
  stopifnot(is.matrix(one$est), is.matrix(one$stderr), length(one$n) == 1)
  one
}

# -------- build summary.stats list --------
summary.stats <- list()
study_names <- sub("\\.rda$", "", basename(files))  # default: filename w/o extension
# keep unique names
if (any(duplicated(study_names))) {
  study_names <- make.unique(study_names)
}

for (i in seq_along(files)) {
  one <- load_one_study(files[i])
  summary.stats[[study_names[i]]] <- one
}

cat("Studies loaded: ", paste(names(summary.stats), collapse=", "), "\n", sep="")

# -------- meta --------
meta.result <- palm.meta.summary(
  summary.stats = summary.stats,
  p.adjust.method = opt$p.adjust.method,
  meta.method = opt$meta.method
)

# save RDA
dir.create(dirname(opt$outPrefix), showWarnings = FALSE, recursive = TRUE)
save(summary.stats, meta.result, file = paste0(opt$outPrefix, "_objects.rda"))
cat("Saved: ", paste0(opt$outPrefix, "_objects.rda"), "\n", sep="")

# -------- write per-feature meta txt (no allpheno) --------
# Your PALM returns meta.result as list of SNP -> data.frame(feature rows)
if (!is.list(meta.result) || length(meta.result) == 0) {
  stop("meta.result is empty or not a list; cannot write txt.")
}

meta_long <- rbindlist(
  lapply(names(meta.result), function(snp) {
    df <- as.data.table(meta.result[[snp]])
    df[, SNP := snp]
    df
  }),
  use.names = TRUE, fill = TRUE
)

# chr/pos from SNP like chr1.123.A.G
meta_long[, SNP := gsub("\\.", ":", SNP)]
parts <- tstrsplit(meta_long$SNP, ":", fixed = TRUE)
meta_long[, CHR := suppressWarnings(as.integer(sub("^chr", "", parts[[1]], ignore.case = TRUE)))]
meta_long[, POS := suppressWarnings(as.integer(parts[[2]]))]

# coef -> est
if ("coef" %in% names(meta_long)) setnames(meta_long, "coef", "est")

# order columns
front <- c("SNP","CHR","POS","feature","est","stderr","pval","qval","pval.het","qval.het")
meta_long <- meta_long[, ..front]

# split output
for (ph in unique(meta_long$feature)) {
  out_ph <- meta_long[feature == ph]
  out_ph[, feature := NULL]  # optional: drop feature column
  out_file <- paste0(opt$outPrefix, "_", ph, ".txt")
  fwrite(out_ph, file = out_file, sep = "\t", quote = FALSE, na = "NA")
  cat("Wrote: ", out_file, "\n", sep = "")
}

cat("Done.\n")
