#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(snpStats)
})

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    stop("Unable to determine script path from commandArgs().")
  }
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
}

script_path <- get_script_path()
simu_root <- normalizePath(dirname(script_path), mustWork = TRUE)

seed <- 20260409L
set.seed(seed)

provided_dir <- file.path(simu_root, "provided")
variant_file <- file.path(provided_dir, "snp_reference.tsv")
feature_file <- file.path(provided_dir, "feature_names.txt")
covariate_file <- file.path(provided_dir, "covariate_names.txt")
sample_file <- file.path(provided_dir, "sample_config.tsv")

input_dir <- file.path(simu_root, "input", "study1")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading simulation metadata from provided/.")
variant_dt <- fread(variant_file)
required_variant_cols <- c("CHR", "SNP", "POS", "A1", "A2", "N", "AF")
missing_variant_cols <- setdiff(required_variant_cols, names(variant_dt))
if (length(missing_variant_cols) > 0L) {
  stop("Missing required columns in snp_reference.tsv: ", paste(missing_variant_cols, collapse = ", "))
}

feature_names <- readLines(feature_file, warn = FALSE)
feature_names <- feature_names[nzchar(feature_names)]
covariate_names <- readLines(covariate_file, warn = FALSE)
covariate_names <- covariate_names[nzchar(covariate_names)]
sample_cfg <- fread(sample_file)
if (!all(c("key", "value") %in% names(sample_cfg))) {
  stop("sample_config.tsv must contain columns 'key' and 'value'.")
}
sample_cfg <- setNames(sample_cfg$value, sample_cfg$key)
sample_n <- as.integer(sample_cfg[["sample_n"]])
sample_prefix <- sample_cfg[["sample_prefix"]]
if (is.na(sample_n) || sample_n <= 0L) {
  stop("sample_n in sample_config.tsv must be a positive integer.")
}
if (is.na(sample_prefix) || !nzchar(sample_prefix)) {
  stop("sample_prefix in sample_config.tsv must be non-empty.")
}
sample_ids <- sprintf(
  paste0(sample_prefix, "_%0", nchar(as.character(sample_n)), "d"),
  seq_len(sample_n)
)

if (!all(c("g_Acinetobacter", "f_Moraxellaceae") %in% feature_names)) {
  stop("feature_names.txt must contain g_Acinetobacter and f_Moraxellaceae.")
}

message("Generating reproducible covariates for ", sample_n, " samples.")
cov_dt <- data.table(IID = sample_ids)
for (j in seq_along(covariate_names)) {
  cov_name <- covariate_names[j]
  if (cov_name == "AGE") {
    cov_dt[, (cov_name) := as.numeric(scale(rnorm(.N, mean = 45, sd = 12)))]
  } else if (cov_name == "SEX") {
    cov_dt[, (cov_name) := as.integer(rbinom(.N, size = 1, prob = 0.48))]
  } else if (grepl("^PC[0-9]+$", cov_name)) {
    cov_dt[, (cov_name) := rnorm(.N, mean = 0, sd = 0.012 / j)]
  } else {
    cov_dt[, (cov_name) := as.numeric(scale(rnorm(.N)))]
  }
}

message("Generating genotype matrix from simulation AF values.")
af <- pmin(pmax(variant_dt$AF, 0.02), 0.98)
geno_vec <- rbinom(sample_n * nrow(variant_dt), size = 2L, prob = rep(af, each = sample_n))
geno <- matrix(geno_vec, nrow = sample_n, ncol = nrow(variant_dt))
rm(geno_vec)
gc(verbose = FALSE)
rownames(geno) <- sample_ids
colnames(geno) <- variant_dt$SNP

causal_snps <- c(
  "chr4:1682869:G:C",
  "chr3:68821275:C:T",
  "chr7:81941667:G:C"
)
missing_causal <- setdiff(causal_snps, colnames(geno))
if (length(missing_causal) > 0L) {
  stop("Missing causal SNP(s): ", paste(missing_causal, collapse = ", "))
}

message("Generating abundance matrix with planted signal in g_Acinetobacter.")
feature_n <- length(feature_names)
size_factor <- rgamma(sample_n, shape = 30, rate = 30)
abd_mat <- matrix(
  0L,
  nrow = sample_n,
  ncol = feature_n,
  dimnames = list(sample_ids, feature_names)
)

feature_prev <- seq(0.30, 0.92, length.out = feature_n)
feature_mean <- exp(seq(log(20), log(260), length.out = feature_n))
age_z <- as.numeric(scale(cov_dt$AGE))
sex_centered <- cov_dt$SEX - mean(cov_dt$SEX)

for (j in seq_len(feature_n)) {
  present <- rbinom(sample_n, size = 1L, prob = feature_prev[j])
  lambda <- feature_mean[j] * size_factor * exp(0.18 * age_z + 0.08 * sex_centered)
  counts <- rnbinom(sample_n, mu = lambda, size = 12 + (j %% 9))
  abd_mat[, j] <- as.integer(present * counts)
}

genetic_score <- as.numeric(scale(
  geno[, causal_snps, drop = FALSE] %*% c(1.8, 1.0, 0.8)
))

target_mu <- exp(6.2 + 1.35 * genetic_score + 0.20 * age_z + 0.15 * sex_centered)
target_counts <- rnbinom(sample_n, mu = target_mu, size = 18)
abd_mat[, "g_Acinetobacter"] <- as.integer(target_counts)

if ("f_Moraxellaceae" %in% colnames(abd_mat)) {
  abd_mat[, "f_Moraxellaceae"] <- abd_mat[, "g_Acinetobacter"] +
    rnbinom(sample_n, mu = 30, size = 20)
}

abd_dt <- data.table(IID = sample_ids)
abd_dt <- cbind(abd_dt, as.data.table(abd_mat))

message("Writing abundance and covariate tables.")
fwrite(abd_dt, file.path(input_dir, "abd.txt"), sep = "\t", quote = FALSE)
fwrite(cov_dt, file.path(input_dir, "cov.txt"), sep = "\t", quote = FALSE)

message("Writing PLINK files.")
snp_matrix <- as(geno, "SnpMatrix")
write.plink(
  file.base = file.path(input_dir, "geno"),
  snps = snp_matrix,
  pedigree = rep(0, sample_n),
  id = sample_ids,
  father = rep(0, sample_n),
  mother = rep(0, sample_n),
  sex = rep(0, sample_n),
  phenotype = rep(-9, sample_n),
  chromosome = variant_dt$CHR,
  genetic.distance = rep(0, nrow(variant_dt)),
  position = variant_dt$POS,
  allele.1 = variant_dt$A1,
  allele.2 = variant_dt$A2
)

signal_summary <- data.table(
  seed = seed,
  sample_n = sample_n,
  snp_n = nrow(variant_dt),
  causal_snp = causal_snps
)
fwrite(signal_summary, file.path(input_dir, "simulation_manifest.txt"), sep = "\t", quote = FALSE)

message("Done.")
message("Seed: ", seed)
message("Samples: ", sample_n)
message("SNPs: ", nrow(variant_dt))
message("Causal SNPs: ", paste(causal_snps, collapse = ", "))
