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

study_cfg <- list(
  study1 = list(
    sample_n = sample_n,
    sample_prefix = paste0(sample_prefix, "1"),
    age_mean = 45,
    age_sd = 12,
    sex_prob = 0.48,
    pc_scale = 1.00,
    genetic_beta = c(0.78, 0.44, 0.36),
    signal_genetic = 0.58,
    signal_intercept = 6.2,
    signal_age = 0.20,
    signal_sex = 0.15,
    drop_features = character(0),
    snp_drop_every = 0L,
    snp_drop_offset = 0L
  ),
  study2 = list(
    sample_n = 430L,
    sample_prefix = paste0(sample_prefix, "2"),
    age_mean = 49,
    age_sd = 11,
    sex_prob = 0.44,
    pc_scale = 1.15,
    genetic_beta = c(0.70, 0.40, 0.32),
    signal_genetic = 0.54,
    signal_intercept = 6.0,
    signal_age = 0.16,
    signal_sex = 0.12,
    drop_features = c("g_Negativicoccus", "g_Atopobium", "f_Tannerellaceae"),
    snp_drop_every = 29L,
    snp_drop_offset = 7L
  ),
  study3 = list(
    sample_n = 610L,
    sample_prefix = paste0(sample_prefix, "3"),
    age_mean = 42,
    age_sd = 13,
    sex_prob = 0.51,
    pc_scale = 0.90,
    genetic_beta = c(0.86, 0.48, 0.40),
    signal_genetic = 0.60,
    signal_intercept = 6.3,
    signal_age = 0.22,
    signal_sex = 0.18,
    drop_features = c("g_Terrisporobacter", "f_Erysipelotrichaceae", "f_Pasteurellaceae"),
    snp_drop_every = 31L,
    snp_drop_offset = 11L
  )
)

causal_snps <- c(
  "chr4:1682869:G:C",
  "chr3:68821275:C:T",
  "chr7:81941667:G:C"
)
missing_causal <- setdiff(causal_snps, variant_dt$SNP)
if (length(missing_causal) > 0L) {
  stop("Missing causal SNP(s): ", paste(missing_causal, collapse = ", "))
}

keep_variant_rows <- function(cfg, variant_dt, causal_snps) {
  keep <- rep(TRUE, nrow(variant_dt))
  if (cfg$snp_drop_every > 0L) {
    idx <- seq.int(from = cfg$snp_drop_offset + 1L, to = nrow(variant_dt), by = cfg$snp_drop_every)
    keep[idx] <- FALSE
  }
  keep[variant_dt$SNP %in% causal_snps] <- TRUE
  keep
}

for (study_name in names(study_cfg)) {
  cfg <- study_cfg[[study_name]]
  sample_ids <- sprintf(
    paste0(cfg$sample_prefix, "_%0", nchar(as.character(cfg$sample_n)), "d"),
    seq_len(cfg$sample_n)
  )
  input_dir <- file.path(simu_root, "input", study_name)
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

  message("Generating reproducible covariates for ", study_name, " with ", cfg$sample_n, " samples.")
  cov_dt <- data.table(IID = sample_ids)
  for (j in seq_along(covariate_names)) {
    cov_name <- covariate_names[j]
    if (cov_name == "AGE") {
      cov_dt[, (cov_name) := as.numeric(scale(rnorm(.N, mean = cfg$age_mean, sd = cfg$age_sd)))]
    } else if (cov_name == "SEX") {
      cov_dt[, (cov_name) := as.integer(rbinom(.N, size = 1, prob = cfg$sex_prob))]
    } else if (grepl("^PC[0-9]+$", cov_name)) {
      cov_dt[, (cov_name) := rnorm(.N, mean = 0, sd = cfg$pc_scale * 0.012 / j)]
    } else {
      cov_dt[, (cov_name) := as.numeric(scale(rnorm(.N)))]
    }
  }

  variant_keep <- keep_variant_rows(cfg, variant_dt, causal_snps)
  variant_sub <- copy(variant_dt[variant_keep])

  message(
    "Generating genotype matrix for ", study_name,
    " using ", nrow(variant_sub), " SNP(s)."
  )
  af <- pmin(pmax(variant_sub$AF, 0.02), 0.98)
  geno_vec <- rbinom(cfg$sample_n * nrow(variant_sub), size = 2L, prob = rep(af, each = cfg$sample_n))
  geno <- matrix(geno_vec, nrow = cfg$sample_n, ncol = nrow(variant_sub))
  rm(geno_vec)
  gc(verbose = FALSE)
  rownames(geno) <- sample_ids
  colnames(geno) <- variant_sub$SNP

  feature_keep <- setdiff(feature_names, cfg$drop_features)
  feature_n <- length(feature_keep)
  size_factor <- rgamma(cfg$sample_n, shape = 30, rate = 30)
  abd_mat <- matrix(
    0L,
    nrow = cfg$sample_n,
    ncol = feature_n,
    dimnames = list(sample_ids, feature_keep)
  )

  message(
    "Generating abundance matrix for ", study_name,
    " using ", feature_n, " feature(s)."
  )
  feature_prev <- seq(0.30, 0.92, length.out = feature_n)
  feature_mean <- exp(seq(log(20), log(260), length.out = feature_n))
  age_z <- as.numeric(scale(cov_dt$AGE))
  sex_centered <- cov_dt$SEX - mean(cov_dt$SEX)

  for (j in seq_len(feature_n)) {
    present <- rbinom(cfg$sample_n, size = 1L, prob = feature_prev[j])
    lambda <- feature_mean[j] * size_factor * exp(0.18 * age_z + 0.08 * sex_centered)
    counts <- rnbinom(cfg$sample_n, mu = lambda, size = 12 + (j %% 9))
    abd_mat[, j] <- as.integer(present * counts)
  }

  genetic_score <- as.numeric(scale(
    geno[, causal_snps, drop = FALSE] %*% cfg$genetic_beta
  ))

  target_mu <- exp(
    cfg$signal_intercept +
      cfg$signal_genetic * genetic_score +
      cfg$signal_age * age_z +
      cfg$signal_sex * sex_centered
  )
  target_counts <- rnbinom(cfg$sample_n, mu = target_mu, size = 18)
  if ("g_Acinetobacter" %in% colnames(abd_mat)) {
    abd_mat[, "g_Acinetobacter"] <- as.integer(target_counts)
  }

  if ("f_Moraxellaceae" %in% colnames(abd_mat)) {
    abd_mat[, "f_Moraxellaceae"] <- target_counts + rnbinom(cfg$sample_n, mu = 30, size = 20)
  }

  abd_dt <- data.table(IID = sample_ids)
  abd_dt <- cbind(abd_dt, as.data.table(abd_mat))

  message("Writing abundance and covariate tables for ", study_name, ".")
  fwrite(abd_dt, file.path(input_dir, "abd.txt"), sep = "\t", quote = FALSE)
  fwrite(cov_dt, file.path(input_dir, "cov.txt"), sep = "\t", quote = FALSE)

  message("Writing PLINK files for ", study_name, ".")
  snp_matrix <- as(geno, "SnpMatrix")
  write.plink(
    file.base = file.path(input_dir, "geno"),
    snps = snp_matrix,
    pedigree = rep(0, cfg$sample_n),
    id = sample_ids,
    father = rep(0, cfg$sample_n),
    mother = rep(0, cfg$sample_n),
    sex = rep(0, cfg$sample_n),
    phenotype = rep(-9, cfg$sample_n),
    chromosome = variant_sub$CHR,
    genetic.distance = rep(0, nrow(variant_sub)),
    position = variant_sub$POS,
    allele.1 = variant_sub$A1,
    allele.2 = variant_sub$A2
  )

  signal_summary <- data.table(
    seed = seed,
    study = study_name,
    sample_n = cfg$sample_n,
    snp_n = nrow(variant_sub),
    feature_n = feature_n,
    dropped_feature = paste(cfg$drop_features, collapse = ","),
    causal_snp = paste(causal_snps, collapse = ",")
  )
  fwrite(signal_summary, file.path(input_dir, "simulation_manifest.txt"), sep = "\t", quote = FALSE)
}

message("Done.")
message("Seed: ", seed)
message("Studies generated: ", paste(names(study_cfg), collapse = ", "))
message("Causal SNPs retained in all studies: ", paste(causal_snps, collapse = ", "))
