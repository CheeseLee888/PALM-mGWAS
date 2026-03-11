#' Run PALM summary statistics for PLINK genotype inputs
#'
#' Reads PLINK bed/bim/fam files, applies a pre-fitted PALM null model, and
#' writes per-phenotype summary statistics (one file per phenotype).
#'
#' @param genoPrefix Path prefix to PLINK files (without extension); expects
#'   `<genoPrefix>.bed/.bim/.fam` to exist.
#' @param NULLmodelFile Path to `.rda` containing the fitted null model
#'   object named `modglmm`.
#' @param PALMOutputFile Output prefix for per-phenotype result files; each
#'   phenotype file will be written as `<PALMOutputFile>_<pheno>.txt`.
#' @param chrom Optional chromosome filter (numeric or string like `"chr1"`).
#'   Use `NULL` to keep all chromosomes.
#' @param correct Passed to `PALM::palm.get.summary()`; defaults to `"NULL"`.
#' @param useCluster Logical; if `TRUE`, uses FID from `.fam` as cluster
#'   information when available.
#'
#' @return Invisibly returns a character vector of written file paths.
#' @export
getSummary <- function(genoPrefix,
                       NULLmodelFile,
                       PALMOutputFile,
                       chrom = NULL,
                       correct = "NULL",
                       useCluster = TRUE) {
  if (!requireNamespace("PALM", quietly = TRUE)) {
    stop("Package 'PALM' is required but not installed.")
  }
  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("Package 'snpStats' is required but not installed.")
  }

  if (missing(genoPrefix) || !nzchar(genoPrefix)) {
    stop("'genoPrefix' must be provided (PLINK prefix without extension).")
  }
  if (missing(NULLmodelFile) || !nzchar(NULLmodelFile)) {
    stop("'NULLmodelFile' must be provided.")
  }
  if (!file.exists(NULLmodelFile)) {
    stop("NULL model file not found: ", NULLmodelFile)
  }
  if (missing(PALMOutputFile) || !nzchar(PALMOutputFile)) {
    stop("'PALMOutputFile' must be provided.")
  }

  env <- new.env()
  load(NULLmodelFile, envir = env) # load modglmm
  if (!exists("modglmm", envir = env)) {
    stop("Object 'modglmm' not found in ", NULLmodelFile)
  }
  modglmm <- env$modglmm
  message("Loaded NULL model from ", NULLmodelFile)

  # normalize chrom input: treat "" or "NULL" as NULL
  if (!is.null(chrom) && is.character(chrom)) {
    if (!nzchar(chrom) || toupper(chrom) == "NULL") chrom <- NULL
  }
  if (is.null(chrom)) {
    message("Chromosome filter disabled: using all SNPs.")
  } else {
    message("Chromosome filter enabled: requested chromosome ", chrom)
  }
  if (is.null(correct)) {
    message("Compositional correction disabled: correct=NULL")
  } else {
    message("Compositional correction enabled: correct=", correct)
  }
  message("Cluster option requested: useCluster=", useCluster)

  # read genotype data and make it a data.frame
  bed <- paste0(genoPrefix, ".bed")
  bim <- paste0(genoPrefix, ".bim")
  fam <- paste0(genoPrefix, ".fam")
  for (f in c(bed, bim, fam)) if (!file.exists(f)) stop("Missing PLINK file: ", f)

  # read fam to get cluster info (make sure IDs align with abdFile's)
  cluster <- NULL
  if (useCluster) {
    message("Reading PLINK .fam file for cluster info.")
    fam_data <- utils::read.table(fam, stringsAsFactors = FALSE)
    colnames(fam_data) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
    cluster <- fam_data$FID
    names(cluster) <- fam_data$IID
    ## Handle FID = 0 case
    if (all(cluster == 0)) {
      message("All FID values are 0. No valid cluster information detected. Setting useCluster = FALSE.")
      useCluster <- FALSE
      cluster <- NULL
    }
  }
  # ----------------------------

  plink <- snpStats::read.plink(bed, bim, fam)

  G <- plink$genotypes

  iid <- rownames(G)
  if (is.null(iid)) stop("No rownames (IID) found in genotype matrix from read.plink().")

  # Convert to numeric 0/1/2/NA matrix, then to data.frame
  geno <- as(G, "numeric") # returns matrix with 0/1/2 and NA
  rownames(geno) <- iid
  colnames(geno) <- colnames(G)
  message("Loaded genotype matrix: ", nrow(geno), " samples x ", ncol(geno), " SNPs before chromosome filtering.")

  # Subset for quick testing (every 10th SNP)
  # geno <- geno[, seq(1, ncol(geno), by = 10), drop = FALSE]

  # --------------------------
  # Subset by chromosome if specified
  # --------------------------
  if (!is.null(chrom) && toupper(chrom) != "NULL") {
    message("Subsetting genotype data for chromosome: ", chrom)
    chrom <- sub("^chr", "", chrom, ignore.case = TRUE)

    # plink$map has chromosome info from .bim
    chr_map <- as.character(plink$map$chromosome)
    chr_map <- sub("^chr", "", chr_map, ignore.case = TRUE)
    keep <- which(chr_map == chrom)

    if (length(keep) == 0L) stop("No SNPs found for --chrom=", chrom)
    geno <- geno[, keep, drop = FALSE]
    message("Genotype matrix after chromosome filtering: ", nrow(geno), " samples x ", ncol(geno), " SNPs.")
  }


  # # Save dosage matrix
  # write.table(
  #   geno,
  #   file = file.path(PALMOutputFile, "geno_allele2.txt"),
  #   sep = "\t",
  #   quote = FALSE,
  #   col.names = NA
  # )
  # cat("Wrote full allele.2 dosage matrix to: ",
  #     file.path(PALMOutputFile, "geno_allele2_012_full.txt"), "\n"
  # )

  # ---------- run palm.get.summary ----------
  if (!useCluster) {
    message("No cluster provided; running palm.get.summary without cluster.")
    res <- PALM::palm.get.summary(
      null.obj = modglmm,
      covariate.interest = geno,
      correct = correct
    )
  } else {
    message("Cluster provided; running palm.get.summary with cluster (FID in PLINK file).")
    # print some cluster IDs
    uclust <- unique(cluster)
    message("Unique cluster IDs (FID): ", length(uclust))
    max_print <- 10
    if (length(uclust) <= max_print) {
      print(uclust)
    } else {
      print(uclust[1:max_print])
      message("... (", length(uclust) - max_print, " more clusters omitted)")
    }

    tab <- table(cluster)
    message(
      "Cluster size (min/median/max): ",
      min(tab), "/", as.numeric(stats::median(tab)), "/", max(tab)
    )

    res <- PALM::palm.get.summary(
      null.obj = modglmm,
      covariate.interest = geno,
      correct = correct,
      cluster = cluster
    )
  }

  # ----------------------------
  # Write results to files (one per phenotype)
  # ----------------------------
  message("Writing results in split files by phenotype.")

  stopifnot(length(res) >= 1)
  study_names <- names(res)
  if (is.null(study_names) || any(study_names == "")) study_names <- paste0("Study", seq_along(res))
  names(res) <- study_names

  res_df_list <- lapply(study_names, function(d) {
    est_df <- as.data.frame(res[[d]]$est, check.names = FALSE)
    se_df <- as.data.frame(res[[d]]$stderr, check.names = FALSE)

    colnames(est_df) <- paste0(d, ".est.", colnames(est_df))
    colnames(se_df) <- paste0(d, ".stderr.", colnames(se_df))

    n_df <- data.frame(tmp = res[[d]]$n)
    colnames(n_df) <- paste0(d, ".n")

    cbind(est_df, se_df, n_df)
  })

  res <- res_df_list[[1]]
  if (length(res_df_list) > 1) {
    for (i in 2:length(res_df_list)) {
      # feature rows align (usually rownames are feature IDs); if not aligned, merge by rownames
      res <- cbind(res, res_df_list[[i]])
    }
  }

  # inherit rownames（feature IDs）
  rownames(res) <- rownames(res_df_list[[1]])


  # Automatically infer the study prefix (usually "Study")
  prefix <- sub("\\.est\\..*$", "", grep("\\.est\\.", colnames(res), value = TRUE)[1])
  if (is.na(prefix) || prefix == "") prefix <- "Study"

  est_pat <- paste0("^", prefix, "\\.est\\.")
  stderr_pat <- paste0("^", prefix, "\\.stderr\\.")

  est_cols <- grep(est_pat, colnames(res), value = TRUE)
  stderr_cols <- grep(stderr_pat, colnames(res), value = TRUE)

  if (length(est_cols) == 0 || length(stderr_cols) == 0) {
    stop(
      "Cannot find est/stderr columns in res. Example colnames(res): ",
      paste(head(colnames(res), 5), collapse = ", ")
    )
  }

  snp_est <- sub(est_pat, "", est_cols)
  snp_stderr <- sub(stderr_pat, "", stderr_cols)
  common_snp <- intersect(snp_est, snp_stderr)
  if (length(common_snp) == 0) stop("No matched SNPs between est and stderr columns.")

  # keep SNP order as in est columns
  common_snp <- snp_est[snp_est %in% common_snp]
  est_map <- setNames(est_cols, snp_est)
  stderr_map <- setNames(stderr_cols, snp_stderr)

  if (is.null(rownames(res)) || any(rownames(res) == "")) {
    stop("res has no rownames (phenotype names). Please ensure rownames(res)=pheno names.")
  }

  # PALMOutputFile can be a "directory" or "prefix"; here we treat it as a directory for clarity
  out_dir <- dirname(PALMOutputFile)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  out_paths <- character(0)
  for (pheno in rownames(res)) {
    out <- data.frame(
      SNP = common_snp,
      est = as.numeric(res[pheno, est_map[common_snp], drop = TRUE]),
      stderr = as.numeric(res[pheno, stderr_map[common_snp], drop = TRUE]),
      check.names = FALSE
    )

    ## compute p-value
    # out$stat  <- out$est / out$stderr
    # out$pval <- 2 * pnorm(-abs(out$stat))
    out$pval <- 1 - pchisq((out$est / out$stderr)^2, df = 1)


    ## ----------------------------
    ## Post-process SNP format and derive CHR/POS
    ## ----------------------------

    # 1) Replace "." with ":" in SNP (e.g., chr1.1119172.G.A -> chr1:1119172:G:A)
    out$SNP <- gsub("\\.", ":", out$SNP)

    # 2) Parse CHR and POS from SNP (expects chr:pos:...)
    parts <- data.table::tstrsplit(out$SNP, ":", fixed = TRUE)
    if (length(parts) < 2) {
      stop("SNP format invalid after conversion; expected at least chr:pos:... (e.g., chr1:12345:A:G)")
    }

    chr_raw <- parts[[1]]
    pos_raw <- parts[[2]]

    # 3) Clean chromosome labels and convert to integer codes
    chr_clean <- sub("^chr", "", chr_raw, ignore.case = TRUE)
    chr_clean[chr_clean %in% c("X", "x")] <- "23"
    chr_clean[chr_clean %in% c("Y", "y")] <- "24"
    chr_clean[chr_clean %in% c("MT", "Mt", "mt", "M", "m")] <- "25"

    out$CHR <- suppressWarnings(as.integer(chr_clean))
    out$POS <- suppressWarnings(as.integer(pos_raw))

    # (Optional) Reorder columns
    out <- out[, c("SNP", "CHR", "POS", "est", "stderr", "pval")]


    out_file <- paste0(PALMOutputFile, "_", pheno, ".txt")

    write.table(out,
      file = out_file, sep = "\t",
      quote = FALSE, row.names = FALSE, col.names = TRUE
    )
    message("Wrote per-pheno file: ", out_file)
    out_paths <- c(out_paths, out_file)
  }

  message("Done. PALM summary results saved with prefix: ", PALMOutputFile)
  invisible(out_paths)
}
