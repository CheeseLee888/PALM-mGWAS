#' Run PALM summary statistics for genotype inputs
#'
#' Reads genotype data from PLINK or VCF input, applies a pre-fitted PALM null
#' model, and writes per-phenotype summary statistics (one file per phenotype).
#' VCF input is read directly.
#'
#' @param genoFile Path to genotype input. If the path ends with
#'   `.vcf`, `.vcf.gz`, or `.vcf.bgz`, VCF input is assumed. Otherwise it is
#'   treated as a PLINK prefix without extension.
#' @param NULLmodelFile Path to `.rda` containing the fitted null model
#'   object named `modglmm`.
#' @param PALMOutputFile Output prefix for per-phenotype result files; each
#'   phenotype file will be written as `<PALMOutputFile>_<pheno>.txt`.
#' @param vcfField Optional VCF FORMAT field override. By default the reader
#'   auto-detects and prefers `"DS"` when present, otherwise falls back to
#'   `"GT"`. Supported explicit values are `"DS"` and `"GT"`.
#' @param chrom Optional chromosome filter (numeric or string like `"chr1"`).
#'   Use `NULL` to keep all chromosomes.
#' @param featureColList Optional feature IDs to keep from the Step1 null model.
#'   Accepts either a character vector or a single comma-separated string.
#'   Use `NULL` (default) to keep all modeled features.
#' @param minMAF Optional minimum minor allele frequency threshold used to
#'   filter SNPs before running `PALM::palm.get.summary()`. Defaults to `0.05`.
#'   Use `0` to disable MAF filtering.
#' @param minMAC Optional minimum minor allele count threshold used to filter
#'   SNPs before running `PALM::palm.get.summary()`. Defaults to `5`. Use `0`
#'   to disable MAC filtering.
#' @param outputSnpFile Optional output path for SNP sample count and allele
#'   frequency computed from the Step2 genotype matrix after NULL-model sample
#'   alignment and optional chromosome subsetting, but before minMAF/minMAC
#'   filtering. Use `NULL` (default) to skip writing this file.
#' @param correct Passed to `PALM::palm.get.summary()`; defaults to `"NULL"`.
#' @param useCluster Logical; if `TRUE`, uses FID from `.fam` as cluster
#'   information when available. Defaults to `FALSE`.
#'
#' @return Invisibly returns a character vector of written file paths.
#' @export
getSummary <- function(genoFile,
                       NULLmodelFile,
                       PALMOutputFile,
                       vcfField = NULL,
                       chrom = NULL,
                       featureColList = NULL,
                       minMAF = 0.05,
                       minMAC = 5,
                       outputSnpFile = NULL,
                       correct = "NULL",
                       useCluster = FALSE) {
  if (!requireNamespace("PALM", quietly = TRUE)) {
    stop("Package 'PALM' is required but not installed.")
  }
  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("Package 'snpStats' is required but not installed.")
  }

  if (missing(genoFile) || !nzchar(genoFile)) {
    stop("'genoFile' must be provided.")
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
  if (!is.numeric(minMAF) || length(minMAF) != 1L || is.na(minMAF) || minMAF < 0 || minMAF > 0.5) {
    stop("'minMAF' must be a single numeric value between 0 and 0.5.")
  }
  if (!is.numeric(minMAC) || length(minMAC) != 1L || is.na(minMAC) || minMAC < 0 || minMAC != as.integer(minMAC)) {
    stop("'minMAC' must be a single non-negative integer value.")
  }
  output_dir <- dirname(PALMOutputFile)
  if (!output_dir %in% c("", ".")) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  genoFormat <- infer_geno_format(genoFile)
  message("Genotype input format inferred from genoFile: ", genoFormat)

  env <- new.env()
  load(NULLmodelFile, envir = env) # load modglmm
  if (!exists("modglmm", envir = env)) {
    stop("Object 'modglmm' not found in ", NULLmodelFile)
  }
  modglmm <- env$modglmm
  message("Loaded NULL model from ", NULLmodelFile)

  normalize_col_list <- function(x, arg_name) {
    if (is.null(x)) {
      return(NULL)
    }
    if (length(x) == 1L) {
      x <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
    } else {
      x <- trimws(as.character(x))
    }
    x <- x[nzchar(x)]
    if (!length(x)) {
      stop("'", arg_name, "' did not contain any valid column names.")
    }
    x
  }

  subset_null_model_features <- function(null_obj, feature_ids) {
    if (is.null(feature_ids)) {
      return(null_obj)
    }
    feature_by_study <- lapply(names(null_obj), function(d) colnames(null_obj[[d]]$Y_I))
    names(feature_by_study) <- names(null_obj)
    feature_union <- sort(unique(unlist(feature_by_study, use.names = FALSE)))
    missing_features <- setdiff(feature_ids, feature_union)
    if (length(missing_features) > 0L) {
      stop(
        "Requested feature(s) not found in NULL model: ",
        paste(missing_features, collapse = ", ")
      )
    }

    null_subset <- null_obj
    for (d in names(null_subset)) {
      keep <- intersect(feature_ids, colnames(null_subset[[d]]$Y_I))
      if (!length(keep)) {
        stop(
          "Requested feature(s) are absent from study '", d,
          "' in the NULL model after prevalence filtering."
        )
      }
      null_subset[[d]]$Y_I <- null_subset[[d]]$Y_I[, keep, drop = FALSE]
      null_subset[[d]]$Y_R <- null_subset[[d]]$Y_R[, keep, drop = FALSE]
    }
    null_subset
  }

  extract_null_sample_ids <- function(null_obj) {
    if (!is.list(null_obj) || length(null_obj) < 1L) {
      return(NULL)
    }
    first_study <- null_obj[[1]]
    if (!is.list(first_study)) {
      return(NULL)
    }
    for (candidate in c("Z", "Y_I", "Y_R")) {
      mat <- first_study[[candidate]]
      ids <- rownames(mat)
      if (!is.null(ids) && length(ids) > 0L) {
        return(ids)
      }
    }
    NULL
  }

  null_sample_ids <- extract_null_sample_ids(modglmm)
  if (is.null(null_sample_ids)) {
    message("Could not infer sample IDs from NULL model; using genotype rows as-is.")
  } else {
    message("NULL model sample count: ", length(null_sample_ids))
  }

  featureColList <- normalize_col_list(featureColList, "featureColList")
  if (is.null(featureColList)) {
    message("Feature filter disabled: using all modeled features.")
  } else {
    message("Feature filter enabled: requested ", length(featureColList), " feature(s).")
    modglmm <- subset_null_model_features(modglmm, featureColList)
    message(
      "NULL model after feature filtering: ",
      paste(
        sprintf("%s=%d", names(modglmm), vapply(modglmm, function(x) ncol(x$Y_I), integer(1))),
        collapse = "; "
      ),
      " feature(s) per study."
    )
  }

  # normalize chrom input: treat "" or "NULL" as NULL
  if (!is.null(chrom) && is.character(chrom)) {
    if (!nzchar(chrom) || toupper(chrom) == "NULL") chrom <- NULL
  }
  if (is.null(chrom)) {
    message("Chromosome filter disabled: using all SNPs.")
  } else {
    message("Chromosome filter enabled: requested chromosome ", chrom)
  }
  if (minMAF > 0) {
    message("MAF filter enabled: minMAF=", minMAF)
  } else {
    message("MAF filter disabled: minMAF=0")
  }
  if (minMAC > 0) {
    message("MAC filter enabled: minMAC=", minMAC)
  } else {
    message("MAC filter disabled: minMAC=0")
  }
  if (is.null(correct)) {
    message("Compositional correction disabled: correct=NULL")
  } else {
    message("Compositional correction enabled: correct=", correct)
  }
  message("Cluster option requested: useCluster=", useCluster)

  if (!identical(genoFormat, "plink") && isTRUE(useCluster)) {
    stop("`useCluster=TRUE` is only supported for native PLINK input.")
  }
  cluster <- NULL
  chr_map <- NULL
  if (identical(genoFormat, "vcf")) {
    vcf_input <- read_vcf_genotypes(genoFile, vcfField = vcfField)
    geno <- vcf_input$geno
    chr_map <- as.character(vcf_input$map$chromosome)
  } else {
    geno_input <- prepare_plink_input(
      genoFile = genoFile,
      tempLabel = "summary"
    )
    if (length(geno_input$cleanup) > 0L) {
      on.exit(unlink(geno_input$cleanup, force = TRUE), add = TRUE)
    }
    genoPrefix <- geno_input$prefix

    bed <- paste0(genoPrefix, ".bed")
    bim <- paste0(genoPrefix, ".bim")
    fam <- paste0(genoPrefix, ".fam")
    for (f in c(bed, bim, fam)) if (!file.exists(f)) stop("Missing PLINK file: ", f)

    if (useCluster) {
      message("Reading PLINK .fam file for cluster info.")
      fam_data <- utils::read.table(fam, stringsAsFactors = FALSE)
      colnames(fam_data) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
      cluster <- fam_data$FID
      names(cluster) <- fam_data$IID
      if (all(cluster == 0)) {
        message("All FID values are 0. No valid cluster information detected. Setting useCluster = FALSE.")
        useCluster <- FALSE
        cluster <- NULL
      }
    }

    plink <- snpStats::read.plink(bed, bim, fam)
    G <- plink$genotypes
    iid <- rownames(G)
    if (is.null(iid)) stop("No rownames (IID) found in genotype matrix from read.plink().")

    geno <- as(G, "numeric")
    rownames(geno) <- iid
    colnames(geno) <- colnames(G)
    chr_map <- as.character(plink$map$chromosome)
    message("Loaded genotype matrix: ", nrow(geno), " samples x ", ncol(geno), " SNPs before chromosome filtering.")
  }

  if (!is.null(null_sample_ids)) {
    missing_null_ids <- setdiff(null_sample_ids, rownames(geno))
    if (length(missing_null_ids) > 0L) {
      stop(
        "Some NULL-model samples are missing from genotype input: ",
        paste(utils::head(missing_null_ids, 5), collapse = ", ")
      )
    }
    geno <- geno[null_sample_ids, , drop = FALSE]
    if (!is.null(cluster)) {
      cluster <- cluster[null_sample_ids]
    }
    message("Genotype matrix after aligning to NULL model samples: ", nrow(geno), " samples x ", ncol(geno), " SNPs.")
  }

  # Subset for quick testing (every 10th SNP)
  # geno <- geno[, seq(1, ncol(geno), by = 10), drop = FALSE]

  # --------------------------
  # Subset by chromosome if specified
  # --------------------------
  if (!is.null(chrom) && toupper(chrom) != "NULL") {
    message("Subsetting genotype data for chromosome: ", chrom)
    chrom <- sub("^chr", "", chrom, ignore.case = TRUE)

    chr_map <- sub("^chr", "", chr_map, ignore.case = TRUE)
    keep <- which(chr_map == chrom)

    if (length(keep) == 0L) stop("No SNPs found for --chrom=", chrom)
    geno <- geno[, keep, drop = FALSE]
    message("Genotype matrix after chromosome filtering: ", nrow(geno), " samples x ", ncol(geno), " SNPs.")
  }

  if (!is.null(outputSnpFile) && nzchar(outputSnpFile)) {
    message(
      "Generating SNPInfo from the Step2 genotype matrix after NULL-model sample alignment",
      if (!is.null(chrom) && toupper(chrom) != "NULL") " and chromosome filtering" else "",
      ", before minMAF/minMAC filtering. Output path: ",
      outputSnpFile
    )
    snp_stats <- snp_info_from_geno_matrix(geno, snp_ids = colnames(geno))
    dir.create(dirname(outputSnpFile), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(
      snp_stats,
      file = outputSnpFile,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )
    message("SNPInfo finished: ", nrow(snp_stats), " SNP(s) written to ", outputSnpFile)
  } else {
    message("SNPInfo skipped: outputSnpFile is NULL.")
  }

  if (minMAF > 0 || minMAC > 0) {
    allele_freq <- colMeans(geno, na.rm = TRUE) / 2
    maf <- pmin(allele_freq, 1 - allele_freq)
    called_alleles <- 2 * colSums(!is.na(geno))
    mac <- pmin(colSums(geno, na.rm = TRUE), called_alleles - colSums(geno, na.rm = TRUE))
    keep <- rep(TRUE, ncol(geno))
    if (minMAF > 0) {
      keep <- keep & !is.na(maf) & maf >= minMAF
    }
    if (minMAC > 0) {
      keep <- keep & !is.na(mac) & mac >= minMAC
    }
    keep <- which(keep)
    if (length(keep) == 0L) {
      stop("No SNPs remain after applying SNP filters (--minMAF=", minMAF, ", --minMAC=", minMAC, ").")
    }
    message(
      "Genotype matrix after MAF/MAC filtering: ", nrow(geno), " samples x ", length(keep),
      " SNPs (removed ", ncol(geno) - length(keep), ")."
    )
    geno <- geno[, keep, drop = FALSE]
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
  message("Writing ", nrow(res), " per-phenotype result file(s).")

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
    out_paths <- c(out_paths, out_file)
  }

  message("Done. Wrote ", length(out_paths), " file(s) with prefix: ", PALMOutputFile)
  invisible(out_paths)
}
