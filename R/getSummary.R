#' Run PALM summary statistics for genotype inputs
#'
#' Reads genotype data from PLINK, VCF, or BGEN input, applies a pre-fitted
#' PALM null model, and writes per-phenotype summary statistics (one file per
#' phenotype). Non-PLINK inputs are converted to a temporary PLINK dataset
#' before testing.
#'
#' @param genoFile Path to genotype input. If the path ends with
#'   `.vcf`, `.vcf.gz`, or `.vcf.bgz`, VCF input is assumed. If it ends with
#'   `.bgen`, BGEN input is assumed. Otherwise it is treated as a PLINK prefix
#'   without extension.
#' @param NULLmodelFile Path to `.rda` containing the fitted null model
#'   object named `modglmm`.
#' @param PALMOutputFile Output prefix for per-phenotype result files; each
#'   phenotype file will be written as `<PALMOutputFile>_<pheno>.txt`.
#' @param vcfField VCF FORMAT field to import when converting VCF input.
#'   Supported values are `"DS"` and `"GT"`.
#' @param alleleOrder Allele order for BGEN conversion. Supported values are
#'   `"ref-first"`, `"ref-last"`, and `"ref-unknown"`. Defaults to
#'   `"ref-last"` for BGEN so imported allele coding matches the existing
#'   PLINK-based step2 convention.
#' @param plinkPath Path to the `plink` executable used for PLINK/VCF
#'   conversion fallback. Defaults to `"plink"`.
#' @param plink2Path Path to the `plink2` executable used for BGEN conversion
#'   and preferred VCF conversion. Defaults to `"plink2"`.
#' @param keepTemp Logical; if `TRUE`, keep temporary converted PLINK files for
#'   VCF/BGEN inputs. Defaults to `FALSE`.
#' @param chrom Optional chromosome filter (numeric or string like `"chr1"`).
#'   Use `NULL` to keep all chromosomes.
#' @param minMAF Optional minimum minor allele frequency threshold used to
#'   filter SNPs before running `PALM::palm.get.summary()`. Use `0` (default)
#'   to disable MAF filtering.
#' @param correct Passed to `PALM::palm.get.summary()`; defaults to `"NULL"`.
#' @param useCluster Logical; if `TRUE`, uses FID from `.fam` as cluster
#'   information when available.
#'
#' @return Invisibly returns a character vector of written file paths.
#' @export
getSummary <- function(genoFile,
                       NULLmodelFile,
                       PALMOutputFile,
                       vcfField = "DS",
                       alleleOrder = NULL,
                       plinkPath = "plink",
                       plink2Path = "plink2",
                       keepTemp = FALSE,
                       chrom = NULL,
                       minMAF = 0,
                       correct = "NULL",
                       useCluster = TRUE) {
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
  if (!is.logical(keepTemp) || length(keepTemp) != 1L || is.na(keepTemp)) {
    stop("'keepTemp' must be TRUE or FALSE.")
  }
  output_dir <- dirname(PALMOutputFile)
  if (!output_dir %in% c("", ".")) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  infer_geno_format <- function(path) {
    lower <- tolower(path)
    if (grepl("\\.(vcf|vcf\\.gz|vcf\\.bgz)$", lower)) {
      return("vcf")
    }
    if (grepl("\\.bgen$", lower)) {
      return("bgen")
    }
    if (file.exists(paste0(path, ".bed")) &&
        file.exists(paste0(path, ".bim")) &&
        file.exists(paste0(path, ".fam"))) {
      return("plink")
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
    message("Running command: ", command, " ", paste(args, collapse = " "))
    status <- system2(command, args = args)
    if (!identical(status, 0L)) {
      stop("Command failed with exit status ", status, ": ", command)
    }
  }

  prepare_plink_input <- function(genoFile,
                                  genoFormat,
                                  vcfField,
                                  alleleOrder,
                                  plinkPath,
                                  plink2Path,
                                  keepTemp,
                                  PALMOutputFile) {
    if (identical(genoFormat, "plink")) {
      bed <- paste0(genoFile, ".bed")
      bim <- paste0(genoFile, ".bim")
      fam <- paste0(genoFile, ".fam")
      for (f in c(bed, bim, fam)) {
        if (!file.exists(f)) {
          stop("Missing PLINK file: ", f)
        }
      }
      return(list(prefix = genoFile, format = "plink", cleanup = character(0)))
    }

    if (!file.exists(genoFile)) {
      stop("Genotype input file not found: ", genoFile)
    }

    temp_prefix <- file.path(
      dirname(PALMOutputFile),
      paste0(
        basename(PALMOutputFile),
        "_tmp_",
        genoFormat,
        "_plink"
      )
    )
    cleanup <- paste0(temp_prefix, c(".bed", ".bim", ".fam", ".log", ".nosex"))

    if (identical(genoFormat, "vcf")) {
      vcfField <- toupper(trimws(as.character(vcfField)[1]))
      if (!vcfField %in% c("DS", "GT")) {
        stop("'vcfField' must be either 'DS' or 'GT'.")
      }

      plink2_exec <- find_executable(plink2Path)
      if (nzchar(plink2_exec)) {
        args <- c("--vcf", genoFile)
        if (identical(vcfField, "DS")) {
          args <- c(args, "dosage=DS")
        }
        args <- c(args, "--double-id", "--keep-allele-order", "--make-bed", "--out", temp_prefix)
        run_command_checked(plink2_exec, args)
      } else {
        plink_exec <- find_executable(plinkPath)
        if (!nzchar(plink_exec)) {
          stop("Neither '", plink2Path, "' nor '", plinkPath, "' was found in PATH; cannot convert VCF input.")
        }
        args <- c("--vcf", genoFile)
        if (identical(vcfField, "DS")) {
          args <- c(args, "dosage=DS")
        }
        args <- c(args, "--double-id", "--keep-allele-order", "--make-bed", "--out", temp_prefix)
        run_command_checked(plink_exec, args)
      }
    } else if (identical(genoFormat, "bgen")) {
      plink2_exec <- find_executable(plink2Path)
      if (!nzchar(plink2_exec)) {
        stop("BGEN input requires plink2, but executable '", plink2Path, "' was not found in PATH.")
      }
      if (is.null(alleleOrder) || !nzchar(as.character(alleleOrder)[1]) || toupper(as.character(alleleOrder)[1]) == "NULL") {
        alleleOrder <- "ref-last"
      }
      alleleOrder <- trimws(as.character(alleleOrder)[1])
      if (!alleleOrder %in% c("ref-first", "ref-last", "ref-unknown")) {
        stop("'alleleOrder' must be 'ref-first', 'ref-last', or 'ref-unknown' for BGEN input.")
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
        message("Detected BGEN sample file: ", sample_file)
      } else {
        message("No BGEN sample file detected; plink2 will rely on sample IDs embedded in the .bgen file.")
      }

      bgi_candidates <- c(
        sub("\\.bgen$", ".bgi", genoFile, ignore.case = TRUE),
        paste0(genoFile, ".bgi")
      )
      bgi_candidates <- unique(bgi_candidates[file.exists(bgi_candidates)])
      if (length(bgi_candidates) > 0L) {
        message("Detected BGEN index file: ", bgi_candidates[1])
      } else {
        message("No standalone BGEN .bgi file detected next to input; relying on plink2 defaults.")
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
      message("Keeping temporary converted PLINK files with prefix: ", temp_prefix)
    }

    list(prefix = temp_prefix, format = genoFormat, cleanup = cleanup)
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
  if (is.null(correct)) {
    message("Compositional correction disabled: correct=NULL")
  } else {
    message("Compositional correction enabled: correct=", correct)
  }
  message("Cluster option requested: useCluster=", useCluster)

  if (!identical(genoFormat, "plink") && isTRUE(useCluster)) {
    message("Clustering from PLINK .fam FID is only available for native PLINK input. Setting useCluster=FALSE.")
    useCluster <- FALSE
  }

  geno_input <- prepare_plink_input(
    genoFile = genoFile,
    genoFormat = genoFormat,
    vcfField = vcfField,
    alleleOrder = alleleOrder,
    plinkPath = plinkPath,
    plink2Path = plink2Path,
    keepTemp = keepTemp,
    PALMOutputFile = PALMOutputFile
  )
  if (!keepTemp && length(geno_input$cleanup) > 0L) {
    on.exit(unlink(geno_input$cleanup, force = TRUE), add = TRUE)
  }
  genoPrefix <- geno_input$prefix

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

  if (minMAF > 0) {
    allele_freq <- colMeans(geno, na.rm = TRUE) / 2
    maf <- pmin(allele_freq, 1 - allele_freq)
    keep <- which(!is.na(maf) & maf >= minMAF)
    if (length(keep) == 0L) {
      stop("No SNPs remain after applying --minMAF=", minMAF)
    }
    message(
      "Genotype matrix after MAF filtering: ", nrow(geno), " samples x ", length(keep),
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
