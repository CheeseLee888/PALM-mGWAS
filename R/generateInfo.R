# step0 information helpers

#' Compute feature-level prevalence and average proportion
#'
#' @param abd_file Path to abundance table (tab-separated, first column sample ID).
#' @return A data frame with FeatureID, Prevalence, AvgProportion.
#' @export
feature_info <- function(abd_file) {
  if (missing(abd_file) || is.null(abd_file) || !nzchar(abd_file)) {
    stop("`abd_file` is required")
  }
  if (!file.exists(abd_file)) {
    stop("abd_file does not exist: ", abd_file)
  }

  abd_data <- utils::read.table(abd_file, header = TRUE, sep = "\t")
  if (ncol(abd_data) < 2) {
    stop("abd_file must contain at least one feature column")
  }

  gene_mat <- as.matrix(abd_data[, -1, drop = FALSE])
  gene_mat <- apply(gene_mat, 2, as.numeric)
  n_sample <- nrow(gene_mat)

  prevalence <- colSums(gene_mat > 0, na.rm = TRUE) / n_sample

  row_sum <- rowSums(gene_mat, na.rm = TRUE)
  row_sum[row_sum == 0] <- NA
  gene_prop <- gene_mat / row_sum
  avg_proportion <- colMeans(gene_prop, na.rm = TRUE)

  data.frame(
    FeatureID = colnames(gene_mat),
    Prevalence = prevalence,
    AvgProportion = avg_proportion,
    row.names = NULL,
    check.names = FALSE
  )
}


#' Compute sequencing depth per sample
#'
#' @param abd_file Path to abundance table (tab-separated, first column sample ID).
#' @return A data frame with SampleID and SeqDepth.
#' @export
seqdepth_info <- function(abd_file) {
  if (missing(abd_file) || is.null(abd_file) || !nzchar(abd_file)) {
    stop("`abd_file` is required")
  }
  if (!file.exists(abd_file)) {
    stop("abd_file does not exist: ", abd_file)
  }

  abd <- utils::read.table(abd_file, header = TRUE, sep = "\t")
  if (ncol(abd) < 2) {
    stop("abd_file must contain at least one feature column")
  }

  sample_id <- abd[[1]]
  abd_mat <- as.matrix(abd[, -1, drop = FALSE])
  rownames(abd_mat) <- sample_id

  data.frame(
    SampleID = sample_id,
    SeqDepth = rowSums(abd_mat, na.rm = TRUE),
    row.names = NULL,
    check.names = FALSE
  )
}


#' Parse SNP IDs of the form "chr1:1119172:G:A"
#'
#' @param snp_ids Character vector of SNP IDs.
#' @return Data frame with CHR, SNP, POS, A1, A2.
#' @export
parse_snp_ids <- function(snp_ids) {
  parts <- strsplit(snp_ids, ":", fixed = TRUE)
  lens <- vapply(parts, length, integer(1))
  if (any(lens < 4)) {
    bad <- snp_ids[lens < 4][1]
    stop("SNP id format error. Expect 'chr<CHR>:<POS>:<A1>:<A2>'. Bad SNP: ", bad)
  }

  chr_raw <- vapply(parts, `[[`, character(1), 1)
  pos_raw <- vapply(parts, `[[`, character(1), 2)
  a1 <- vapply(parts, `[[`, character(1), 4)
  a2 <- vapply(parts, `[[`, character(1), 3)

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


#' Compute SNP counts and allele frequency from genotype input
#'
#' Accepts a PLINK prefix directly, or converts VCF/BGEN input to a temporary
#' PLINK dataset before reading SNP counts and allele frequencies.
#'
#' @param genoFile Genotype input path. If the value ends with `.vcf`,
#'   `.vcf.gz`, or `.vcf.bgz`, VCF input is assumed. If it ends with `.bgen`,
#'   BGEN input is assumed. Otherwise it is treated as a PLINK prefix without
#'   `.bed/.bim/.fam`.
#' @param vcf_field VCF FORMAT field to import when converting VCF input.
#'   Supported values are `"DS"` and `"GT"`. Defaults to `"DS"`.
#' @param allele_order Allele order for BGEN conversion. Supported values are
#'   `"ref-first"`, `"ref-last"`, and `"ref-unknown"`. Defaults to
#'   `"ref-last"` for BGEN input.
#' @param plink_path Path to the `plink` executable used as a VCF conversion
#'   fallback. Defaults to `"plink"`.
#' @param plink2_path Path to the `plink2` executable used for BGEN conversion
#'   and preferred VCF conversion. Defaults to `"plink2"`.
#' @param keep_temp Logical; if `TRUE`, keep temporary converted PLINK files for
#'   VCF/BGEN input. Defaults to `FALSE`.
#' @return Data frame with CHR, SNP, POS, A1, A2, N, AF.
#' @export
snp_info <- function(genoFile,
                     vcf_field = "DS",
                     allele_order = NULL,
                     plink_path = "plink",
                     plink2_path = "plink2",
                     keep_temp = FALSE) {
  if (missing(genoFile) || is.null(genoFile) || !nzchar(genoFile)) {
    stop("`genoFile` is required")
  }
  if (!is.logical(keep_temp) || length(keep_temp) != 1L || is.na(keep_temp)) {
    stop("`keep_temp` must be TRUE or FALSE.")
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

  prepare_plink_input <- function(geno_path,
                                  geno_format,
                                  vcf_field,
                                  allele_order,
                                  plink_path,
                                  plink2_path,
                                  keep_temp) {
    if (identical(geno_format, "plink")) {
      bed <- paste0(geno_path, ".bed")
      bim <- paste0(geno_path, ".bim")
      fam <- paste0(geno_path, ".fam")
      for (f in c(bed, bim, fam)) {
        if (!file.exists(f)) {
          stop("Missing PLINK file: ", f)
        }
      }
      return(list(prefix = geno_path, cleanup = character(0)))
    }

    if (!file.exists(geno_path)) {
      stop("Genotype input file not found: ", geno_path)
    }

    temp_prefix <- file.path(
      tempdir(),
      paste0(
        sanitize_local_name(basename(geno_path)),
        "_info_tmp_",
        geno_format,
        "_plink"
      )
    )
    cleanup <- paste0(temp_prefix, c(".bed", ".bim", ".fam", ".log", ".nosex"))

    if (identical(geno_format, "vcf")) {
      vcf_field <- toupper(trimws(as.character(vcf_field)[1]))
      if (!vcf_field %in% c("DS", "GT")) {
        stop("`vcf_field` must be either 'DS' or 'GT'.")
      }

      plink2_exec <- find_executable(plink2_path)
      if (nzchar(plink2_exec)) {
        args <- c("--vcf", geno_path)
        if (identical(vcf_field, "DS")) {
          args <- c(args, "dosage=DS")
        }
        args <- c(args, "--double-id", "--keep-allele-order", "--make-bed", "--out", temp_prefix)
        run_command_checked(plink2_exec, args)
      } else {
        plink_exec <- find_executable(plink_path)
        if (!nzchar(plink_exec)) {
          stop("Neither '", plink2_path, "' nor '", plink_path, "' was found in PATH; cannot convert VCF input.")
        }
        args <- c("--vcf", geno_path)
        if (identical(vcf_field, "DS")) {
          args <- c(args, "dosage=DS")
        }
        args <- c(args, "--double-id", "--keep-allele-order", "--make-bed", "--out", temp_prefix)
        run_command_checked(plink_exec, args)
      }
    } else if (identical(geno_format, "bgen")) {
      plink2_exec <- find_executable(plink2_path)
      if (!nzchar(plink2_exec)) {
        stop("BGEN input requires plink2, but executable '", plink2_path, "' was not found in PATH.")
      }
      if (is.null(allele_order) || !nzchar(as.character(allele_order)[1]) || toupper(as.character(allele_order)[1]) == "NULL") {
        allele_order <- "ref-last"
      }
      allele_order <- trimws(as.character(allele_order)[1])
      if (!allele_order %in% c("ref-first", "ref-last", "ref-unknown")) {
        stop("`allele_order` must be 'ref-first', 'ref-last', or 'ref-unknown' for BGEN input.")
      }

      sample_candidates <- c(
        sub("\\.bgen$", ".sample", geno_path, ignore.case = TRUE),
        sub("\\.bgen$", ".samples", geno_path, ignore.case = TRUE),
        paste0(geno_path, ".sample"),
        paste0(geno_path, ".samples")
      )
      sample_candidates <- unique(sample_candidates[file.exists(sample_candidates)])
      sample_file <- character(0)
      if (length(sample_candidates) > 0L) {
        sample_file <- sample_candidates[1]
        message("Detected BGEN sample file: ", sample_file)
      } else {
        message("No BGEN sample file detected; plink2 will rely on sample IDs embedded in the .bgen file.")
      }

      args <- c("--bgen", geno_path, allele_order)
      if (length(sample_file) == 1L && nzchar(sample_file)) {
        args <- c(args, "--sample", sample_file)
      }
      args <- c(args, "--make-bed", "--out", temp_prefix)
      run_command_checked(plink2_exec, args)
    } else {
      stop("Unsupported inferred genotype format: ", geno_format)
    }

    for (f in paste0(temp_prefix, c(".bed", ".bim", ".fam"))) {
      if (!file.exists(f)) {
        stop("Conversion to temporary PLINK files failed; missing output: ", f)
      }
    }

    if (keep_temp) {
      message("Keeping temporary converted PLINK files with prefix: ", temp_prefix)
    }

    list(prefix = temp_prefix, cleanup = cleanup)
  }

  geno_format <- infer_geno_format(genoFile)
  message("Genotype input format inferred from input path: ", geno_format)
  geno_input <- prepare_plink_input(
    geno_path = genoFile,
    geno_format = geno_format,
    vcf_field = vcf_field,
    allele_order = allele_order,
    plink_path = plink_path,
    plink2_path = plink2_path,
    keep_temp = keep_temp
  )
  if (!keep_temp && length(geno_input$cleanup) > 0L) {
    on.exit(unlink(geno_input$cleanup, force = TRUE), add = TRUE)
  }

  bed <- paste0(geno_input$prefix, ".bed")
  bim <- paste0(geno_input$prefix, ".bim")
  fam <- paste0(geno_input$prefix, ".fam")
  for (f in c(bed, bim, fam)) {
    if (!file.exists(f)) {
      stop("Missing PLINK file: ", f)
    }
  }

  plink <- snpStats::read.plink(bed, bim, fam)

  G <- as(plink$genotypes, "numeric")
  snps <- colnames(G)
  if (is.null(snps)) {
    stop("No SNP names found in PLINK genotypes.")
  }

  N <- colSums(!is.na(G))
  MAC <- colSums(G, na.rm = TRUE)
  AF <- rep(NA_real_, length(snps))
  ok <- N > 0
  AF[ok] <- (MAC[ok] / N[ok]) / 2

  info <- parse_snp_ids(snps)
  info$N <- as.integer(N)
  info$AF <- as.numeric(AF)

  info[, c("CHR", "SNP", "POS", "A1", "A2", "N", "AF")]
}


#' Run all step0 info generation and write outputs
#'
#' @param genoFile Genotype input path. Accepts a PLINK prefix directly, or a
#'   VCF/BGEN file that will be converted to a temporary PLINK dataset.
#' @param abd_file Path to abundance table.
#' @param output_snp Output path for SNP info.
#' @param output_feature Output path for feature info.
#' @param output_seqdepth Output path for sequencing depth info.
#' @param vcf_field VCF FORMAT field to import when `genoFile` is VCF input.
#'   Defaults to `"DS"`.
#' @param allele_order Allele order for BGEN conversion. Defaults to
#'   `"ref-last"` when `genoFile` is BGEN input.
#' @param keep_temp Logical; if `TRUE`, keep temporary PLINK files created for
#'   VCF/BGEN inputs. Defaults to `FALSE`.
#' @return Invisibly returns a list of output paths.
#' @export
run_info <- function(genoFile,
                     abd_file,
                     output_snp,
                     output_feature,
                     output_seqdepth,
                     vcf_field = "DS",
                     allele_order = NULL,
                     keep_temp = FALSE) {
  missing_args <- vapply(
    list(genoFile, abd_file, output_snp, output_feature, output_seqdepth),
    function(x) is.null(x) || !nzchar(x),
    logical(1)
  )
  if (any(missing_args)) {
    stop("All input and output paths must be provided.")
  }

  feature <- feature_info(abd_file)
  seqdepth <- seqdepth_info(abd_file)
  snp <- snp_info(
    genoFile = genoFile,
    vcf_field = vcf_field,
    allele_order = allele_order,
    keep_temp = keep_temp
  )

  dir.create(dirname(output_snp), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(output_feature), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(output_seqdepth), recursive = TRUE, showWarnings = FALSE)

  data.table::fwrite(snp, file = output_snp, sep = "\t", quote = FALSE, na = "NA")
  data.table::fwrite(feature, file = output_feature, sep = "\t", quote = FALSE, na = "NA")
  data.table::fwrite(seqdepth, file = output_seqdepth, sep = "\t", quote = FALSE, na = "NA")

  invisible(list(
    snp = output_snp,
    feature = output_feature,
    seqdepth = output_seqdepth
  ))
}
