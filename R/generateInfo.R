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
  geno_input <- prepare_plink_input(
    genoFile = genoFile,
    vcfField = vcf_field,
    alleleOrder = allele_order,
    plinkPath = plink_path,
    plink2Path = plink2_path,
    keepTemp = keep_temp,
    tempLabel = "info_tmp"
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
#' @param plink_path Path to the `plink` executable used as a VCF conversion
#'   fallback. Defaults to `"plink"`.
#' @param plink2_path Path to the `plink2` executable used for BGEN conversion
#'   and preferred VCF conversion. Defaults to `"plink2"`.
#' @return Invisibly returns a list of output paths.
#' @export
run_info <- function(genoFile,
                     abd_file,
                     output_snp,
                     output_feature,
                     output_seqdepth,
                     vcf_field = "DS",
                     allele_order = NULL,
                     keep_temp = FALSE,
                     plink_path = "plink",
                     plink2_path = "plink2") {
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
    plink_path = plink_path,
    plink2_path = plink2_path,
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
