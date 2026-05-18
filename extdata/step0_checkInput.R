#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PALMGWAS)
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--abdFile",  type = "character",
              help = "Abundance table path"),
  make_option("--covFile",  type = "character",
              help = "Covariate table path"),
  make_option("--abdAlignedFile", type = "character",
              help = "Output path for aligned abundance table"),
  make_option("--covAlignedFile", type = "character",
              help = "Output path for aligned covariate table"),
  make_option("--covarColList", type = "character", default = "NULL",
              help = "Optional comma-separated covariate columns used in Step1; if NULL, all non-ID covariate columns are used. Samples missing these columns are removed [default %default]"),
  make_option("--depthCol", type = "character", default = "NULL",
              help = "Optional covariate column name used as sequencing depth [default %default]"),
  make_option("--timeIDCol", type = "character", default = "NULL",
              help = "Optional time ID column required when abdFile contains repeated subject IDs [default %default]"),
  make_option("--clusterCol", type = "character", default = "NULL",
              help = "Optional covFile column used as family/pedigree cluster ID [default %default]"),
  make_option("--depth.filter", type = "double", default = 0,
              help = "Row-level depth threshold; rows with depth <= threshold are removed before ID matching [default %default]"),
  make_option("--genoFile", type = "character",
              help = "Genotype input: PLINK prefix or VCF(.vcf/.vcf.gz/.vcf.bgz)"),
  make_option("--seqdepthInfoFile", type = "character", default = "NULL",
              help = "Optional output file for sequencing depth info used by Step0 filtering [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$abdFile) || is.null(opt$covFile) || is.null(opt$genoFile) ||
    is.null(opt$abdAlignedFile) || is.null(opt$covAlignedFile)) {
  stop("abdFile, covFile, genoFile, abdAlignedFile, and covAlignedFile must all be provided.")
}

normalize_col_list <- function(x, arg_name) {
  if (is.null(x) || !nzchar(x) || toupper(x) == "NULL") {
    return(NULL)
  }
  cols <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  cols <- cols[nzchar(cols)]
  if (!length(cols)) {
    stop("'", arg_name, "' did not contain any valid column names.")
  }
  cols
}

default_covariate_cols <- function(df, time_col = NULL) {
  cols <- colnames(df)
  if (length(cols) <= 1L) {
    return(character(0))
  }
  setdiff(cols[-1], time_col)
}

read_input_table <- function(path, label) {
  df <- fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(df) < 1) stop("File has no columns: ", path)
  ids <- as.character(df[[1]])
  if (anyNA(ids) || any(!nzchar(ids))) stop("Missing/empty subject ID detected in ", label, " file: ", path)
  df
}

make_pair_key <- function(subject_id, time_id) {
  paste(subject_id, time_id, sep = "\r")
}

validate_time_col <- function(df, time_col, label) {
  if (is.null(time_col) || !(time_col %in% colnames(df))) {
    stop(label, " contains repeated subject IDs, so --timeIDCol must name a time ID column present in that file.")
  }
  time_id <- as.character(df[[time_col]])
  if (anyNA(time_id) || any(!nzchar(time_id))) {
    stop("Missing/empty time ID detected in ", label, " column '", time_col, "'.")
  }
  subject_id <- as.character(df[[1]])
  key <- make_pair_key(subject_id, time_id)
  duplicated_key <- unique(key[duplicated(key)])
  if (length(duplicated_key) > 0L) {
    stop(label, " contains duplicated (subject ID, time ID) rows.")
  }
  key
}

match_longitudinal_covariates <- function(abd_df, cov_df, time_col) {
  abd_ids <- as.character(abd_df[[1]])
  cov_ids <- as.character(cov_df[[1]])
  abd_repeated <- anyDuplicated(abd_ids) > 0L
  cov_repeated <- anyDuplicated(cov_ids) > 0L

  if (!abd_repeated) {
    if (cov_repeated) {
      stop("covFile has repeated subject IDs but abdFile does not; cannot infer which covariate row to use.")
    }
    ord <- match(abd_ids, cov_ids)
    if (anyNA(ord)) {
      stop("Some abundance subject IDs are missing from covFile: ", paste(utils::head(abd_ids[is.na(ord)], 5), collapse = ", "))
    }
    return(list(abd = abd_df, cov = cov_df[ord, , drop = FALSE], repeated = FALSE))
  }

  abd_key <- validate_time_col(abd_df, time_col, "abdFile")
  if (cov_repeated) {
    cov_key <- validate_time_col(cov_df, time_col, "covFile")
    missing_cov_key <- setdiff(abd_key, cov_key)
    extra_cov_key <- setdiff(cov_key, abd_key)
    if (length(missing_cov_key) > 0L) {
      stop(
        "covFile is missing repeated-measure covariate rows for ",
        length(missing_cov_key),
        " abundance (subject ID, time ID) pair(s)."
      )
    }
    if (length(extra_cov_key) > 0L) {
      message("Dropping ", length(extra_cov_key), " covariate-only (subject ID, time ID) row(s).")
    }
    ord <- match(abd_key, cov_key)
    return(list(abd = abd_df, cov = cov_df[ord, , drop = FALSE], repeated = TRUE))
  }

  duplicated_cov <- unique(cov_ids[duplicated(cov_ids)])
  if (length(duplicated_cov) > 0L) {
    stop("Internal error: covariate subject-level branch received duplicated subject IDs.")
  }
  missing_subject_cov <- setdiff(unique(abd_ids), cov_ids)
  if (length(missing_subject_cov) > 0L) {
    stop("Some abundance subject IDs are missing from subject-level covFile: ", paste(utils::head(missing_subject_cov, 5), collapse = ", "))
  }
  ord <- match(abd_ids, cov_ids)
  list(abd = abd_df, cov = cov_df[ord, , drop = FALSE], repeated = TRUE)
}

read_fam_iid <- function(prefix) {
  fam <- paste0(prefix, ".fam")
  if (!file.exists(fam)) stop("Missing .fam file: ", fam)
  fam_df <- fread(fam, data.table = FALSE, header = FALSE)
  if (ncol(fam_df) < 2) stop("Invalid .fam (need >=2 cols): ", fam)
  ids <- as.character(fam_df[[2]])  # IID
  if (anyDuplicated(ids)) stop("Duplicated IID in .fam: ", fam)
  ids
}

read_fam_cluster <- function(prefix) {
  fam <- paste0(prefix, ".fam")
  if (!file.exists(fam)) stop("Missing .fam file: ", fam)
  fam_df <- fread(fam, data.table = FALSE, header = FALSE)
  if (ncol(fam_df) < 2) stop("Invalid .fam (need >=2 cols): ", fam)
  ids <- as.character(fam_df[[2]])
  if (anyNA(ids) || any(!nzchar(ids))) stop("Missing/empty IID detected in .fam: ", fam)
  if (anyDuplicated(ids)) stop("Duplicated IID in .fam: ", fam)
  cluster <- as.character(fam_df[[1]])
  names(cluster) <- ids
  cluster
}

read_geno_iid <- function(geno_file) {
  geno_format <- PALMGWAS:::infer_geno_format(geno_file)
  if (identical(geno_format, "vcf")) {
    ids <- PALMGWAS:::read_vcf_header(geno_file)$sample_ids
    if (anyDuplicated(ids)) stop("Duplicated IID in VCF header: ", geno_file)
    return(list(ids = ids, cleanup = character(0)))
  }

  geno_input <- PALMGWAS:::prepare_plink_input(
    genoFile = geno_file,
    tempLabel = "check_tmp"
  )
  list(ids = read_fam_iid(geno_input$prefix), cleanup = geno_input$cleanup)
}

order_rows_by_subject_ref <- function(ids, ref_ids) {
  unlist(lapply(ref_ids, function(id) which(ids == id)), use.names = FALSE)
}

cat("Checking IIDs...\n")

abd_df <- read_input_table(opt$abdFile, "abundance")
cov_df <- read_input_table(opt$covFile, "covariate")
cat(
  "Loaded input tables: abd=", nrow(abd_df),
  " row(s), cov=", nrow(cov_df), " row(s).\n",
  sep = ""
)
if (is.null(opt$seqdepthInfoFile) || !nzchar(opt$seqdepthInfoFile) || toupper(opt$seqdepthInfoFile) == "NULL") {
  opt$seqdepthInfoFile <- NULL
}
opt$covarColList <- normalize_col_list(opt$covarColList, "covarColList")
opt$depthCol <- normalize_col_list(opt$depthCol, "depthCol")
if (!is.null(opt$depthCol) && length(opt$depthCol) != 1L) {
  stop("'depthCol' must specify exactly one column name.")
}
opt$timeIDCol <- normalize_col_list(opt$timeIDCol, "timeIDCol")
if (!is.null(opt$timeIDCol) && length(opt$timeIDCol) != 1L) {
  stop("'timeIDCol' must specify exactly one column name.")
}
opt$clusterCol <- normalize_col_list(opt$clusterCol, "clusterCol")
if (!is.null(opt$clusterCol) && length(opt$clusterCol) != 1L) {
  stop("'clusterCol' must specify exactly one column name.")
}
if (!is.numeric(opt$depth.filter) || length(opt$depth.filter) != 1L || is.na(opt$depth.filter) || opt$depth.filter < 0) {
  stop("--depth.filter must be a single non-negative numeric value.")
}

if (is.null(opt$covarColList)) {
  opt$covarColList <- default_covariate_cols(cov_df, opt$timeIDCol)
  if (length(opt$covarColList) > 0L) {
    cat("covarColList not provided: defaulting to all covariate columns in covFile.\n")
  } else {
    cat("covarColList not provided: covFile has no non-ID covariate columns.\n")
  }
}

required_cov_cols <- unique(c(opt$covarColList, opt$depthCol, opt$clusterCol))
missing_cov_cols <- setdiff(required_cov_cols, colnames(cov_df))
if (length(missing_cov_cols) > 0L) {
  stop("Required covariate column(s) not found in covFile: ", paste(missing_cov_cols, collapse = ", "))
}
if (length(required_cov_cols) > 0L) {
  cat("Step1-required covariate/depth columns: ", paste(required_cov_cols, collapse = ", "), "\n", sep = "")
  missing_counts <- colSums(is.na(cov_df[, required_cov_cols, drop = FALSE]))
  missing_summary <- paste(names(missing_counts), missing_counts, sep = "=")
  cat("Missing counts by required covariate/depth column: ", paste(missing_summary, collapse = ", "), "\n", sep = "")
} else {
  cat("No Step1-required covariate/depth columns requested for missingness filtering.\n")
}

filtered <- FALSE

abd_non_id_cols <- setdiff(colnames(abd_df), c(colnames(abd_df)[1], opt$timeIDCol))
if (length(abd_non_id_cols) == 0L) {
  stop("abdFile must contain at least one numeric abundance feature column after the subject ID column.")
}
if (length(abd_non_id_cols) > 0L) {
  bad_numeric_cols <- character(0)
  for (col in abd_non_id_cols) {
    raw_values <- abd_df[[col]]
    numeric_values <- suppressWarnings(as.numeric(raw_values))
    bad_values <- is.na(numeric_values) & !is.na(raw_values)
    if (any(bad_values)) {
      bad_numeric_cols <- c(bad_numeric_cols, col)
    } else {
      abd_df[[col]] <- numeric_values
    }
  }
  if (length(bad_numeric_cols) > 0L) {
    stop("Abundance feature column(s) in abdFile cannot be safely converted to numeric: ", paste(unique(bad_numeric_cols), collapse = ", "))
  }
  keep_abd_complete <- stats::complete.cases(abd_df[, abd_non_id_cols, drop = FALSE])
  removed_n <- sum(!keep_abd_complete)
  cat("Removing ", removed_n, " abundance row(s) with missing values in abundance columns.\n", sep = "")
  if (removed_n > 0L) {
    abd_df <- abd_df[keep_abd_complete, , drop = FALSE]
    filtered <- TRUE
  }
  if (nrow(abd_df) == 0L) {
    stop("No abundance rows remain after removing rows with missing abundance values.")
  }
}

abd_repeated_subjects <- anyDuplicated(as.character(abd_df[[1]])) > 0L
cat("Abundance repeated subject IDs after abundance filtering: ", abd_repeated_subjects, ".\n", sep = "")
matched <- match_longitudinal_covariates(abd_df, cov_df, opt$timeIDCol)
abd_df <- matched$abd
cov_df <- matched$cov
cat(
  if (matched$repeated) {
    "Matched covariates to longitudinal abundance rows.\n"
  } else {
    "Matched covariates to one row per subject.\n"
  }
)

if (length(required_cov_cols) > 0L) {
  keep_cov_complete <- stats::complete.cases(cov_df[, required_cov_cols, drop = FALSE])
  removed_n <- sum(!keep_cov_complete)
  cat("Removing ", removed_n, " row(s) with missing values in Step1-required covariates/depth columns.\n", sep = "")
  if (removed_n > 0L) {
    cov_df <- cov_df[keep_cov_complete, , drop = FALSE]
    abd_df <- abd_df[keep_cov_complete, , drop = FALSE]
    filtered <- TRUE
  }
  if (nrow(cov_df) == 0L || nrow(abd_df) == 0L) {
    stop("No rows remain after removing rows with missing Step1-required covariate/depth values.")
  }
}

abd_ids_all <- as.character(abd_df[[1]])
if (is.null(opt$depthCol)) {
  depth <- rowSums(as.matrix(abd_df[, abd_non_id_cols, drop = FALSE]), na.rm = TRUE)
  names(depth) <- abd_ids_all
  cat("Depth filter source: row sums of abdFile.\n")
} else {
  depth <- suppressWarnings(as.numeric(cov_df[[opt$depthCol]]))
  if (any(is.na(depth) & !is.na(cov_df[[opt$depthCol]]))) {
    stop("Requested depth column in covFile cannot be safely converted to numeric: ", opt$depthCol)
  }
  names(depth) <- as.character(cov_df[[1]])
  cat("Depth filter source: covFile column '", opt$depthCol, "'.\n", sep = "")
}
cat(
  "Depth summary before filtering: min=", min(depth),
  ", median=", stats::median(depth),
  ", max=", max(depth), ".\n",
  sep = ""
)
if (opt$depth.filter > 0) {
  keep_depth <- depth > opt$depth.filter
  removed_n <- length(depth) - sum(keep_depth)
  cat("Applying depth.filter=", opt$depth.filter, ": keeping ", sum(keep_depth), " row(s), removing ", removed_n, ".\n", sep = "")
  if (!any(keep_depth)) {
    stop("No rows remain after applying --depth.filter=", opt$depth.filter)
  }
  abd_df <- abd_df[keep_depth, , drop = FALSE]
  cov_df <- cov_df[keep_depth, , drop = FALSE]
  depth <- depth[keep_depth]
  filtered <- filtered || removed_n > 0L
}

abd_ids <- as.character(abd_df[[1]])
cov_ids <- as.character(cov_df[[1]])
cat(
  "After row-level filtering: abd=", length(abd_ids),
  " row(s), cov=", length(cov_ids), " row(s).\n",
  sep = ""
)
if (!identical(abd_ids, cov_ids)) {
  stop("Internal error: abdFile and covFile row subject IDs are not matched after covariate processing.")
}
cat("Matched abd/cov row count: ", length(abd_ids), ".\n", sep = "")

# 2) Use genotype sample order as the reference order; genotype may contain extra samples
changed <- FALSE
geno_format <- PALMGWAS:::infer_geno_format(opt$genoFile)
geno_input <- read_geno_iid(
  geno_file = opt$genoFile
)
if (length(geno_input$cleanup) > 0L) {
  on.exit(unlink(geno_input$cleanup, force = TRUE), add = TRUE)
}
geno_ids <- geno_input$ids
cat("Genotype sample count: ", length(geno_ids), ".\n", sep = "")

cluster <- NULL
cluster_source <- NULL
if (!is.null(opt$clusterCol)) {
  cluster <- as.character(cov_df[[opt$clusterCol]])
  names(cluster) <- abd_ids
  cluster_source <- paste0("covFile column '", opt$clusterCol, "'")
  cat("Cluster filter source: ", cluster_source, ".\n", sep = "")
}

if (!is.null(cluster) && !abd_repeated_subjects) {
  cluster_values <- cluster[abd_ids]
  missing_cluster <- is.na(cluster_values) | !nzchar(trimws(cluster_values))
  removed_n <- sum(missing_cluster)
  cat(
    "Removing ", removed_n,
        " subject(s) with missing cluster values",
    if (!is.null(cluster_source)) paste0(" from ", cluster_source) else "",
    ".\n",
    sep = ""
  )
  if (removed_n > 0L) {
    keep_ids <- abd_ids[!missing_cluster]
    if (!length(keep_ids)) {
        stop("No subjects remain after removing subjects with missing cluster values.")
    }
    abd_df <- abd_df[as.character(abd_df[[1]]) %in% keep_ids, , drop = FALSE]
    cov_df <- cov_df[as.character(cov_df[[1]]) %in% keep_ids, , drop = FALSE]
    abd_ids <- as.character(abd_df[[1]])
    cov_ids <- as.character(cov_df[[1]])
    filtered <- TRUE
  }
} else {
  cat("No cluster missingness filtering requested, or repeated-subject clustering will use subject ID directly in Step2.1.\n")
}

abd_ids <- as.character(abd_df[[1]])
cov_ids <- as.character(cov_df[[1]])
subject_ids <- unique(abd_ids)
missing_in_geno <- setdiff(subject_ids, geno_ids)
if (length(missing_in_geno) > 0) {
  stop("Subject ID set mismatch: filtered abd/cov subjects missing from genotype input: ",
       paste(utils::head(missing_in_geno, 5), collapse = ", "))
}

ref_ids <- geno_ids[geno_ids %in% subject_ids]
cat(
  "Retained subjects present in genotype: ", length(ref_ids),
  ". Extra genotype-only subjects: ", length(geno_ids) - length(ref_ids), ".\n",
  sep = ""
)

row_ord <- order_rows_by_subject_ref(abd_ids, ref_ids)
if (!identical(row_ord, seq_len(nrow(abd_df)))) {
  cat("Order mismatch: reordering abdFile/covFile rows to match filtered genotype subject order...\n")
  abd_df <- abd_df[row_ord, , drop = FALSE]
  cov_df <- cov_df[row_ord, , drop = FALSE]
  changed <- TRUE
}

abd_ids2 <- as.character(abd_df[[1]])
cov_ids2 <- as.character(cov_df[[1]])

if (!identical(abd_ids2, cov_ids2)) stop("After reorder, abdFile/covFile subject ID order still mismatched.")
if (!identical(unique(abd_ids2), ref_ids)) stop("After reorder, abdFile subject order still mismatched to genotype.")

if (!is.null(opt$seqdepthInfoFile)) {
  cat("Generating DepthInfo from the final Step0 filtered/aligned sample set...\n")
  final_depth <- depth[row_ord]
  if (anyNA(final_depth)) {
    stop("Internal error: final retained row(s) missing depth values.")
  }
  seqdepth_df <- PALMGWAS:::seqdepth_info_from_values(abd_ids2, final_depth)
  dir.create(dirname(opt$seqdepthInfoFile), recursive = TRUE, showWarnings = FALSE)
  fwrite(seqdepth_df, file = opt$seqdepthInfoFile, sep = "\t", quote = FALSE, na = "NA", col.names = TRUE)
  cat(
    "DepthInfo finished: ", nrow(seqdepth_df),
    " final sample(s) written to ", opt$seqdepthInfoFile, ".\n",
    sep = ""
  )
} else {
  cat("DepthInfo skipped: --seqdepthInfoFile is NULL.\n")
}

dir.create(dirname(opt$abdAlignedFile), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$covAlignedFile), recursive = TRUE, showWarnings = FALSE)
abd_out <- abd_df
cov_out <- cov_df
if (!is.null(opt$timeIDCol)) {
  abd_out <- abd_out[, setdiff(colnames(abd_out), opt$timeIDCol), drop = FALSE]
  cov_out <- cov_out[, setdiff(colnames(cov_out), opt$timeIDCol), drop = FALSE]
}
fwrite(abd_out, file = opt$abdAlignedFile, sep = "\t", quote = FALSE, na = "NA", col.names = TRUE)
fwrite(cov_out, file = opt$covAlignedFile, sep = "\t", quote = FALSE, na = "NA", col.names = TRUE)
cat(
  "Final aligned row count: abd=", nrow(abd_out),
  ", cov=", nrow(cov_out), ".\n",
  sep = ""
)

if (changed || filtered) {
  cat("Done: wrote filtered/aligned abd file to ", opt$abdAlignedFile, "\n", sep = "")
  cat("Done: wrote filtered/aligned cov file to ", opt$covAlignedFile, "\n", sep = "")
} else {
  cat("ID check passed: wrote aligned copies without modifying original input files.\n")
}
