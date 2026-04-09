# genotype input helpers shared by step0/info/step2

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

run_command_checked <- function(command, args) {
  message("Running command: ", command, " ", paste(args, collapse = " "))
  status <- system2(command, args = args)
  if (!identical(status, 0L)) {
    stop("Command failed with exit status ", status, ": ", command)
  }
}

open_geno_text_connection <- function(path) {
  lower <- tolower(path)
  if (grepl("\\.(gz|bgz)$", lower)) {
    return(gzfile(path, open = "rt"))
  }
  file(path, open = "rt")
}

parse_vcf_gt_dosage <- function(gt) {
  gt <- strsplit(gt, ":", fixed = TRUE)[[1]][1]
  if (!nzchar(gt) || grepl("\\.", gt, fixed = TRUE)) {
    return(NA_real_)
  }
  alleles <- strsplit(gt, "[/|]")[[1]]
  if (length(alleles) != 2L || any(!nzchar(alleles))) {
    return(NA_real_)
  }
  alleles_num <- suppressWarnings(as.numeric(alleles))
  if (any(is.na(alleles_num))) {
    return(NA_real_)
  }
  if (any(!alleles_num %in% c(0, 1))) {
    return(NA_real_)
  }
  2 - sum(alleles_num)
}

parse_vcf_sample_field <- function(sample_value, field_idx, vcfField) {
  parts <- strsplit(sample_value, ":", fixed = TRUE)[[1]]
  if (field_idx > length(parts)) {
    return(NA_real_)
  }
  value <- parts[field_idx]
  if (!nzchar(value) || value == ".") {
    return(NA_real_)
  }
  if (identical(vcfField, "DS")) {
    out <- suppressWarnings(as.numeric(value))
    if (is.na(out)) {
      return(NA_real_)
    }
    return(2 - out)
  }
  parse_vcf_gt_dosage(value)
}

read_vcf_header <- function(genoFile) {
  con <- open_geno_text_connection(genoFile)
  on.exit(close(con), add = TRUE)

  repeat {
    line <- readLines(con, n = 1L)
    if (!length(line)) {
      stop("VCF header line beginning with '#CHROM' not found: ", genoFile)
    }
    if (startsWith(line, "#CHROM")) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) < 9L) {
        stop("Invalid VCF header: expected at least 9 columns in ", genoFile)
      }
      return(list(
        fields = fields,
        sample_ids = if (length(fields) > 9L) fields[10:length(fields)] else character(0)
      ))
    }
  }
}

read_vcf_genotypes <- function(genoFile, vcfField = NULL) {
  if (is.null(vcfField) || !nzchar(trimws(as.character(vcfField)[1])) ||
      toupper(trimws(as.character(vcfField)[1])) == "NULL") {
    vcfField <- "AUTO"
  } else {
    vcfField <- toupper(trimws(as.character(vcfField)[1]))
  }
  if (!vcfField %in% c("AUTO", "DS", "GT")) {
    stop("'vcfField' must be 'AUTO', 'DS', or 'GT'.")
  }

  header <- read_vcf_header(genoFile)
  sample_ids <- header$sample_ids
  nsamples <- length(sample_ids)
  if (nsamples == 0L) {
    stop("VCF contains no sample columns: ", genoFile)
  }

  con <- open_geno_text_connection(genoFile)
  on.exit(close(con), add = TRUE)

  repeat {
    line <- readLines(con, n = 1L)
    if (!length(line)) {
      stop("VCF header line beginning with '#CHROM' not found: ", genoFile)
    }
    if (startsWith(line, "#CHROM")) {
      break
    }
  }

  geno_cols <- list()
  chr_vec <- character(0)
  snp_ids <- character(0)
  pos_vec <- integer(0)
  a1_vec <- character(0)
  a2_vec <- character(0)
  skipped_multiallelic <- 0L
  skipped_missing_format <- 0L
  fallback_to_gt <- 0L

  repeat {
    lines <- readLines(con, n = 1000L)
    if (!length(lines)) {
      break
    }
    for (line in lines) {
      if (!nzchar(line) || startsWith(line, "#")) {
        next
      }
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) < 9L + nsamples) {
        stop("Malformed VCF record with too few columns in ", genoFile)
      }
      alt <- fields[5]
      if (!nzchar(alt) || alt == "." || grepl(",", alt, fixed = TRUE)) {
        skipped_multiallelic <- skipped_multiallelic + 1L
        next
      }

      format_fields <- strsplit(fields[9], ":", fixed = TRUE)[[1]]
      field_idx <- NA_integer_
      field_name <- vcfField
      if (identical(vcfField, "AUTO")) {
        field_idx <- match("DS", format_fields)
        field_name <- "DS"
        if (is.na(field_idx)) {
          field_idx <- match("GT", format_fields)
          field_name <- "GT"
        }
      } else {
        field_idx <- match(vcfField, format_fields)
        if (is.na(field_idx) && identical(vcfField, "DS")) {
          field_idx <- match("GT", format_fields)
          if (!is.na(field_idx)) {
            field_name <- "GT"
            fallback_to_gt <- fallback_to_gt + 1L
          }
        }
      }
      if (is.na(field_idx)) {
        skipped_missing_format <- skipped_missing_format + 1L
        next
      }

      dosage <- vapply(
        fields[10:(9 + nsamples)],
        parse_vcf_sample_field,
        numeric(1),
        field_idx = field_idx,
        vcfField = field_name
      )

      snp_id <- fields[3]
      if (!nzchar(snp_id) || snp_id == ".") {
        snp_id <- paste(fields[1], fields[2], fields[4], alt, sep = ":")
      }

      geno_cols[[length(geno_cols) + 1L]] <- dosage
      chr_vec <- c(chr_vec, fields[1])
      snp_ids <- c(snp_ids, snp_id)
      pos_vec <- c(pos_vec, suppressWarnings(as.integer(fields[2])))
      a1_vec <- c(a1_vec, alt)
      a2_vec <- c(a2_vec, fields[4])
    }
  }

  if (!length(geno_cols)) {
    stop(
      "No usable VCF variants were loaded from ", genoFile,
      ". Skipped multiallelic=", skipped_multiallelic,
      ", skipped missing ", vcfField, " field=", skipped_missing_format, "."
    )
  }

  geno <- do.call(cbind, geno_cols)
  rownames(geno) <- sample_ids
  colnames(geno) <- make.unique(snp_ids, sep = "_dup")

  map <- data.frame(
    chromosome = chr_vec,
    snp.name = colnames(geno),
    position = pos_vec,
    allele.1 = a1_vec,
    allele.2 = a2_vec,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  message(
    "Loaded VCF genotype matrix directly: ", nrow(geno), " samples x ", ncol(geno), " SNPs",
    " (skipped multiallelic=", skipped_multiallelic,
    ", skipped missing ", if (identical(vcfField, "AUTO")) "DS/GT" else vcfField, " field=", skipped_missing_format,
    if (fallback_to_gt > 0L) paste0(", fallback DS->GT=", fallback_to_gt) else "",
    ")."
  )

  list(geno = geno, map = map, sample_ids = sample_ids)
}

prepare_plink_input <- function(genoFile,
                                vcfField = "DS",
                                alleleOrder = NULL,
                                plinkPath = "plink",
                                plink2Path = "plink2",
                                keepTemp = FALSE,
                                tempPrefix = NULL,
                                tempLabel = "plink") {
  if (missing(genoFile) || is.null(genoFile) || !nzchar(genoFile)) {
    stop("'genoFile' must be provided.")
  }
  if (!is.logical(keepTemp) || length(keepTemp) != 1L || is.na(keepTemp)) {
    stop("'keepTemp' must be TRUE or FALSE.")
  }

  genoFormat <- infer_geno_format(genoFile)
  message("Genotype input format inferred from genoFile: ", genoFormat)

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

  if (is.null(tempPrefix) || !nzchar(tempPrefix)) {
    tempPrefix <- file.path(
      tempdir(),
      paste0(
        sanitize_local_name(basename(genoFile)),
        "_",
        tempLabel,
        "_",
        genoFormat,
        "_plink"
      )
    )
  }
  cleanup <- paste0(tempPrefix, c(".bed", ".bim", ".fam", ".log", ".nosex"))

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
      args <- c(args, "--double-id", "--keep-allele-order", "--make-bed", "--out", tempPrefix)
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
      args <- c(args, "--double-id", "--keep-allele-order", "--make-bed", "--out", tempPrefix)
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
    args <- c(args, "--make-bed", "--out", tempPrefix)
    run_command_checked(plink2_exec, args)
  } else {
    stop("Unsupported inferred genotype format: ", genoFormat)
  }

  for (f in paste0(tempPrefix, c(".bed", ".bim", ".fam"))) {
    if (!file.exists(f)) {
      stop("Conversion to temporary PLINK files failed; missing output: ", f)
    }
  }

  if (keepTemp) {
    message("Keeping temporary converted PLINK files with prefix: ", tempPrefix)
  }

  list(prefix = tempPrefix, format = genoFormat, cleanup = cleanup)
}
