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
