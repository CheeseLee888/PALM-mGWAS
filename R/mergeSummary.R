#' Merge per-chromosome Step2 result files into per-feature files
#'
#' @description
#' Reads Step2 result files produced from chromosome-split Step2.1 runs and
#' merges files sharing the same feature into one all-chromosome file per
#' feature. Repeated headers are removed automatically.
#'
#' @param inputPrefix Step2 input prefix. Files are expected at either
#'   `<inputPrefix>_<feature>.txt` or, for chromosome-split runs,
#'   `<inputPrefix-with-chr-pattern>_<feature>.txt` such as
#'   `step2_chr1_g_Bifidobacterium.txt`.
#' @param outputPrefix Output prefix for merged files. When `NULL`, the function
#'   derives an all-chromosome prefix from `inputPrefix`.
#'
#' @return Invisibly returns the merged file paths.
#' @export
mergeSummary <- function(inputPrefix, outputPrefix = NULL) {
  if (missing(inputPrefix) || !nzchar(inputPrefix)) {
    stop("'inputPrefix' must be provided.")
  }

  if (is.null(outputPrefix) || !nzchar(outputPrefix)) {
    input_dir <- dirname(inputPrefix)
    input_base <- basename(inputPrefix)
    if (grepl("^(.*_chr)[^_]+$", input_base)) {
      output_base <- sub("^(.*_chr)[^_]+$", "\\1allchr", input_base)
    } else {
      output_base <- paste0(input_base, "_allchr")
    }
    outputPrefix <- file.path(input_dir, output_base)
  }

  input_dir <- dirname(inputPrefix)
  input_base <- basename(inputPrefix)
  output_dir <- dirname(outputPrefix)
  output_base <- basename(outputPrefix)

  if (grepl("^(.*_chr)[^_]+$", input_base)) {
    find_base <- sub("^(.*_chr)[^_]+$", "\\1", input_base)
    pattern <- paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", find_base), "[^_]+_(.*)[.]txt$")
    extract_chrom <- function(path) {
      sub(
        paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", find_base), "([^_]+)_.*[.]txt$"),
        "\\1",
        basename(path)
      )
    }
    extract_feature <- function(path) {
      sub(pattern, "\\1", basename(path))
    }
  } else {
    find_base <- input_base
    pattern <- paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", find_base), "_(.*)[.]txt$")
    extract_chrom <- function(path) {
      NA_character_
    }
    extract_feature <- function(path) {
      sub(pattern, "\\1", basename(path))
    }
  }

  files <- list.files(
    input_dir,
    pattern = pattern,
    full.names = TRUE
  )
  features <- vapply(files, extract_feature, character(1))
  chroms <- vapply(files, extract_chrom, character(1))
  chrom_is_num <- grepl("^[0-9]+$", chroms)
  chrom_num <- suppressWarnings(as.integer(chroms))
  chrom_rank <- ifelse(chrom_is_num, chrom_num, .Machine$integer.max)
  files <- files[order(features, chrom_rank, chroms, files)]

  if (!length(files)) {
    existing_merged <- list.files(
      output_dir,
      pattern = paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", output_base), "_.*[.]txt$"),
      full.names = TRUE
    )
    if (length(existing_merged)) {
      message("step2.2: merged files already exist under output prefix; leaving them in place.")
      return(invisible(unname(existing_merged[order(existing_merged)])))
    }
    stop("No Step2 files found for inputPrefix: ", inputPrefix)
  }

  features <- vapply(files, extract_feature, character(1))
  if (!length(features) || anyNA(features) || any(!nzchar(features))) {
    stop("Could not resolve feature names from Step2 files for inputPrefix: ", inputPrefix)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (identical(normalizePath(outputPrefix, winslash = "/", mustWork = FALSE),
                normalizePath(inputPrefix, winslash = "/", mustWork = FALSE)) &&
      startsWith(output_base, find_base) &&
      all(startsWith(basename(files), paste0(output_base, "_")))) {
    message("step2.2: input already appears merged; leaving files in place.")
    out_paths <- file.path(output_dir, paste0(output_base, "_", unique(features), ".txt"))
    return(invisible(unname(out_paths)))
  }

  split_files <- split(files, features)
  out_paths <- character(length(split_files))
  names(out_paths) <- names(split_files)

  for (feature in names(split_files)) {
    feature_files <- split_files[[feature]]
    out_file <- file.path(output_dir, paste0(output_base, "_", feature, ".txt"))

    merged <- lapply(seq_along(feature_files), function(i) {
      lines <- readLines(feature_files[[i]], warn = FALSE)
      if (i > 1L && length(lines) > 0L) {
        lines <- lines[-1L]
      }
      lines
    })
    writeLines(unlist(merged, use.names = FALSE), out_file)

    message(
      "step2.2: merged ", length(feature_files),
      " chromosome file(s) -> ", out_file
    )
    out_paths[[feature]] <- out_file
  }

  invisible(unname(out_paths))
}
