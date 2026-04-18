#' Merge per-chromosome Step2 result files into per-feature files
#'
#' @description
#' Reads chromosome-split Step2.1 result files and merges files sharing the
#' same feature into one all-chromosome file per feature. Repeated headers are
#' removed automatically.
#'
#' @param inputPrefix Base Step2 prefix. Under the current naming convention,
#'   an optional trailing underscore is ignored. Step2.1 outputs are expected at `<inputPrefix>_chr1_<feature>.txt` through
#'   `<inputPrefix>_chr22_<feature>.txt`, or already-merged files at
#'   `<inputPrefix>_allchr_<feature>.txt`.
#' @param outputPrefix Base output prefix for merged files. When `NULL`, the
#'   function derives the base prefix from `inputPrefix`, normalizes a trailing underscore when present, then writes merged
#'   files to `<outputPrefix>_allchr_<feature>.txt`.
#'
#' @return Invisibly returns the merged file paths.
#' @export
mergeSummary <- function(inputPrefix, outputPrefix = NULL) {
  if (missing(inputPrefix) || !nzchar(inputPrefix)) {
    stop("'inputPrefix' must be provided.")
  }

  normalize_step2_base <- function(prefix) {
    base <- basename(prefix)
    base <- sub("_+$", "", base)
    sub("^(.*)_(chr([1-9]|1[0-9]|2[0-2])|allchr)$", "\\1", base)
  }

  input_dir <- dirname(inputPrefix)
  base_prefix <- normalize_step2_base(inputPrefix)
  base_regex <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", base_prefix)
  chr_pattern <- paste0("^", base_regex, "_chr([1-9]|1[0-9]|2[0-2])_(.*)[.]txt$")

  if (is.null(outputPrefix) || !nzchar(outputPrefix)) {
    outputPrefix <- file.path(input_dir, base_prefix)
  }

  output_base_prefix <- normalize_step2_base(outputPrefix)
  output_dir <- dirname(outputPrefix)
  output_base <- paste0(output_base_prefix, "_allchr")

  files <- list.files(
    input_dir,
    pattern = chr_pattern,
    full.names = TRUE
  )

  extract_chrom <- function(path) {
    as.integer(sub(chr_pattern, "\\1", basename(path)))
  }
  extract_feature <- function(path) {
    sub(chr_pattern, "\\2", basename(path))
  }

  if (length(files)) {
    features <- vapply(files, extract_feature, character(1))
    chrom_num <- vapply(files, extract_chrom, integer(1))
    files <- files[order(features, chrom_num, files)]
  }

  if (!length(files)) {
    merged_pattern <- paste0(
      "^",
      gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", output_base),
      "_.*[.]txt$"
    )
    existing_merged <- list.files(
      output_dir,
      pattern = merged_pattern,
      full.names = TRUE
    )
    if (length(existing_merged)) {
      message("step2.2: merged files already exist under output prefix; leaving them in place.")
      return(invisible(unname(existing_merged[order(existing_merged)])))
    }
    stop(
      "No chromosome-split Step2 files found for base prefix: ", inputPrefix,
      ". Expected files named like ", base_prefix, "_chr1_<feature>.txt"
    )
  }

  features <- vapply(files, extract_feature, character(1))
  if (!length(features) || anyNA(features) || any(!nzchar(features))) {
    stop("Could not resolve feature names from Step2 files for inputPrefix: ", inputPrefix)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
