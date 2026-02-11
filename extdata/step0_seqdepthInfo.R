#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--abdFile",
    type = "character", default = "",
    help = "Path to the microbiome abundance data file"
  ),
  make_option("--outputFile",
    type = "character", default = "",
    help = "Path to the output file for sample information"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate abdFile before reading
if (is.null(opt$abdFile) || !nzchar(opt$abdFile)) {
  stop("The --abdFile argument is missing or empty. Please provide a valid file path.")
}
if (!file.exists(opt$abdFile)) {
  stop("The specified abdFile does not exist: ", opt$abdFile)
}

# Summarize microbiome abundance data
abd <- read.table(opt$abdFile, header = TRUE, sep = "\t")
sample_id <- abd[[1]]
abd_mat <- as.matrix(abd[, -1])
rownames(abd_mat) <- sample_id

# Calculate sequencing depth for each sample
SeqDepth <- rowSums(abd_mat, na.rm = TRUE)

# Create a data frame for sample information
seqdepth_info <- data.frame(
  SampleID = sample_id,
  SeqDepth = SeqDepth,
  row.names = NULL
)

# Write sample information to the output file
if (is.null(opt$outputFile) || !nzchar(opt$outputFile)) {
  stop("The --outputFile argument is missing or empty. Please provide a valid file path.")
}
write.table(seqdepth_info, file = opt$outputFile, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Sequencing depth information written to: ", opt$outputFile, "\n")