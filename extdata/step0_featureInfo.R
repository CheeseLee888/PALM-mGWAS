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
    help = "Path to the output file for feature information"
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
abd_data <- read.table(opt$abdFile, header = TRUE, sep = "\t")

# Calculate prevalence and average proportion for each feature
# abd_data: data.frame
# first column is IID, remaining are gene counts

gene_mat <- as.matrix(abd_data[, -1])
gene_mat <- apply(gene_mat, 2, as.numeric)

n_sample <- nrow(gene_mat)

# -----------------------
# 1) Prevalence: None zero proportion
# -----------------------
prevalence <- colSums(gene_mat > 0, na.rm = TRUE) / n_sample

# -----------------------
# 2) AvgProportion
# normalize by row, then column mean
# -----------------------
row_sum <- rowSums(gene_mat, na.rm = TRUE)

# Avoid division by zero
row_sum[row_sum == 0] <- NA

gene_prop <- gene_mat / row_sum

avg_proportion <- colMeans(gene_prop, na.rm = TRUE)

# -----------------------
# 3) feature_info
# -----------------------
feature_info <- data.frame(
  FeatureID     = colnames(gene_mat),
  Prevalence    = prevalence,
  AvgProportion = avg_proportion,
  row.names = NULL
)

# Write feature information to the output file
if (is.null(opt$outputFile) || !nzchar(opt$outputFile)) {
  stop("The --outputFile argument is missing or empty. Please provide a valid file path.")
}
dir.create(dirname(opt$outputFile), recursive = TRUE, showWarnings = FALSE)
write.table(feature_info, file = opt$outputFile, sep = "\t", quote = FALSE, row.names = FALSE)


# Total samples
# total_samples <- nrow(fam_data)

# Append this information to the feature info file
# cat("Total samples: ", total_samples, "\n", file = opt$outputFile, append = TRUE)

cat("Feature information written to: ", opt$outputFile, "\n")