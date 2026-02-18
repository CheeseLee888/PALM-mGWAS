suppressPackageStartupMessages({
  library(PALM)
})
fitNULL <- function(abdFile,
                    covFile,
                    outputPrefix) {

  abd <- read_firstcol_as_rownames(abdFile)
  abd <- as.matrix(abd)
  
  if(is.null(covFile)) {
    # depth is calculated as the row sums of rel.abd
    cat("Fitting PALM null model without covariates.\n")
    modglmm <- PALM::palm.null.model(
      rel.abd = abd,
      prev.filter = 0
    )
  }else {
    cov <- read_firstcol_as_rownames(covFile)
    modglmm <- PALM::palm.null.model(
      rel.abd = abd,
      covariate.adjust = cov,
      prev.filter = 0
    )
  }
  
  save(modglmm, file = paste0(outputPrefix, ".rda"))
  cat("Done. PALM null model saved to", paste0(outputPrefix, ".rda"), "\n")
}