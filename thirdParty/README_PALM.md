## PALM

`PALM` is a R package and it implements a quasi-Poisson-based framework designed for robust, scalable, and reproducible identification of covariate-associated microbial features in large-scale microbiome association studies and meta-analyses. The package can perform single-study association analysis and meta-analysis across multiple studies.
See following items for more details:

* [`PALM` Manual](https://github.com/ZjpWei/PALM_package/blob/main/doc/PALM_0.1.0.pdf).

* [`PALM` Vignette](https://ZjpWei.github.io/PALM_package/PALM_vignette.html).

* Article: Wei et al. *Fast and reliable association discovery in large-scale microbiome studies and meta-analyses using PALM*. Submitted.

## Author

Zhoujingpeng Wei @[Tang Lab](https://tangzheng1.github.io/tanglab/)

Department of Biostatistics and Medical Informatics, University of Wisconsin-Madison

## System requirements

The `PALM` package (version 0.1.0) should be compatible with Windows, Mac, and Linux operating systems.

Before setting up the package, users should have R version 4.3.0 or higher.

The package depends on the following R packages: `Mass`, `metafor`, `abess`, `dplyr`, `brglm2` and `Rcpp`.

## Installation guide

Install package from github.
```{r}
devtools::install_github("ZjpWei/PALM_package")
```

Install package by [source code](https://github.com/ZjpWei/PALM_package/blob/main/PALM_0.1.0.tar.gz)
```{r}
install.packages("./PALM_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Quick start guide

The following are minimal examples on functionalities of `PALM`. For more detailes, please refer to its [vignette](https://htmlpreview.github.io/?https://github.com/ZjpWei/PALM_package/blob/main/doc/PALM_vignette.html).

* Load package and data
```{r}
library("PALM")
data("CRC_data", package = "PALM")
CRC_abd <- CRC_data$CRC_abd
CRC_meta <- CRC_data$CRC_meta

########## Generate summary statistics ##########
rel.abd <- list()
covariate.interest <- list()
for(d in unique(CRC_meta$Study)){
  rel.abd[[d]] <- CRC_abd[CRC_meta$Sample_ID[CRC_meta$Study == d],]
  disease <- as.numeric(CRC_meta$Group[CRC_meta$Study == d] == "CRC")
  names(disease) <- CRC_meta$Sample_ID[CRC_meta$Study == d]
  covariate.interest[[d]] <- matrix(disease, ncol = 1, dimnames = list(names(disease), "disease"))
}
```

* Perform meta-analysis
```{r}
meta.result <- palm(rel.abd = rel.abd, covariate.interest = covariate.interest)
```

## Issues tracker

Please use the [issues tracker](https://github.com/ZjpWei/PALM_package/issues) to report any bugs or give any suggestions. Please feel free to contact me [Weizhou](mailto:zwei74@wisc.edu?subject=[GitHub]%20Melody%20paper%20analysis).

## License

This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3
