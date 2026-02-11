FROM ubuntu:20.04

LABEL maintainer="cheese.lee.888@gmail.com"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV OMP_NUM_THREADS=1
# Force step_2 to use 1 single thread. More threads are ineffective
ENV PIXI_NO_PATH_UPDATE=1

RUN apt-get update -qq && \
    apt-get install -y pigz && \
    apt-get -y install curl time && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
# Include files into docker image
COPY . .

RUN curl -fsSL https://pixi.sh/install.sh | sh && \
    mv /root/.pixi/bin/pixi /bin && pixi install && \
    rm -rf /root/.cache && \
    pixi run Rscript -e 'install.packages("lintools", repos="https://cloud.r-project.org")' && \
    pixi run Rscript -e 'install.packages("remotes", repos="https://cloud.r-project.org"); remotes::install_github("barkasn/fastSave")' && \
    pixi run Rscript -e 'install.packages(c("GMMAT"), repos="https://cloud.r-project.org")'


RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    wget ca-certificates unzip && \
    rm -rf /var/lib/apt/lists/*

# Install PLINK 1.9
RUN wget -qO /tmp/plink.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip && \
    unzip -q /tmp/plink.zip -d /usr/local/bin && \
    chmod +x /usr/local/bin/plink && \
    rm -f /tmp/plink.zip

# --- install PALM dependencies ---
RUN pixi run Rscript - << 'EOF'
deps <- c("MASS", "metafor", "abess", "dplyr", "brglm2", "Rcpp")

ip <- rownames(installed.packages())
missing <- setdiff(deps, ip)

if (length(missing)) {
  message("Installing PALM dependencies from CRAN: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

ip2 <- rownames(installed.packages())
still_missing <- setdiff(deps, ip2)
if (length(still_missing)) {
  stop("Failed to install PALM dependencies: ", paste(still_missing, collapse = ", "))
} else {
  message("All PALM dependencies installed: ", paste(deps, collapse = ", "))
}
EOF

# --- install PALM-mbQTL dependencies ---
RUN pixi run Rscript - << 'EOF'
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(
  c("gdsfmt", "SNPRelate", "SeqArray", "GenomeInfoDbData"),
  ask        = FALSE,
  update     = FALSE,
  dependencies = TRUE
)

ip <- rownames(installed.packages())
needed <- c("gdsfmt", "SNPRelate", "SeqArray", "GenomeInfoDbData")
missing <- setdiff(needed, ip)

if (length(missing)) {
  stop("Missing Bioconductor packages after install: ",
       paste(missing, collapse = ", "))
}
EOF

# --- install PALM package itself and perform explicit checks ---
RUN pixi run Rscript - << 'EOF'
tarball <- "thirdParty/PALM_0.1.0.tar.gz"
if (!file.exists(tarball)) {
  stop("PALM source tarball not found at: ", tarball,
       ". Please make sure it is included in the repo.")
}

message("Installing PALM from ", tarball, " ...")
install.packages(tarball,
                 repos = NULL,
                 type = "source",
                 dependencies = FALSE)

ip <- rownames(installed.packages())
if (!"PALM" %in% ip) {
  stop("PALM installation did not succeed: PALM not in installed.packages().")
}

message("Testing library(PALM) ...")
library(PALM)
message("PALM loaded successfully.")
EOF

# delete any precompiled files to ensure clean build
RUN rm -f src/*.so src/*.o || true

## 2. Install SAIGEQTL package (key)
RUN pixi run R CMD INSTALL .

# Check all key packages at once
RUN pixi run Rscript - << 'EOF'
ip <- rownames(installed.packages())

needed <- c(
  "lintools", "remotes", "fastSave", "GMMAT",
  "gdsfmt", "SNPRelate", "SeqArray", "GenomeInfoDbData",
  "SAIGEQTL",
  "MASS", "metafor", "abess", "dplyr", "brglm2", "Rcpp",
  "PALM"
)
missing <- setdiff(needed, ip)
if (length(missing)) {
  stop("Missing required R packages: ", paste(missing, collapse = ", "))
}

for (pkg in c("SAIGEQTL", "GMMAT", "SNPRelate", "SeqArray",
              "PALM", "abess")) {
  message("Testing library(", pkg, ") ...")
  library(pkg, character.only = TRUE)
}
message("All key packages installed and loadable.")
EOF

# Create Singularity-compatible wrapper scripts
RUN for script in step0_generateGRM.R step0_featureInfo.R step0_snpInfo.R step0_seqdepthInfo.R step0_checkInput.R step1_fitNULL.R step2_scoreTest.R step1_palm.R step2_palm.R step3_plot.R; do \
      echo '#!/bin/bash' > /usr/local/bin/$script && \
      echo 'cd /app' >> /usr/local/bin/$script && \
      echo 'exec pixi run Rscript /app/extdata/'$script' "$@"' >> /usr/local/bin/$script && \
      chmod +x /usr/local/bin/$script; \
    done

RUN rm -r src/ || true
RUN rm -r R/   || true
