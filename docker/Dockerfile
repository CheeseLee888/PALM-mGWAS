FROM ubuntu:20.04

LABEL maintainer="cheese.lee.888@gmail.com"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV OMP_NUM_THREADS=1
ENV PIXI_NO_PATH_UPDATE=1

# ---- system deps ----
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
      pigz curl time wget ca-certificates unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

# ---- pixi env ----
RUN curl -fsSL https://pixi.sh/install.sh | sh && \
    mv /root/.pixi/bin/pixi /usr/local/bin/pixi && \
    pixi install && \
    rm -rf /root/.cache

# ---- PLINK 1.9 ----
RUN wget -qO /tmp/plink.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip && \
    unzip -q /tmp/plink.zip -d /usr/local/bin && \
    chmod +x /usr/local/bin/plink && \
    rm -f /tmp/plink.zip

# ---- deps ----
RUN pixi run Rscript -e 'install.packages("abess", repos="https://cloud.r-project.org")'

# ---- install PALM (only if PALMmGWAS depends on it) ----
# RUN pixi run Rscript - << 'EOF'
# tarball <- "thirdParty/PALM_0.1.0.tar.gz"
# if (!file.exists(tarball)) stop("PALM tarball not found: ", tarball)

# install.packages(tarball, repos = NULL, type = "source", dependencies = FALSE)
# library(PALM)
# message("PALM installed and loadable.")
# EOF

# Install PALM with explicit BLAS/LAPACK link flags (fix dgesvx_)
RUN pixi run env PKG_LIBS="-L/app/.pixi/envs/default/lib -lopenblas -llapack -lblas" \
    R CMD INSTALL thirdParty/PALM_0.1.0.tar.gz

# smoke test
RUN pixi run Rscript -e 'library(PALM); cat("PALM loaded OK\n")'

# ---- install PALMmGWAS ----
RUN pixi run R CMD INSTALL .

# ---- wrapper scripts (optional) ----
RUN for script in   step0_checkInput.R step1_null.R step2_summary.R step3_meta.R \
                    Info.R Visualization.R; do \
      echo '#!/bin/bash' > /usr/local/bin/$script && \
      echo 'cd /app' >> /usr/local/bin/$script && \
      echo 'exec pixi run Rscript /app/extdata/'$script' "$@"' >> /usr/local/bin/$script && \
      chmod +x /usr/local/bin/$script; \
    done

# optional: shrink image by removing package source
# RUN rm -rf R/ man/ src/ tests/ vignettes/ || true
