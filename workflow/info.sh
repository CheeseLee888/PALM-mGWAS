#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
mkdir -p "${outputFolder}"

echo "Start: Generate SNP/feature/seqdepth information."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/Info.R \
    --genoPrefix=${genoPrefix} \
    --abdFile=${abdFile} \
    --outputSnpFile=${outputSnpFile} \
    --outputFeatureFile=${outputFeatureFile} \
    --outputSeqDepthFile=${outputSeqDepthFile}
echo "Finish: Generate SNP/feature/seqdepth information."
