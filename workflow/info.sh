#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
mkdir -p "${outputFolder}"

echo "Start: Generate SNP information."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step0_snpInfo.R \
    --genoFile=${genoFile} \
    --outputFile=${outputSnpFile}
echo "Finish: Generate SNP information."

echo "Start: Generate feature information."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step0_featureInfo.R \
    --abdFile=${abdFile} \
    --outputFile=${outputFeatureFile}
echo "Finish: Generate feature information."

echo "Start: Generate sequencing depth information."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step0_seqdepthInfo.R \
    --abdFile=${abdFile} \
    --outputFile=${outputSeqDepthFile}
echo "Finish: Generate sequencing depth information."