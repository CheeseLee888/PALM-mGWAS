#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
echo "Start: Generate plots."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step4_plot.R \
    --inFile=${plotInFile} \
    --outdir=${plotOutDir}
echo "Finish: Generate plots."