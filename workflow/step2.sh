#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################

# step2: score test
echo "Start: Performe score test."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step2_palm.R \
    --inFile=${genoFile} \
    --NULLmodelFile=${palm1_step1_prefix}.rda \
    --PALMOutputFile=${palm1_step2_prefix} \
    --chrom=${chrom} \
    --useOriginalFormat=${useOriginalFormat}
echo "Finish: Performe score test."