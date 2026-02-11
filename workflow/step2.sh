#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
if [[ "${PALMmethod}" != 1 && "${PALMmethod}" != 2 ]]; then
    echo "PALMmethod must be specified as 1 (PALM) or 2 (PALM-mbQTL)."
    exit 1
fi

if [[ "${PALMmethod}" == 1 ]]; then
    echo "Running PALM method..."
else
    echo "Running PALM-mbQTL method..."
fi

# step2: score test for phenoCol
echo "Start: Performe score test."
if [[ "${PALMmethod}" == 1 ]]; then
    pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step2_palm.R \
        --inFile=${genoFile} \
        --NULLmodelFile=${palm1_step1_prefix}.rda \
        --PALMOutputFile=${palm1_step2_prefix} \
        --chrom=${chrom}
else
    pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step2_scoreTest.R \
        --inFile=${genoFile} \
        --NULLmodelFile=${palm2_step1_prefix}.rda \
        --PALMOutputFile=${palm2_step2_prefix} \
        --chrom=${chrom}
fi
echo "Finish: Performe score test."