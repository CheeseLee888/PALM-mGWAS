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

mkdir -p "${outputFolder}"

# step0: check input files
echo "Start: Check input files."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step0_checkInput.R \
    --abdFile=${abdFile} \
    --covFile=${covFile} \
    --genoFile=${genoFile}
echo "Finish: Check input files."    

# step1: fit null model for all phenotypes
echo "Start: Fit null model."
if [[ "${PALMmethod}" == 1 ]]; then
    pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step1_palm.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --outputPrefix=${palm1_step1_prefix}
else
    # Remove PALM2
fi