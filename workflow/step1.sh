#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
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
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step1_palm.R \
    --abdFile=${abdFile} \
    --covFile=${covFile} \
    --outputPrefix=${palm1_step1_prefix}
echo "Finish: Fit null model."