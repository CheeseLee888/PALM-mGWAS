#!/usr/bin/env bash
set -eo pipefail
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
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step1_null.R \
    --abdFile=${abdFile} \
    --phenoColList=${phenoColList:-NULL} \
    --covFile=${covFile} \
    --covarColList=${covarColList:-NULL} \
    --depthFile=${depthFile:-NULL} \
    --depth.filter=${depth_filter:-0} \
    --prev.filter=${prev_filter:-0.1} \
    --NULLmodelFile=${NULLmodelFile}
echo "Finish: Fit null model."
