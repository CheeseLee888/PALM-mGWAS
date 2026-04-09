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
    --abdAlignedFile=${abdAlignedFile} \
    --covAlignedFile=${covAlignedFile} \
    --covarColList=${covarColList:-NULL} \
    --depthCol=${depthCol:-NULL} \
    --depth.filter=${depth_filter:-0} \
    --genoFile=${genoFile} \
    --vcfField=${vcfField:-DS} \
    --alleleOrder=${alleleOrder:-NULL} \
    --keepTemp=${keepTemp:-FALSE} \
    --outputSeqDepthFile=${outputSeqDepthFile:-NULL}
echo "Finish: Check input files."    

# step1: fit null model for all phenotypes
echo "Start: Fit null model."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step1_null.R \
    --abdFile=${abdAlignedFile} \
    --covFile=${covAlignedFile} \
    --covarColList=${covarColList:-NULL} \
    --depthCol=${depthCol:-NULL} \
    --prev.filter=${prev_filter:-0.1} \
    --outputFeatureFile=${outputFeatureFile:-NULL} \
    --NULLmodelFile=${NULLmodelFile}
echo "Finish: Fit null model."
