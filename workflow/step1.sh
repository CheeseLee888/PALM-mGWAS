#!/usr/bin/env bash
set -eo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
mkdir -p "${outputFolder}"

# Prefer renamed Step1 variables while remaining compatible with existing config.sh files.
FeatureInfoFile="${FeatureInfoFile:-${outputFeatureFile:-NULL}}"
FeatureNameListFile="${FeatureNameListFile:-${outputFeatureListFile:-NULL}}"
NULLObjPrefix="${NULLObjPrefix:-${NULLmodelFile:-}}"
SeqDepthInfoFile="${SeqDepthInfoFile:-NULL}"
NULLObjPrefix="${NULLObjPrefix%.rda}"

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
    --SeqDepthInfoFile=${SeqDepthInfoFile:-NULL}
echo "Finish: Check input files."    

# step1: fit null model for all phenotypes
echo "Start: Fit null model."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step1_null.R \
    --abdFile=${abdAlignedFile} \
    --covFile=${covAlignedFile} \
    --covarColList=${covarColList:-NULL} \
    --depthCol=${depthCol:-NULL} \
    --prev.filter=${prev_filter:-0.1} \
    --FeatureInfoFile=${FeatureInfoFile:-NULL} \
    --FeatureNameListFile=${FeatureNameListFile:-NULL} \
    --NULLObjPrefix=${NULLObjPrefix}
echo "Finish: Fit null model."
