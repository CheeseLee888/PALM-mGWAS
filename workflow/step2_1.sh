#!/usr/bin/env bash
set -eo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
NULLObjPrefix="${NULLObjPrefix:-${NULLmodelFile:-}}"
NULLObjPrefix="${NULLObjPrefix%.rda}"
NULLmodelFile="${NULLmodelFile:-${NULLObjPrefix}.rda}"

# step2.1: score test without compositional correction
echo "Start: Perform step2.1 score test."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step2_1_summary.R \
    --genoFile=${genoFile} \
    --NULLmodelFile=${NULLmodelFile} \
    --PALMOutputFile=${palm1_step2_prefix} \
    --chrom=${chrom:-NULL} \
    --featureColList=${featureColList:-NULL} \
    --minMAF=${minMAF:-0.05} \
    --minMAC=${minMAC:-5} \
    --outputSnpFile=${outputSnpFile:-NULL} \
    --useCluster=${useCluster:-FALSE}
echo "Finish: Perform step2.1 score test."
