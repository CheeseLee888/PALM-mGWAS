#!/usr/bin/env bash
set -eo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################

# step2: score test
echo "Start: Perform score test."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step2_summary.R \
    --genoFile=${genoFile} \
    --NULLmodelFile=${NULLmodelFile} \
    --PALMOutputFile=${palm1_step2_prefix} \
    --chrom=${chrom:-NULL} \
    --minMAF=${minMAF:-0.05} \
    --minMAC=${minMAC:-5} \
    --outputSnpFile=${outputSnpFile:-NULL} \
    --correct=${correct:-NULL} \
    --useCluster=${useCluster:-FALSE}
echo "Finish: Perform score test."
