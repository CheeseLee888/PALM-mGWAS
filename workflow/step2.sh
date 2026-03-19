#!/usr/bin/env bash
set -eo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################

# step2: score test
echo "Start: Performe score test."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step2_summary.R \
    --genoFile=${genoFile} \
    --vcfField=${vcfField:-DS} \
    --alleleOrder=${alleleOrder:-NULL} \
    --keepTemp=${keepTemp:-FALSE} \
    --NULLmodelFile=${NULLmodelFile} \
    --PALMOutputFile=${palm1_step2_prefix} \
    --chrom=${chrom:-NULL} \
    --minMAF=${minMAF:-0.1} \
    --minMAC=${minMAC:-5} \
    --correct=${correct:-NULL} \
    --useCluster=${useCluster:-FALSE}
echo "Finish: Performe score test."
