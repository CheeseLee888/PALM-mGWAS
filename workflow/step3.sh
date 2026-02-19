#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################

# step3: meta-analysis
echo "Start: Performe meta analysis."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step3_meta.R \
    --studyDirList=${studyDirList} \
    --inPrefix=${inPrefix} \
    --outDir=${outDir}
echo "Finish: Performe meta analysis."