#!/usr/bin/env bash
set -eo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################

# step2.2: compositional correction on one Step2 scope
echo "Start: Perform step2.2 compositional correction."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step2_2_correction.R \
    --inputPrefix=${step2InputPrefix} \
    --chrom=${step2CorrectChrom:-NULL} \
    --overwriteOutput=${overwriteOutput:-TRUE}
echo "Finish: Perform step2.2 compositional correction."
