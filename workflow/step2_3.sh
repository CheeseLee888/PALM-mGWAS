#!/usr/bin/env bash
set -eo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################

# step2.3: compositional correction on merged Step2 outputs
echo "Start: Perform step2.3 compositional correction."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step2_3_correct.R \
    --inputPrefix=${step2InputPrefix} \
    --outputPrefix=${step2OutputPrefix:-NULL} \
    --NULLmodelFile=${NULLmodelFile:-NULL}
echo "Finish: Perform step2.3 compositional correction."
