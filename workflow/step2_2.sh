#!/usr/bin/env bash
set -eo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################

# step2.2: compositional correction on completed step2 outputs
echo "Start: Perform step2.2 compositional correction."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step2_2_correct.R \
    --inputPrefix=${step2InputPrefix} \
    --outputPrefix=${step2OutputPrefix:-NULL} \
    --NULLmodelFile=${NULLmodelFile:-NULL} \
    --correct=${correct:-median}
echo "Finish: Perform step2.2 compositional correction."
