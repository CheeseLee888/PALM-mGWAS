#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################

# step3: meta-analysis
echo "Start: Performe meta analysis."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step3_meta.R \
    --studyDirFile=${studyDirFile} \
    --inputPrefix=${step3InputPrefix} \
    --chrom=${step3MetaChrom:-NULL} \
    --featureColList=${step3FeatureColList:-NULL} \
    --metaPrefix=${metaPrefix}
echo "Finish: Performe meta analysis."
