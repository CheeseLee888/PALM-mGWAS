#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
echo "Start: Visualize results."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step_visualization.R \
    --metaDir=${metaDir} \
    --plotDir=${plotDir} \
    --pheno=${pheno} \
    --snp=${snp}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step_visualization.R \
    --metaDir=${metaDir} \
    --plotDir=${plotDir} \
    --pheno=${pheno}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step_visualization.R \
    --metaDir=${metaDir} \
    --plotDir=${plotDir} \
    --snp=${snp}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step_visualization.R \
    --metaDir=${metaDir} \
    --plotDir=${plotDir}
echo "Finish: Visualize results."