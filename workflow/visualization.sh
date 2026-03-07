#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
echo "Start: Visualize results."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/Visualization.R \
    --metaDir=${metaDir:-${outputFolder}/meta} \
    --pattern=${visPattern:-step3_meta_.*\\.txt$} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --pheno=${pheno:-NA} \
    --snp=${snp:-NA} \
    --showMeta=${showMeta:-TRUE} \
    --showHet=${showHet:-TRUE} \
    --width=${plotWidth:-NA} \
    --height=${plotHeight:-NA}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/Visualization.R \
    --metaDir=${metaDir:-${outputFolder}/meta} \
    --pattern=${visPattern:-step3_meta_.*\\.txt$} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --pheno=${pheno:-NA} \
    --width=${plotWidth:-NA} \
    --height=${plotHeight:-NA}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/Visualization.R \
    --metaDir=${metaDir:-${outputFolder}/meta} \
    --pattern=${visPattern:-step3_meta_.*\\.txt$} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --snp=${snp:-NA} \
    --pCut=${pCut:-NA} \
    --showMeta=${showMeta:-TRUE} \
    --showHet=${showHet:-TRUE} \
    --width=${plotWidth:-NA} \
    --height=${plotHeight:-NA}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/Visualization.R \
    --metaDir=${metaDir:-${outputFolder}/meta} \
    --pattern=${visPattern:-step3_meta_.*\\.txt$} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --printCut=${printCut:-1e-8} \
    --width=${plotWidth:-NA} \
    --height=${plotHeight:-NA}
echo "Finish: Visualize results."
