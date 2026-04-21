#!/usr/bin/env bash
set -eo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
echo "Start: Step4 visualization."
step4MetaDir="${step4MetaDir:-${metaDir:-${outputFolder}/meta}}"

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step4_visualization.R \
    --metaDir=${step4MetaDir} \
    --plotPrefix=${plotPrefix:-${outputFolder}/plot/} \
    --feature=${feature:-NA} \
    --snp=${snp:-NA} \
    --showMeta=${showMeta:-TRUE} \
    --showHet=${showHet:-TRUE} \
    --width=${plotWidth:-NA} \
    --height=${plotHeight:-NA}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step4_visualization.R \
    --metaDir=${step4MetaDir} \
    --plotPrefix=${plotPrefix:-${outputFolder}/plot/} \
    --feature=${feature:-NA} \
    --plotMinP=${plotMinP:-NA} \
    --width=${plotWidth:-NA} \
    --height=${plotHeight:-NA}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step4_visualization.R \
    --metaDir=${step4MetaDir} \
    --plotPrefix=${plotPrefix:-${outputFolder}/plot/} \
    --snp=${snp:-NA} \
    --pCut=${pCut1:-5e-8} \
    --showMeta=${showMeta:-TRUE} \
    --showHet=${showHet:-TRUE} \
    --width=${plotWidth:-NA} \
    --height=${plotHeight:-NA}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step4_visualization.R \
    --metaDir=${step4MetaDir} \
    --plotPrefix=${plotPrefix:-${outputFolder}/plot/} \
    --pCut=${pCut2:-5e-8} \
    --plotMinP=${plotMinP:-NA} \
    --width=${plotWidth:-NA} \
    --height=${plotHeight:-NA}
echo "Finish: Step4 visualization."
