#!/usr/bin/env bash
set -eo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
echo "Start: Step4 visualization."
step4Pattern="${step4Pattern:-${visPattern:-step3_meta_.*\\.txt$}}"
step4MetaDir="${step4MetaDir:-${metaDir:-${outputFolder}/meta}}"

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step4_visualization.R \
    --metaDir=${step4MetaDir} \
    --pattern=${step4Pattern} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --pheno=${pheno:-NA} \
    --snp=${snp:-NA} \
    --showMeta=${showMeta:-TRUE} \
    --showHet=${showHet:-TRUE} \
    --width=${plotWidth:-NA} \
    --height=${plotHeight:-NA}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step4_visualization.R \
    --metaDir=${step4MetaDir} \
    --pattern=${step4Pattern} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --pheno=${pheno:-NA} \
    --plotMinP=${plotMinP:-NA} \
    --width=${plotWidth:-NA} \
    --height=${plotHeight:-NA}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step4_visualization.R \
    --metaDir=${step4MetaDir} \
    --pattern=${step4Pattern} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --snp=${snp:-NA} \
    --pCut=${pCut1:-5e-8} \
    --showMeta=${showMeta:-TRUE} \
    --showHet=${showHet:-TRUE} \
    --width=${plotWidth:-NA} \
    --height=${plotHeight:-NA}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step4_visualization.R \
    --metaDir=${step4MetaDir} \
    --pattern=${step4Pattern} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --pCut=${pCut2:-5e-8} \
    --plotMinP=${plotMinP:-NA} \
    --width=${plotWidth:-NA} \
    --height=${plotHeight:-NA}
echo "Finish: Step4 visualization."
