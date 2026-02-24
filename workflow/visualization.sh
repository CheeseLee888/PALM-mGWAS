#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
echo "Start: Visualize results."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step_visualization.R \
    --metaDir=${metaDir:-${outputFolder}/meta} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --pheno=${pheno:-NA} \
    --snp=${snp:-NA} \
    --pCut=${pCut:-1e-5} \
    --sigOnlyPheno=${sigOnlyPheno:-FALSE} \
    --ciMult=${ciMult:-1.96} \
    --showMeta=${showMeta:-TRUE} \
    --xlim=${xlim:-NA} \
    --width=${plotWidth:-10} \
    --height=${plotHeight:-6} \
    --dpi=${plotDpi:-300}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step_visualization.R \
    --metaDir=${metaDir:-${outputFolder}/meta} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --pheno=${pheno:-NA} \
    --pCut=${pCut:-1e-5} \
    --sigOnlyPheno=${sigOnlyPheno:-FALSE} \
    --ciMult=${ciMult:-1.96} \
    --showMeta=${showMeta:-TRUE} \
    --xlim=${xlim:-NA} \
    --width=${plotWidth:-10} \
    --height=${plotHeight:-6} \
    --dpi=${plotDpi:-300}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step_visualization.R \
    --metaDir=${metaDir:-${outputFolder}/meta} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --snp=${snp:-NA} \
    --pCut=${pCut:-1e-5} \
    --sigOnlyPheno=${sigOnlyPheno:-FALSE} \
    --ciMult=${ciMult:-1.96} \
    --showMeta=${showMeta:-TRUE} \
    --xlim=${xlim:-NA} \
    --width=${plotWidth:-10} \
    --height=${plotHeight:-6} \
    --dpi=${plotDpi:-300}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step_visualization.R \
    --metaDir=${metaDir:-${outputFolder}/meta} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --pCut=${pCut:-1e-5} \
    --printCut=${printCut:-1e-8} \
    --sigOnlyPheno=${sigOnlyPheno:-FALSE} \
    --ciMult=${ciMult:-1.96} \
    --showMeta=${showMeta:-TRUE} \
    --xlim=${xlim:-NA} \
    --width=${plotWidth:-10} \
    --height=${plotHeight:-6} \
    --dpi=${plotDpi:-300}
echo "Finish: Visualize results."
