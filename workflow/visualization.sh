#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
echo "Start: Visualize results."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/Visualization.R \
    --metaDir=${metaDir:-${outputFolder}/meta} \
    --pattern=${pattern:-step3_meta_.*\\.txt$} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --pheno=${pheno:-NA} \
    --snp=${snp:-NA} \
    --pCut=${pCut:-1e-5} \
    --sigOnlyPheno=${sigOnlyPheno:-FALSE} \
    --showMeta=${showMeta:-TRUE} \
    --showHet=${showHet:-TRUE}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/Visualization.R \
    --metaDir=${metaDir:-${outputFolder}/meta} \
    --pattern=${pattern:-step3_meta_.*\\.txt$} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --pheno=${pheno:-NA} \
    --pCut=${pCut:-1e-5} \
    --sigOnlyPheno=${sigOnlyPheno:-FALSE} \
    --showMeta=${showMeta:-TRUE} \
    --showHet=${showHet:-TRUE}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/Visualization.R \
    --metaDir=${metaDir:-${outputFolder}/meta} \
    --pattern=${pattern:-step3_meta_.*\\.txt$} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --snp=${snp:-NA} \
    --pCut=${pCut:-1e-5} \
    --sigOnlyPheno=${sigOnlyPheno:-FALSE} \
    --showMeta=${showMeta:-TRUE} \
    --showHet=${showHet:-TRUE}

pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/Visualization.R \
    --metaDir=${metaDir:-${outputFolder}/meta} \
    --pattern=${pattern:-step3_meta_.*\\.txt$} \
    --plotDir=${plotDir:-${outputFolder}/plot} \
    --pCut=${pCut:-1e-5} \
    --printCut=${printCut:-1e-8} \
    --sigOnlyPheno=${sigOnlyPheno:-FALSE} \
    --showMeta=${showMeta:-TRUE} \
    --showHet=${showHet:-TRUE}
echo "Finish: Visualize results."
