#!/bin/bash
#SBATCH --job-name=mgwas-visualization
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
echo "Start: Visualize results."
apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    Visualization.R \
    --metaDir=${metaDir} \
    --plotDir=${plotDir} \
    --pheno=${pheno} \
    --snp=${snp} \
    --showHet=${showHet:-TRUE}

apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    Visualization.R \
    --metaDir=${metaDir} \
    --plotDir=${plotDir} \
    --pheno=${pheno} \
    --showHet=${showHet:-TRUE}

apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    Visualization.R \
    --metaDir=${metaDir} \
    --plotDir=${plotDir} \
    --snp=${snp} \
    --showHet=${showHet:-TRUE}

apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    Visualization.R \
    --metaDir=${metaDir} \
    --plotDir=${plotDir} \
    --showHet=${showHet:-TRUE}
echo "Finish: Visualize results."
