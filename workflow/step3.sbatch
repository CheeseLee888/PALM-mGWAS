#!/bin/bash
#SBATCH --job-name=mgwas-step3
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################

# step3: meta-analysis
echo "Start: Performe meta analysis."
apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    step3_meta.R \
    --studyDirList=${studyDirList} \
    --inPrefix=${inPrefix} \
    --outDir=${outDir}
echo "Finish: Performe meta analysis."