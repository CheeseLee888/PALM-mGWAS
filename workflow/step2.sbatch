#!/bin/bash
#SBATCH --job-name=palm-step2
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/scratch/group/ztang2/pli297/PALMmbQTL/logs/%x_%A_%a.out
#SBATCH --error=/mnt/scratch/group/ztang2/pli297/PALMmbQTL/logs/%x_%A_%a.err
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
if [[ "${PALMmethod}" != 1 && "${PALMmethod}" != 2 ]]; then
    echo "PALMmethod must be specified as 1 (PALM) or 2 (PALM-mbQTL)."
    exit 1
fi

if [[ "${PALMmethod}" == 1 ]]; then
    echo "Running PALM method..."
else
    echo "Running PALM-mbQTL method..."
fi

# step2: score test for phenoCol
echo "Start: Performe score test."
if [[ "${PALMmethod}" == 1 ]]; then
    apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
        step2_palm.R \
        --inFile=${genoFile} \
        --NULLmodelFile=${palm1_step1_prefix}.rda \
        --PALMOutputFile=${palm1_step2_prefix} \
        --chrom=${chrom}
else
    # Remove PALM2
fi
echo "Finish: Performe score test."