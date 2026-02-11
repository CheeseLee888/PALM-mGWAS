#!/bin/bash
#SBATCH --job-name=palm-step1
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

mkdir -p "${outputFolder}"

# step0: check input files
echo "Start: Check input files."
apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    step0_checkInput.R \
    --abdFile=${abdFile} \
    --covFile=${covFile} \
    --genoFile=${genoFile}
echo "Finish: Check input files."    

# step1: fit null model for all phenotypes
echo "Start: Fit null model."
if [[ "${PALMmethod}" == 1 ]]; then
    apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
        step1_palm.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --outputPrefix=${palm1_step1_prefix}
else
    # Remove PALM2
fi
echo "Finish: Fit null model."