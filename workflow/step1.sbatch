#!/bin/bash
#SBATCH --job-name=mgwas-step1
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
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
apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    step1_null.R \
    --abdFile=${abdFile} \
    --covFile=${covFile} \
    --outputPrefix=${palm1_step1_prefix}
echo "Finish: Fit null model."