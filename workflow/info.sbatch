#!/bin/bash
#SBATCH --job-name=palm-info
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/scratch/group/ztang2/pli297/PALMmbQTL/logs/%x_%A_%a.out
#SBATCH --error=/mnt/scratch/group/ztang2/pli297/PALMmbQTL/logs/%x_%A_%a.err
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
mkdir -p "${outputFolder}"

echo "Start: Generate SNP/feature/seqdepth information."
apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    Info.R \
    --genoFile=${genoFile} \
    --abdFile=${abdFile} \
    --outputSnpFile=${outputSnpFile} \
    --outputFeatureFile=${outputFeatureFile} \
    --outputSeqDepthFile=${outputSeqDepthFile}
echo "Finish: Generate SNP/feature/seqdepth information."
