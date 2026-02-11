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

echo "Start: Generate SNP information."
apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    step0_snpInfo.R \
    --genoFile=${genoFile} \
    --outputFile=${outputSnpFile}
echo "Finish: Generate SNP information."

echo "Start: Generate feature information."
apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    step0_featureInfo.R \
    --abdFile=${abdFile} \
    --outputFile=${outputFeatureFile}
echo "Finish: Generate feature information."

echo "Start: Generate sequencing depth information."
apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    step0_seqdepthInfo.R \
    --abdFile=${abdFile} \
    --outputFile=${outputSeqDepthFile}
echo "Finish: Generate sequencing depth information."