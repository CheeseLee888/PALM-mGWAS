#!/bin/bash
#SBATCH --job-name=step0
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/step0_%A_%a.out
#SBATCH --error=logs/step0_%A_%a.err
# Step0 in the fixed simulation workflow:
# align abundance/covariate inputs for one study.

set -euo pipefail
WORK="${WORK:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIF="${SIF:-${WORK}/PALMmGWAS.sif}"
cd "${WORK}"

STUDY="${STUDY:-study1}"
COVAR_COL_LIST="${COVAR_COL_LIST:-AGE,SEX,PC01,PC02,PC03,PC04,PC05}"
DEPTH_COL="${DEPTH_COL:-NULL}"
DEPTH_FILTER="${DEPTH_FILTER:-0}"
SEQDEPTH_INFO_FILE="${SEQDEPTH_INFO_FILE:-NULL}"

INPUT_DIR="${WORK}/input/${STUDY}"
OUTPUT_DIR="${WORK}/output/${STUDY}"

ABD_FILE="${INPUT_DIR}/abd.txt"
COV_FILE="${INPUT_DIR}/cov.txt"
GENO_PREFIX="${INPUT_DIR}/geno"
ABD_ALIGNED_FILE="${INPUT_DIR}/abd_aligned.txt"
COV_ALIGNED_FILE="${INPUT_DIR}/cov_aligned.txt"

echo "step0_align_inputs: WORK=${WORK}"
echo "step0_align_inputs: SIF=${SIF}"
echo "step0_align_inputs: STUDY=${STUDY}"
echo "step0_align_inputs: running Step0 alignment"
echo "step0_align_inputs: covarColList=${COVAR_COL_LIST}"
echo "step0_align_inputs: depthCol=${DEPTH_COL}, depth.filter=${DEPTH_FILTER}"
echo "step0_align_inputs: SeqDepthInfoFile=${SEQDEPTH_INFO_FILE}"

apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}" \
  step0_checkInput.R \
  --abdFile="${ABD_FILE}" \
  --covFile="${COV_FILE}" \
  --abdAlignedFile="${ABD_ALIGNED_FILE}" \
  --covAlignedFile="${COV_ALIGNED_FILE}" \
  --covarColList="${COVAR_COL_LIST}" \
  --depthCol="${DEPTH_COL}" \
  --depth.filter="${DEPTH_FILTER}" \
  --genoFile="${GENO_PREFIX}" \
  --SeqDepthInfoFile="${SEQDEPTH_INFO_FILE}"

echo "step0_align_inputs: finished"
echo "step0_align_inputs: aligned outputs"
echo "  ${ABD_ALIGNED_FILE}"
echo "  ${COV_ALIGNED_FILE}"
