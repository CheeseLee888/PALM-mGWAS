#!/bin/bash
#SBATCH --job-name=step1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/step1_%A_%a.out
#SBATCH --error=logs/step1_%A_%a.err
# Step1 in the fixed simulation workflow:
# fit the null model and write feature info for one study.

set -euo pipefail
WORK="${WORK:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIF="${SIF:-${WORK}/PALMGWAS.sif}"
cd "${WORK}"

STUDY="${STUDY:-study1}"
COV_LIST="${COV_LIST:-AGE,SEX,PC01,PC02,PC03,PC04,PC05}"
DEPTH_COL="${DEPTH_COL:-NULL}"
PREV_FILTER="${PREV_FILTER:-0.1}"
OUTPUT_FEATURE_FILE="${OUTPUT_FEATURE_FILE:-NULL}"

INPUT_DIR="${WORK}/input/${STUDY}"
OUTPUT_DIR="${WORK}/output/${STUDY}"

ABD_ALIGNED_FILE="${INPUT_DIR}/abd_aligned.txt"
COV_ALIGNED_FILE="${INPUT_DIR}/cov_aligned.txt"
nullModelPrefix="${OUTPUT_DIR}/step1_allpheno"

echo "step1_fit_null_model: WORK=${WORK}"
echo "step1_fit_null_model: SIF=${SIF}"
echo "step1_fit_null_model: STUDY=${STUDY}"
echo "step1_fit_null_model: covarColList=${COV_LIST}"
echo "step1_fit_null_model: depthCol=${DEPTH_COL}, prev.filter=${PREV_FILTER}"
echo "step1_fit_null_model: featureInfoFile=${OUTPUT_FEATURE_FILE}"

apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}" \
  step1_null.R \
  --abdFile="${ABD_ALIGNED_FILE}" \
  --covFile="${COV_ALIGNED_FILE}" \
  --covarColList="${COV_LIST}" \
  --depthCol="${DEPTH_COL}" \
  --prev.filter="${PREV_FILTER}" \
  --featureInfoFile="${OUTPUT_FEATURE_FILE}" \
  --nullModelPrefix="${nullModelPrefix}"

echo "step1_fit_null_model: finished"
echo "step1_fit_null_model: outputs"
echo "  ${nullModelPrefix}.rda"
echo "  ${OUTPUT_FEATURE_FILE}"
