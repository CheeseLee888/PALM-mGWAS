#!/bin/bash
#SBATCH --job-name=simu_step1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/simu_step1_%A_%a.out
#SBATCH --error=logs/simu_step1_%A_%a.err
#
# Usage:
#   sbatch --export=ALL,ENV_FILE=/path/to/your.env simulation/step1_fit_null_model.sbatch
#
# Simulation root layout:
#   /mnt/scratch/group/ztang2/pli297/simulation/
#     ├── input/
#     ├── output/
#     └── PALMmGWAS.sif
#
# Optional environment overrides:
#   ENV_FILE=/path/to/your.env
#   SIMU_ROOT=/path/to/simulation
#   SIF=/path/to/PALMmGWAS.sif
#   STUDY=study1
#   COVAR_COL_LIST=AGE,SEX,PC01,PC02,PC03,PC04,PC05
#   DEPTH_COL=NULL
#   PREV_FILTER=0.1
#   OUTPUT_FEATURE_FILE=NULL

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi
DEFAULT_SIMU_ROOT="/mnt/scratch/group/ztang2/pli297/simulation"
SIMU_ROOT="${SIMU_ROOT:-${PWD:-${DEFAULT_SIMU_ROOT}}}"
SIF="${SIF:-${SIMU_ROOT}/PALMmGWAS.sif}"
cd "${SIMU_ROOT}"

STUDY="${STUDY:-study1}"
COVAR_COL_LIST="${COVAR_COL_LIST:-AGE,SEX,PC01,PC02,PC03,PC04,PC05}"
DEPTH_COL="${DEPTH_COL:-NULL}"
PREV_FILTER="${PREV_FILTER:-0.1}"
OUTPUT_FEATURE_FILE="${OUTPUT_FEATURE_FILE:-NULL}"

INPUT_DIR="${SIMU_ROOT}/input/${STUDY}"
OUTPUT_DIR="${SIMU_ROOT}/output/${STUDY}"

ABD_ALIGNED_FILE="${INPUT_DIR}/abd_aligned.txt"
COV_ALIGNED_FILE="${INPUT_DIR}/cov_aligned.txt"
NULLMODEL_FILE="${OUTPUT_DIR}/step1_allpheno.rda"
FEATURE_LIST_FILE="${OUTPUT_DIR}/feature_list.txt"

mkdir -p "${SIMU_ROOT}/logs" "${OUTPUT_DIR}"

if [[ ! -f "${ABD_ALIGNED_FILE}" || ! -f "${COV_ALIGNED_FILE}" ]]; then
  echo "Missing Step0-aligned input files. Run simulation/step0_align_inputs.sbatch first." >&2
  exit 1
fi

echo "step1_fit_null_model: SIMU_ROOT=${SIMU_ROOT}"
echo "step1_fit_null_model: SIF=${SIF}"
echo "step1_fit_null_model: STUDY=${STUDY}"
echo "step1_fit_null_model: covarColList=${COVAR_COL_LIST}"
echo "step1_fit_null_model: depthCol=${DEPTH_COL}, prev.filter=${PREV_FILTER}"
echo "step1_fit_null_model: outputFeatureFile=${OUTPUT_FEATURE_FILE}"

apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  step1_null.R \
  --abdFile="${ABD_ALIGNED_FILE}" \
  --covFile="${COV_ALIGNED_FILE}" \
  --covarColList="${COVAR_COL_LIST}" \
  --depthCol="${DEPTH_COL}" \
  --prev.filter="${PREV_FILTER}" \
  --outputFeatureFile="${OUTPUT_FEATURE_FILE}" \
  --outputFeatureListFile="${FEATURE_LIST_FILE}" \
  --NULLmodelFile="${NULLMODEL_FILE}"

echo "step1_fit_null_model: finished"
echo "step1_fit_null_model: outputs"
echo "  ${NULLMODEL_FILE}"
echo "  ${FEATURE_LIST_FILE}"
