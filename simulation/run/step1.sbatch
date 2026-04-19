#!/bin/bash
#SBATCH --job-name=simu_step1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/simu_step1_%A_%a.out
#SBATCH --error=logs/simu_step1_%A_%a.err
#
# Usage:
#   sbatch --export=ALL,ENV_FILE=/path/to/your.env simulation/run/step1.sbatch
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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_SIMU_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SIMU_ROOT="${SIMU_ROOT:-${DEFAULT_SIMU_ROOT}}"
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

echo "step1_fit_null_model: SIMU_ROOT=${SIMU_ROOT}"
echo "step1_fit_null_model: SIF=${SIF}"
echo "step1_fit_null_model: STUDY=${STUDY}"
echo "step1_fit_null_model: covarColList=${COVAR_COL_LIST}"
echo "step1_fit_null_model: depthCol=${DEPTH_COL}, prev.filter=${PREV_FILTER}"
echo "step1_fit_null_model: FeatureInfoFile=${OUTPUT_FEATURE_FILE}"

apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  step1_null.R \
  --abdFile="${ABD_ALIGNED_FILE}" \
  --covFile="${COV_ALIGNED_FILE}" \
  --covarColList="${COVAR_COL_LIST}" \
  --depthCol="${DEPTH_COL}" \
  --prev.filter="${PREV_FILTER}" \
  --FeatureInfoFile="${OUTPUT_FEATURE_FILE}" \
  --FeatureNameListFile="${FEATURE_LIST_FILE}" \
  --NULLObjPrefix="${NULLMODEL_FILE}"

echo "step1_fit_null_model: finished"
echo "step1_fit_null_model: outputs"
echo "  ${NULLMODEL_FILE}"
echo "  ${FEATURE_LIST_FILE}"
