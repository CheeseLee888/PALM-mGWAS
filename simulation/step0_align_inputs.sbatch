#!/bin/bash
#SBATCH --job-name=simu_step0
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/simu_step0_%j.out
#SBATCH --error=logs/simu_step0_%j.err
#
# Usage:
#   sbatch --export=ALL,ENV_FILE=/path/to/your.env simulation/step0_align_inputs.sbatch
#
# Simulation root layout:
#   /mnt/scratch/group/ztang2/pli297/simulation/
#     ├── input/
#     ├── output/
#     ├── provided/
#     └── PALMmGWAS.sif
#
# Optional environment overrides:
#   ENV_FILE=/path/to/your.env
#   SIMU_ROOT=/path/to/simulation
#   SIF=/path/to/PALMmGWAS.sif
#   STUDY=study1
#   COVAR_COL_LIST=AGE,SEX,PC01,PC02,PC03,PC04,PC05
#   DEPTH_COL=NULL
#   DEPTH_FILTER=0
#   OUTPUT_SEQDEPTH_FILE=NULL

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
DEPTH_FILTER="${DEPTH_FILTER:-0}"
OUTPUT_SEQDEPTH_FILE="${OUTPUT_SEQDEPTH_FILE:-NULL}"

INPUT_DIR="${SIMU_ROOT}/input/${STUDY}"
OUTPUT_DIR="${SIMU_ROOT}/output/${STUDY}"

ABD_FILE="${INPUT_DIR}/abd.txt"
COV_FILE="${INPUT_DIR}/cov.txt"
GENO_PREFIX="${INPUT_DIR}/geno"
ABD_ALIGNED_FILE="${INPUT_DIR}/abd_aligned.txt"
COV_ALIGNED_FILE="${INPUT_DIR}/cov_aligned.txt"

mkdir -p "${SIMU_ROOT}/logs" "${OUTPUT_DIR}"

if [[ ! -f "${ABD_FILE}" || ! -f "${COV_FILE}" || ! -f "${GENO_PREFIX}.bed" || ! -f "${GENO_PREFIX}.bim" || ! -f "${GENO_PREFIX}.fam" ]]; then
  echo "Missing simulation input files. Run Rscript simulation/generate_simulation_data.R first." >&2
  exit 1
fi

echo "step0_align_inputs: SIMU_ROOT=${SIMU_ROOT}"
echo "step0_align_inputs: SIF=${SIF}"
echo "step0_align_inputs: STUDY=${STUDY}"
echo "step0_align_inputs: running Step0 alignment"
echo "step0_align_inputs: covarColList=${COVAR_COL_LIST}"
echo "step0_align_inputs: depthCol=${DEPTH_COL}, depth.filter=${DEPTH_FILTER}"
echo "step0_align_inputs: outputSeqDepthFile=${OUTPUT_SEQDEPTH_FILE}"

apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  step0_checkInput.R \
  --abdFile="${ABD_FILE}" \
  --covFile="${COV_FILE}" \
  --abdAlignedFile="${ABD_ALIGNED_FILE}" \
  --covAlignedFile="${COV_ALIGNED_FILE}" \
  --covarColList="${COVAR_COL_LIST}" \
  --depthCol="${DEPTH_COL}" \
  --depth.filter="${DEPTH_FILTER}" \
  --genoFile="${GENO_PREFIX}" \
  --outputSeqDepthFile="${OUTPUT_SEQDEPTH_FILE}"

echo "step0_align_inputs: finished"
echo "step0_align_inputs: aligned outputs"
echo "  ${ABD_ALIGNED_FILE}"
echo "  ${COV_ALIGNED_FILE}"
