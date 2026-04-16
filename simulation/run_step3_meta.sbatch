#!/bin/bash
#SBATCH --job-name=simu_step3_meta
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/simu_step3_meta_%A_%a.out
#SBATCH --error=logs/simu_step3_meta_%A_%a.err

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
META_DIR="${META_DIR:-/mnt/scratch/group/ztang2/pli297/simulation/output/meta}"
FEATURE_LIST_FILE="${FEATURE_LIST_FILE:-${SIMU_ROOT}/output/study1/feature_list.txt}"
cd "${SIMU_ROOT}"
mkdir -p "${SIMU_ROOT}/logs"

echo "run_step3_meta.sbatch: simulation root = ${SIMU_ROOT}"
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  if [[ ! -f "${FEATURE_LIST_FILE}" ]]; then
    echo "run_step3_meta.sbatch: missing feature list ${FEATURE_LIST_FILE}" >&2
    exit 1
  fi
  FEATURE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FEATURE_LIST_FILE}")
  if [[ -z "${FEATURE}" ]]; then
    echo "run_step3_meta.sbatch: empty feature for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
  fi
  echo "run_step3_meta.sbatch: array task ${SLURM_ARRAY_TASK_ID}, feature=${FEATURE}"
  FEATURE="${FEATURE}" bash "${SIMU_ROOT}/run_step3_meta.sh"
else
  bash "${SIMU_ROOT}/run_step3_meta.sh"
fi
