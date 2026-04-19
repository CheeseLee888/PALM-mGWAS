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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_SIMU_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SIMU_ROOT="${SIMU_ROOT:-${DEFAULT_SIMU_ROOT}}"
SIF="${SIF:-${SIMU_ROOT}/PALMmGWAS.sif}"
META_PREFIX="${META_PREFIX:-step3_meta_}"
FEATURE="${FEATURE:-}"

cd "${SIMU_ROOT}"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && -z "${FEATURE}" ]]; then
  FEATURE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FEATURE_LIST_FILE}")
  if [[ -z "${FEATURE}" ]]; then
    echo "run_step3_meta: empty feature for task ${SLURM_ARRAY_TASK_ID}; skipping"
    exit 0
  fi
fi

echo "run_step3_meta: simulation root = ${SIMU_ROOT}"
echo "run_step3_meta: study dir file = ${STUDY_FILE}"
echo "run_step3_meta: meta dir = ${META_DIR}"
echo "run_step3_meta: feature = ${FEATURE}"

CMD=(
  apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}"
  step3_meta.R
  --studyDirFile="${STUDY_FILE}"
  --pattern='step2_allchr_.*\.txt$'
  --metaPrefix="${META_DIR}/${META_PREFIX}"
)

CMD+=(--features="${FEATURE}")

"${CMD[@]}"

echo "run_step3_meta: finished"
