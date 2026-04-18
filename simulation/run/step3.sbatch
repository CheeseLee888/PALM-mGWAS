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
META_DIR="${META_DIR:-/mnt/scratch/group/ztang2/pli297/simulation/output/meta}"
META_PREFIX="${META_PREFIX:-step3_meta_}"
STUDY_FILE="${STUDY_FILE:-${META_DIR}/study_dirs.tsv}"
FEATURE_LIST_FILE="${FEATURE_LIST_FILE:-${SIMU_ROOT}/output/study1/feature_list.txt}"
FEATURE="${FEATURE:-}"

cd "${SIMU_ROOT}"
mkdir -p "${SIMU_ROOT}/logs" "${META_DIR}"

if [[ ! -f "${STUDY_FILE}" ]]; then
  echo "run_step3_meta: missing study dir file ${STUDY_FILE}" >&2
  echo "run_step3_meta: run Step2.2 merge first." >&2
  exit 1
fi

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && -z "${FEATURE}" ]]; then
  if [[ ! -f "${FEATURE_LIST_FILE}" ]]; then
    echo "run_step3_meta: missing feature list ${FEATURE_LIST_FILE}" >&2
    exit 1
  fi
  FEATURE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FEATURE_LIST_FILE}")
  if [[ -z "${FEATURE}" ]]; then
    echo "run_step3_meta: empty feature for task ${SLURM_ARRAY_TASK_ID}; skipping"
    exit 0
  fi
fi

if [[ -z "${FEATURE}" ]]; then
  echo "run_step3_meta: FEATURE is required." >&2
  echo "run_step3_meta: use submit_step3_meta_array.sbatch or set FEATURE=... explicitly." >&2
  exit 1
fi

echo "run_step3_meta: simulation root = ${SIMU_ROOT}"
echo "run_step3_meta: study dir file = ${STUDY_FILE}"
echo "run_step3_meta: meta dir = ${META_DIR}"
echo "run_step3_meta: feature = ${FEATURE}"

matched=0
while IFS=$'\t' read -r _study_id study_dir _rest; do
  [[ -z "${study_dir:-}" ]] && continue
  if [[ -f "${study_dir}/step2_allchr_${FEATURE}.txt" ]]; then
    matched=1
    break
  fi
done < "${STUDY_FILE}"

if [[ "${matched}" -eq 0 ]]; then
  echo "run_step3_meta: no matched step2_allchr file found for feature ${FEATURE}; skipping"
  exit 0
fi

CMD=(
  apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}"
  step3_meta.R
  --studyDirFile="${STUDY_FILE}"
  --pattern='step2_allchr_.*\.txt$'
  --metaDir="${META_DIR}"
  --metaPrefix="${META_PREFIX}"
)

CMD+=(--features="${FEATURE}")

"${CMD[@]}"

echo "run_step3_meta: finished"
