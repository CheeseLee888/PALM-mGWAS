#!/bin/bash
#SBATCH --job-name=simu_step3_submit
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/simu_step3_submit_%A_%a.out
#SBATCH --error=logs/simu_step3_submit_%A_%a.err

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
META_DIR="${META_DIR:-/mnt/scratch/group/ztang2/pli297/simulation/output/meta}"
FEATURE_LIST_FILE="${FEATURE_LIST_FILE:-${SIMU_ROOT}/output/study1/feature_list.txt}"

cd "${SIMU_ROOT}"
mkdir -p "${SIMU_ROOT}/logs" "${META_DIR}"

if [[ ! -f "${META_DIR}/study_dirs.tsv" ]]; then
  echo "submit_step3_meta_array: missing ${META_DIR}/study_dirs.tsv" >&2
  echo "submit_step3_meta_array: run Step2.2 merge first." >&2
  exit 1
fi

if [[ ! -f "${FEATURE_LIST_FILE}" ]]; then
  echo "submit_step3_meta_array: missing feature list ${FEATURE_LIST_FILE}" >&2
  echo "submit_step3_meta_array: run step1_fit_null_model.sbatch first." >&2
  exit 1
fi

N_FEATURE=$(wc -l < "${FEATURE_LIST_FILE}")
if [[ "${N_FEATURE}" -lt 1 ]]; then
  echo "submit_step3_meta_array: feature list is empty: ${FEATURE_LIST_FILE}" >&2
  exit 1
fi

echo "submit_step3_meta_array: simulation root = ${SIMU_ROOT}"
echo "submit_step3_meta_array: meta dir = ${META_DIR}"
echo "submit_step3_meta_array: feature list = ${FEATURE_LIST_FILE}"
echo "submit_step3_meta_array: total features = ${N_FEATURE}"

jid=$(sbatch --parsable \
  --array=1-${N_FEATURE} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},META_DIR=${META_DIR},FEATURE_LIST_FILE=${FEATURE_LIST_FILE} \
  run_step3_meta.sh)
echo "submit_step3_meta_array: submitted Step3 meta array job ${jid} with 1-${N_FEATURE}"
