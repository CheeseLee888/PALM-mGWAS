#!/bin/bash
#SBATCH --job-name=simu_meta_vis_dispatch
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/simu_meta_vis_dispatch_%A_%a.out
#SBATCH --error=logs/simu_meta_vis_dispatch_%A_%a.err

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_SIMU_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SIMU_ROOT="${SIMU_ROOT:-${DEFAULT_SIMU_ROOT}}"
export STUDIES

cd "${SIMU_ROOT}"

IFS=',' read -r -a STUDY_LIST <<< "${STUDIES}"

step2_done_jids=()
for study in "${STUDY_LIST[@]}"; do
  jid_file="${STEP2_DONE_JID_DIR}/${study}.jid"
  jid=$(tr -d '[:space:]' < "${jid_file}")
  step2_done_jids+=("${jid}")
done

step2_done_dep=$(IFS=:; echo "${step2_done_jids[*]}")
STUDY_FILE="${META_DIR}/study_dirs.tsv"

: > "${STUDY_FILE}"
for study in "${STUDY_LIST[@]}"; do
  OUTPUT_DIR="${SIMU_ROOT}/output/${study}"
  printf '%s\t%s\n' "${study}" "${OUTPUT_DIR}" >> "${STUDY_FILE}"
done
echo "dispatch_meta_vis: wrote ${STUDY_FILE}"

n_feature=$(wc -l < "${FEATURE_LIST_FILE}")

jid3=$(sbatch --parsable \
  --dependency=afterok:${step2_done_dep} \
  --array=1-${n_feature} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},META_DIR=${META_DIR},FEATURE_LIST_FILE=${FEATURE_LIST_FILE},STUDY_FILE=${STUDY_FILE} \
  run/step3.sbatch)
echo "dispatch_meta_vis: submitted Step3 meta array job ${jid3} with 1-${n_feature}"

jid4=$(sbatch --parsable \
  --dependency=afterok:${jid3} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},META_DIR=${META_DIR} \
  run/step4.sbatch)
echo "dispatch_meta_vis: submitted Step4 job ${jid4}"
