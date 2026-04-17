#!/bin/bash
#SBATCH --job-name=simu_post_step2_1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/simu_post_step2_1_%A_%a.out
#SBATCH --error=logs/simu_post_step2_1_%A_%a.err

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
STUDIES="${STUDIES:-study1,study2,study3}"
META_DIR="${META_DIR:-${SIMU_ROOT}/output/meta}"
FEATURE_LIST_FILE="${FEATURE_LIST_FILE:-${SIMU_ROOT}/output/study1/feature_list.txt}"
STEP2_JID_DIR="${STEP2_JID_DIR:-${SIMU_ROOT}/logs/step2_jids}"

cd "${SIMU_ROOT}"
mkdir -p "${SIMU_ROOT}/logs" "${META_DIR}"

IFS=',' read -r -a STUDY_LIST <<< "${STUDIES}"
if [[ "${#STUDY_LIST[@]}" -eq 0 ]]; then
  echo "submit_post_step2_pipeline: no studies resolved from STUDIES" >&2
  exit 1
fi

step2_array_jids=()
for study in "${STUDY_LIST[@]}"; do
  jid_file="${STEP2_JID_DIR}/${study}.jid"
  if [[ ! -f "${jid_file}" ]]; then
    echo "submit_post_step2_pipeline: missing Step2.1 jid file ${jid_file}" >&2
    exit 1
  fi
  jid=$(tr -d '[:space:]' < "${jid_file}")
  if [[ -z "${jid}" ]]; then
    echo "submit_post_step2_pipeline: empty Step2.1 jid in ${jid_file}" >&2
    exit 1
  fi
  step2_array_jids+=("${jid}")
done

step2_dep=$(IFS=:; echo "${step2_array_jids[*]}")
echo "submit_post_step2_pipeline: Step2.1 array dependencies = ${step2_dep}"

jid_merge=$(sbatch --parsable \
  --dependency=afterok:${step2_dep} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDIES=${STUDIES},META_DIR=${META_DIR} \
  run_step2_2.sbatch)
echo "submit_post_step2_pipeline: Step2.2 merge job ${jid_merge}"

jid_corr=$(sbatch --parsable \
  --dependency=afterok:${jid_merge} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDIES=${STUDIES} \
  run_step2_3_correction.sbatch)
echo "submit_post_step2_pipeline: Step2.3 correction job ${jid_corr}"

if [[ ! -f "${FEATURE_LIST_FILE}" ]]; then
  echo "submit_post_step2_pipeline: missing feature list ${FEATURE_LIST_FILE}" >&2
  exit 1
fi

n_feature=$(wc -l < "${FEATURE_LIST_FILE}")
if [[ "${n_feature}" -lt 1 ]]; then
  echo "submit_post_step2_pipeline: feature list is empty: ${FEATURE_LIST_FILE}" >&2
  exit 1
fi

jid_meta=$(sbatch --parsable \
  --dependency=afterok:${jid_corr} \
  --array=1-${n_feature} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},META_DIR=${META_DIR},FEATURE_LIST_FILE=${FEATURE_LIST_FILE} \
  run_step3_meta.sh)
echo "submit_post_step2_pipeline: step3 meta array job ${jid_meta} with 1-${n_feature}"

jid_vis=$(sbatch --parsable \
  --dependency=afterok:${jid_meta} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},META_DIR=${META_DIR} \
  visualize_meta_results.sh)
echo "submit_post_step2_pipeline: visualization job ${jid_vis}"
