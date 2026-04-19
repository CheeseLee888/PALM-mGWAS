#!/bin/bash
set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_SIMU_ROOT="${SCRIPT_DIR}"
SIMU_ROOT="${SIMU_ROOT:-${DEFAULT_SIMU_ROOT}}"
STUDIES="${STUDIES:-study1,study2,study3}"
CHROMS="${CHROMS:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"
META_DIR="${META_DIR:-${SIMU_ROOT}/output/meta}"
FEATURE_LIST_FILE="${FEATURE_LIST_FILE:-${SIMU_ROOT}/output/study1/feature_list.txt}"
STEP2_DONE_JID_DIR="${STEP2_DONE_JID_DIR:-${SIMU_ROOT}/logs/step2_done_jids}"
export STUDIES CHROMS

cd "${SIMU_ROOT}"
mkdir -p "${SIMU_ROOT}/logs" "${STEP2_DONE_JID_DIR}"

IFS=',' read -r -a STUDY_LIST <<< "${STUDIES}"

echo "submit_all: simulation root = ${SIMU_ROOT}"
echo "submit_all: studies = ${STUDIES}"
echo "submit_all: step2.1 mode = chrom_all_features"
echo "submit_all: chrom list = ${CHROMS}"

step2_dispatch_jids=()
for study in "${STUDY_LIST[@]}"; do
  jid0=$(sbatch --parsable \
    --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${study} \
    run/step0.sbatch)
  echo "submit_all: ${study} Step0 job ${jid0}"

  jid1=$(sbatch --parsable \
    --dependency=afterok:${jid0} \
    --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${study} \
    run/step1.sbatch)
  echo "submit_all: ${study} Step1 job ${jid1}"

  jid2=$(sbatch --parsable \
    --dependency=afterok:${jid1} \
    --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${study},STEP2_DONE_JID_FILE=${STEP2_DONE_JID_DIR}/${study}.jid \
    run/dispatch_step2.sbatch)
  echo "submit_all: ${study} Step2 dispatcher job ${jid2}"

  step2_dispatch_jids+=("${jid2}")
done

step2_dispatch_dep=$(IFS=:; echo "${step2_dispatch_jids[*]}")
jid_meta_vis=$(sbatch --parsable \
  --dependency=afterok:${step2_dispatch_dep} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},META_DIR=${META_DIR},FEATURE_LIST_FILE=${FEATURE_LIST_FILE},STEP2_DONE_JID_DIR=${STEP2_DONE_JID_DIR} \
  run/dispatch_meta_vis.sbatch)
echo "submit_all: meta_vis dispatcher job ${jid_meta_vis}"
