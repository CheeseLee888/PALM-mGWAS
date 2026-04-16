#!/bin/bash
set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
STUDIES="${STUDIES:-study1,study2,study3}"
STEP2_MODE="${STEP2_MODE:-feature_chrom}"
CHROMS="${CHROMS:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"
META_DIR="${META_DIR:-${SIMU_ROOT}/output/meta}"
FEATURE_LIST_FILE="${FEATURE_LIST_FILE:-${SIMU_ROOT}/output/study1/feature_list.txt}"
STEP2_JID_DIR="${STEP2_JID_DIR:-${SIMU_ROOT}/logs/step2_jids}"

cd "${SIMU_ROOT}"
mkdir -p "${STEP2_JID_DIR}"

IFS=',' read -r -a STUDY_LIST <<< "${STUDIES}"
if [[ "${#STUDY_LIST[@]}" -eq 0 ]]; then
  echo "submit_all: no studies resolved from STUDIES" >&2
  exit 1
fi

echo "submit_all: simulation root = ${SIMU_ROOT}"
echo "submit_all: studies = ${STUDIES}"
echo "submit_all: step2.1 mode = ${STEP2_MODE}"
echo "submit_all: chrom list = ${CHROMS}"

step2_submit_jids=()
for study in "${STUDY_LIST[@]}"; do
  jid0=$(sbatch --parsable \
    --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${study} \
    step0_align_inputs.sbatch)
  echo "submit_all: ${study} Step0 job ${jid0}"

  jid1=$(sbatch --parsable \
    --dependency=afterok:${jid0} \
    --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${study} \
    step1_fit_null_model.sbatch)
  echo "submit_all: ${study} Step1 job ${jid1}"

  jid2=$(sbatch --parsable \
    --dependency=afterok:${jid1} \
    --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${study},STEP2_MODE=${STEP2_MODE},CHROMS=${CHROMS},STEP2_JID_FILE=${STEP2_JID_DIR}/${study}.jid \
    submit_step2_array.sbatch)
  echo "submit_all: ${study} Step2.1 submit job ${jid2}"

  step2_submit_jids+=("${jid2}")
done

step2_submit_dep=$(IFS=:; echo "${step2_submit_jids[*]}")
jid_post=$(sbatch --parsable \
  --dependency=afterok:${step2_submit_dep} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDIES=${STUDIES},META_DIR=${META_DIR},FEATURE_LIST_FILE=${FEATURE_LIST_FILE},STEP2_JID_DIR=${STEP2_JID_DIR} \
  submit_post_step2_pipeline.sbatch)
echo "submit_all: post-Step2.1 submitter job ${jid_post}"
