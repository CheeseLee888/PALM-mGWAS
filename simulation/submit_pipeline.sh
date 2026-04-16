#!/bin/bash
set -euo pipefail

# One-command Slurm submission for simulation analysis only:
# Step0 -> Step1 -> Step2.1 array.
# Step2.1 is submitted by a lightweight bridge job after Step1 finishes.

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
STUDY="${STUDY:-study1}"
STEP2_MODE="${STEP2_MODE:-feature_chrom}"
CHROMS="${CHROMS:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"

cd "${SIMU_ROOT}"

echo "submit_pipeline: simulation root = ${SIMU_ROOT}"
echo "submit_pipeline: study = ${STUDY}"
echo "submit_pipeline: step2.1 mode = ${STEP2_MODE}"
echo "submit_pipeline: chrom list = ${CHROMS}"

jid0=$(sbatch --parsable \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${STUDY} \
  step0_align_inputs.sbatch)
echo "submit_pipeline: submitted Step0 job ${jid0}"

jid1=$(sbatch --parsable \
  --dependency=afterok:${jid0} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${STUDY} \
  step1_fit_null_model.sbatch)
echo "submit_pipeline: submitted Step1 job ${jid1}"

jid2=$(sbatch --parsable \
  --dependency=afterok:${jid1} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${STUDY},STEP2_MODE=${STEP2_MODE} \
  submit_step2_array.sbatch)
echo "submit_pipeline: submitted Step2.1 bridge job ${jid2}"
