#!/bin/bash
#SBATCH --job-name=simu_step2_1_submit
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/simu_step2_1_submit_%A_%a.out
#SBATCH --error=logs/simu_step2_1_submit_%A_%a.err

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
STUDY="${STUDY:-study1}"
STEP2_MODE="${STEP2_MODE:-feature_chrom}"
CHROMS="${CHROMS:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"
STEP2_JID_FILE="${STEP2_JID_FILE:-}"

cd "${SIMU_ROOT}"

FEATURE_FILE="${SIMU_ROOT}/output/${STUDY}/feature_list.txt"
if [[ ! -f "${FEATURE_FILE}" ]]; then
  echo "submit_step2_array: missing feature list: ${FEATURE_FILE}" >&2
  exit 1
fi

N_FEATURE=$(wc -l < "${FEATURE_FILE}")
case "${STEP2_MODE}" in
  feature_allchr)
    TOTAL_TASKS=${N_FEATURE}
    ;;
  feature_chrom)
    IFS=',' read -r -a CHROM_LIST <<< "${CHROMS}"
    N_CHROM=${#CHROM_LIST[@]}
    TOTAL_TASKS=$(( N_FEATURE * N_CHROM ))
    ;;
  *)
    echo "submit_step2_array: invalid STEP2_MODE=${STEP2_MODE}. Use feature_allchr or feature_chrom." >&2
    exit 1
    ;;
esac

echo "submit_step2_array: simulation root = ${SIMU_ROOT}"
echo "submit_step2_array: study = ${STUDY}"
echo "submit_step2_array: step2.1 mode = ${STEP2_MODE}"
echo "submit_step2_array: n_feature = ${N_FEATURE}"
if [[ "${STEP2_MODE}" == "feature_allchr" ]]; then
  echo "submit_step2_array: one array task per feature across all chromosomes"
else
  echo "submit_step2_array: one array task per feature and chromosome"
  echo "submit_step2_array: chrom list = ${CHROMS}"
  echo "submit_step2_array: n_chrom = ${N_CHROM}"
fi
echo "submit_step2_array: total tasks = ${TOTAL_TASKS}"

jid=$(sbatch --parsable \
  --array=1-${TOTAL_TASKS} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${STUDY},STEP2_MODE=${STEP2_MODE},CHROMS=${CHROMS} \
  step2_run_array.sbatch)
echo "submit_step2_array: submitted Step2.1 array job ${jid} with 1-${TOTAL_TASKS}"
if [[ -n "${STEP2_JID_FILE}" ]]; then
  mkdir -p "$(dirname "${STEP2_JID_FILE}")"
  printf '%s\n' "${jid}" > "${STEP2_JID_FILE}"
  echo "submit_step2_array: wrote Step2.1 array job id to ${STEP2_JID_FILE}"
fi
