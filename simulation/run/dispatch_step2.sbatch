#!/bin/bash
#SBATCH --job-name=simu_step2_dispatch
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/simu_step2_dispatch_%A_%a.out
#SBATCH --error=logs/simu_step2_dispatch_%A_%a.err

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_SIMU_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SIMU_ROOT="${SIMU_ROOT:-${DEFAULT_SIMU_ROOT}}"
CHROMS="${CHROMS:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"
export CHROMS

cd "${SIMU_ROOT}"

IFS=',' read -r -a CHROM_LIST <<< "${CHROMS}"
N_CHROM=${#CHROM_LIST[@]}
TOTAL_TASKS=${N_CHROM}

echo "dispatch_step2: study = ${STUDY}"
echo "dispatch_step2: total Step2.1 tasks = ${TOTAL_TASKS} (one chromosome per task)"

jid2_1=$(sbatch --parsable \
  --array=1-${TOTAL_TASKS} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${STUDY} \
  run/step2_1.sbatch)
echo "dispatch_step2: submitted Step2.1 array job ${jid2_1} with 1-${TOTAL_TASKS}"

jid_merge=$(sbatch --parsable \
  --dependency=afterok:${jid2_1} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${STUDY} \
  run/step_merge.sbatch)
echo "dispatch_step2: submitted merge job ${jid_merge}"

jid2_2_correct=$(sbatch --parsable \
  --dependency=afterok:${jid_merge} \
  --export=ALL,ENV_FILE=${ENV_FILE:-},SIMU_ROOT=${SIMU_ROOT},STUDY=${STUDY},OVERWRITE_OUTPUT=TRUE \
  run/step2_2_correct.sbatch)
echo "dispatch_step2: submitted Step2.2 correction job ${jid2_2_correct}"

printf '%s\n' "${jid2_2_correct}" > "${STEP2_DONE_JID_FILE}"
echo "dispatch_step2: wrote Step2.2 jid to ${STEP2_DONE_JID_FILE}"
