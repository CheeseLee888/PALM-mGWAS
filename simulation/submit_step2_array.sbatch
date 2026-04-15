#!/bin/bash
#SBATCH --job-name=simu_step2_submit
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/simu_step2_submit_%j.out
#SBATCH --error=logs/simu_step2_submit_%j.err

set -euo pipefail

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
STUDY="${STUDY:-study1}"
FEATURE_BLOCK="${FEATURE_BLOCK:-10}"
CHROMS="${CHROMS:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"

cd "${SIMU_ROOT}"

FEATURE_FILE="${SIMU_ROOT}/output/${STUDY}/feature_list.txt"
if [[ ! -f "${FEATURE_FILE}" ]]; then
  echo "submit_step2_array: missing feature list: ${FEATURE_FILE}" >&2
  exit 1
fi

N_FEATURE=$(wc -l < "${FEATURE_FILE}")
IFS=',' read -r -a CHROM_LIST <<< "${CHROMS}"
N_CHROM=${#CHROM_LIST[@]}
N_BLOCK=$(( (N_FEATURE + FEATURE_BLOCK - 1) / FEATURE_BLOCK ))
TOTAL_TASKS=$(( N_BLOCK * N_CHROM ))

echo "submit_step2_array: simulation root = ${SIMU_ROOT}"
echo "submit_step2_array: study = ${STUDY}"
echo "submit_step2_array: feature block = ${FEATURE_BLOCK}"
echo "submit_step2_array: n_feature = ${N_FEATURE}"
echo "submit_step2_array: n_block = ${N_BLOCK}"
echo "submit_step2_array: n_chrom = ${N_CHROM}"
echo "submit_step2_array: total tasks = ${TOTAL_TASKS}"

jid=$(sbatch --parsable \
  --array=1-${TOTAL_TASKS} \
  --export=ALL,SIMU_ROOT=${SIMU_ROOT},STUDY=${STUDY},FEATURE_BLOCK=${FEATURE_BLOCK} \
  step2_run_array.sbatch)
echo "submit_step2_array: submitted Step2 array job ${jid} with 1-${TOTAL_TASKS}"
