#!/bin/bash
#SBATCH --job-name=step3
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/step3_%A_%a.out
#SBATCH --error=logs/step3_%A_%a.err
# Step3 in the fixed simulation workflow:
# one Slurm array task per feature, and each task runs all chromosomes.

set -euo pipefail
WORK="${WORK:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIF="${SIF:-${WORK}/PALMGWAS.sif}"
META_DIR="${META_DIR:-${WORK}/output/meta}"
INPUT_PREFIX_FILE="${INPUT_PREFIX_FILE:-${META_DIR}/input_prefixes.tsv}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-${META_DIR}/step3_meta}"
FEATURE_INFO_FILE="${FEATURE_INFO_FILE:-${WORK}/output/study1/info_feature.txt}"
STEP3_CHROM="NULL"

cd "${WORK}"

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
if [[ ! -s "${FEATURE_INFO_FILE}" ]]; then
  echo "run_step3_meta: missing or empty feature info file: ${FEATURE_INFO_FILE}" >&2
  exit 1
fi

TOTAL_TASKS="$(( $(wc -l < "${FEATURE_INFO_FILE}" | tr -d '[:space:]') - 1 ))"
if (( TASK_ID < 1 || TASK_ID > TOTAL_TASKS )); then
  echo "run_step3_meta: task ${TASK_ID} out of range for ${TOTAL_TASKS} features" >&2
  exit 1
fi

FEATURE="$(awk -v n="${TASK_ID}" 'NR == n + 1 {print $1}' "${FEATURE_INFO_FILE}")"
if [[ -z "${FEATURE}" ]]; then
  echo "run_step3_meta: empty feature at task ${TASK_ID} in ${FEATURE_INFO_FILE}" >&2
  exit 1
fi

echo "run_step3_meta: WORK = ${WORK}"
echo "run_step3_meta: task=${TASK_ID}/${TOTAL_TASKS}"
echo "run_step3_meta: input prefix file = ${INPUT_PREFIX_FILE}"
echo "run_step3_meta: meta dir = ${META_DIR}"
echo "run_step3_meta: output prefix = ${OUTPUT_PREFIX}"
echo "run_step3_meta: chrom = ${STEP3_CHROM}"
echo "run_step3_meta: feature info = ${FEATURE_INFO_FILE}"
echo "run_step3_meta: feature = ${FEATURE}"

CMD=(
  apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}"
  step3_meta.R
  --inputPrefixFile="${INPUT_PREFIX_FILE}"
  --chrom="${STEP3_CHROM}"
  --outputPrefix="${OUTPUT_PREFIX}"
  --featureList="${FEATURE}"
)

"${CMD[@]}"

echo "run_step3_meta: finished"
