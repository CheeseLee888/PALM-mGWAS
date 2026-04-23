#!/bin/bash
#SBATCH --job-name=step2_2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/step2_2_%A_%a.out
#SBATCH --error=logs/step2_2_%A_%a.err

set -euo pipefail
WORK="${WORK:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIF="${SIF:-${WORK}/PALMmGWAS.sif}"

cd "${WORK}"

OUTPUT_DIR="${WORK}/output/${STUDY}"
INPUT_PREFIX="${OUTPUT_DIR}/step2"
NULLMODEL_FILE="${OUTPUT_DIR}/step1_allpheno.rda"
OVERWRITE_OUTPUT="TRUE"
STEP2_CORRECT_MODE="${STEP2_CORRECT_MODE:-median}"
TOTAL_TASKS=22
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
if (( TASK_ID < 1 || TASK_ID > TOTAL_TASKS )); then
  echo "run_step2_2: task ${TASK_ID} out of range for ${TOTAL_TASKS} chromosomes" >&2
  exit 1
fi
CHROM="${TASK_ID}"

echo "run_step2_2: WORK = ${WORK}"
echo "run_step2_2: study = ${STUDY}"
echo "run_step2_2: task=${TASK_ID}/${TOTAL_TASKS}"
echo "run_step2_2: chrom = ${CHROM}"
echo "run_step2_2: overwriteOutput = ${OVERWRITE_OUTPUT}"
echo "run_step2_2: correct = ${STEP2_CORRECT_MODE}"

echo "run_step2_2: inputPrefix=${INPUT_PREFIX}"
echo "run_step2_2: NULLmodelFile=${NULLMODEL_FILE}"

apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}" \
  step2_2_correction.R \
  --inputPrefix="${INPUT_PREFIX}" \
  --chrom="${CHROM}" \
  --overwriteOutput="${OVERWRITE_OUTPUT}" \
  --correct="${STEP2_CORRECT_MODE}" \
  --NULLmodelFile="${NULLMODEL_FILE}"
