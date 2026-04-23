#!/bin/bash
#SBATCH --job-name=step3
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/step3_%A_%a.out
#SBATCH --error=logs/step3_%A_%a.err

set -euo pipefail
WORK="${WORK:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIF="${SIF:-${WORK}/PALMmGWAS.sif}"
META_DIR="${META_DIR:-${WORK}/output/meta}"
FEATURE_LIST_FILE="${FEATURE_LIST_FILE:-${WORK}/output/study1/feature_list.txt}"
STUDY_FILE="${STUDY_FILE:-${META_DIR}/study_dirs.tsv}"
FEATURE="${FEATURE:-}"
META_PREFIX="step3_meta_"
STEP3_INPUT_PREFIX="step2"
STEP3_CHROM="NULL"

cd "${WORK}"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && -z "${FEATURE}" ]]; then
  FEATURE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FEATURE_LIST_FILE}")
  if [[ -z "${FEATURE}" ]]; then
    echo "run_step3_meta: empty feature for task ${SLURM_ARRAY_TASK_ID}; skipping"
    exit 0
  fi
fi

echo "run_step3_meta: WORK = ${WORK}"
echo "run_step3_meta: study dir file = ${STUDY_FILE}"
echo "run_step3_meta: meta dir = ${META_DIR}"
echo "run_step3_meta: input prefix = ${STEP3_INPUT_PREFIX}"
echo "run_step3_meta: chrom = ${STEP3_CHROM}"
echo "run_step3_meta: feature = ${FEATURE}"

CMD=(
  apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}"
  step3_meta.R
  --studyDirFile="${STUDY_FILE}"
  --inputPrefix="${STEP3_INPUT_PREFIX}"
  --chrom="${STEP3_CHROM}"
  --metaPrefix="${META_DIR}/${META_PREFIX}"
)

CMD+=(--featureColList="${FEATURE:-NULL}")

"${CMD[@]}"

echo "run_step3_meta: finished"
