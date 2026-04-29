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
STUDY_FILE="${STUDY_FILE:-${META_DIR}/study_dirs.tsv}"
META_PREFIX="step3_meta_"
STEP3_INPUT_PREFIX="step2"
STEP3_CHROM="${STEP3_CHROM:-NULL}"

cd "${WORK}"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && "${STEP3_CHROM}" == "NULL" ]]; then
  STEP3_CHROM="${SLURM_ARRAY_TASK_ID}"
fi

if [[ "${STEP3_CHROM}" != "NULL" && ! "${STEP3_CHROM}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
  echo "run_step3_meta: STEP3_CHROM must be NULL or 1..22; got ${STEP3_CHROM}" >&2
  exit 1
fi

echo "run_step3_meta: WORK = ${WORK}"
echo "run_step3_meta: study dir file = ${STUDY_FILE}"
echo "run_step3_meta: meta dir = ${META_DIR}"
echo "run_step3_meta: input prefix = ${STEP3_INPUT_PREFIX}"
echo "run_step3_meta: chrom = ${STEP3_CHROM}"
echo "run_step3_meta: featureColList = NULL (all features)"

CMD=(
  apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}"
  step3_meta.R
  --studyDirFile="${STUDY_FILE}"
  --inputPrefix="${STEP3_INPUT_PREFIX}"
  --chrom="${STEP3_CHROM}"
  --metaPrefix="${META_DIR}/${META_PREFIX}"
  --featureColList=NULL
)

"${CMD[@]}"

echo "run_step3_meta: finished"
