#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/merge_%A_%a.out
#SBATCH --error=logs/merge_%A_%a.err

set -euo pipefail
WORK="${WORK:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIF="${SIF:-${WORK}/PALMmGWAS.sif}"

cd "${WORK}"

OUTPUT_DIR="${WORK}/output/${STUDY}"
FEATURE_LIST_FILE="${FEATURE_LIST_FILE:-${OUTPUT_DIR}/feature_list.txt}"
FEATURE="${FEATURE:-}"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && -z "${FEATURE}" ]]; then
  if [[ -z "${FEATURE_LIST_FILE:-}" || ! -s "${FEATURE_LIST_FILE}" ]]; then
    echo "run_merge: missing or empty FEATURE_LIST_FILE: ${FEATURE_LIST_FILE:-unset}" >&2
    exit 1
  fi
  FEATURE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FEATURE_LIST_FILE}")
  if [[ -z "${FEATURE}" ]]; then
    echo "run_merge: empty feature for task ${SLURM_ARRAY_TASK_ID}; skipping"
    exit 0
  fi
fi

echo "run_merge: WORK = ${WORK}"
echo "run_merge: study = ${STUDY}"
echo "run_merge: feature = ${FEATURE}"

apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}" \
  step_merge.R \
  --inputPrefix="${OUTPUT_DIR}/step2" \
  --featureColList="${FEATURE:-NULL}" \
  --outputPrefix="${OUTPUT_DIR}/step2"

echo "run_merge: finished"
