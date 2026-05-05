#!/bin/bash
#SBATCH --job-name=step2_1_feat
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/step2_1_%A_%a.out
#SBATCH --error=logs/step2_1_%A_%a.err
# Step2.1 in the fixed simulation workflow:
# one Slurm array task per feature, and each task runs all chromosomes.

set -euo pipefail
WORK="${WORK:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIF="${SIF:-${WORK}/PALMmGWAS.sif}"
cd "${WORK}"

STUDY="${STUDY:-study1}"
MIN_MAF="${MIN_MAF:-0.05}"
MIN_MAC="${MIN_MAC:-5}"
MAX_MISSING="${MAX_MISSING:-0.15}"
IMPUTE_METHOD="${IMPUTE_METHOD:-best_guess}"
USE_CLUSTER="${USE_CLUSTER:-FALSE}"
CLUSTER_FILE="${CLUSTER_FILE:-NULL}"

INPUT_DIR="${WORK}/input/${STUDY}"
OUTPUT_DIR="${WORK}/output/${STUDY}"
GENO_PREFIX="${INPUT_DIR}/geno"
NULLMODEL_FILE="${OUTPUT_DIR}/step1_allpheno.rda"
FEATURE_FILE="${OUTPUT_DIR}/feature_list.txt"
STEP2_SUMMARY_PREFIX="${OUTPUT_DIR}/step2"

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
if [[ ! -s "${FEATURE_FILE}" ]]; then
  echo "step2_1: missing or empty feature list: ${FEATURE_FILE}" >&2
  exit 1
fi

TOTAL_TASKS="$(wc -l < "${FEATURE_FILE}" | tr -d '[:space:]')"
if (( TASK_ID < 1 || TASK_ID > TOTAL_TASKS )); then
  echo "step2_1: task ${TASK_ID} out of range for ${TOTAL_TASKS} features" >&2
  exit 1
fi

FEATURE="$(sed -n "${TASK_ID}p" "${FEATURE_FILE}")"
if [[ -z "${FEATURE}" ]]; then
  echo "step2_1: empty feature at line ${TASK_ID} in ${FEATURE_FILE}" >&2
  exit 1
fi

echo "step2_1: task=${TASK_ID}/${TOTAL_TASKS}"
echo "step2_1: WORK=${WORK}"
echo "step2_1: SIF=${SIF}"
echo "step2_1: step2.1 mode=feature_all_chromosomes"
echo "step2_1: chrom=NULL (all chromosomes)"
echo "step2_1: feature=${FEATURE}"
echo "step2_1: feature list=${FEATURE_FILE}"
echo "step2_1: SummaryPrefix=${STEP2_SUMMARY_PREFIX}"
echo "step2_1: cluster file=${CLUSTER_FILE}"

apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}" \
  step2_1_summary.R \
  --genoFile="${GENO_PREFIX}" \
  --NULLObjPrefix="${NULLMODEL_FILE%.rda}" \
  --SummaryPrefix="${STEP2_SUMMARY_PREFIX}" \
  --chrom=NULL \
  --featureColList="${FEATURE}" \
  --minMAF="${MIN_MAF}" \
  --minMAC="${MIN_MAC}" \
  --maxMissing="${MAX_MISSING}" \
  --impute_method="${IMPUTE_METHOD}" \
  --SnpInfoFile=NULL \
  --useCluster="${USE_CLUSTER}" \
  --clusterFile="${CLUSTER_FILE}"
