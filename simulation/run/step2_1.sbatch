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
SIF="${SIF:-${WORK}/PALMGWAS.sif}"
cd "${WORK}"

STUDY="${STUDY:-study1}"
MIN_MAF="${MIN_MAF:-0.05}"
MIN_MAC="${MIN_MAC:-5}"
MAX_MISSING="${MAX_MISSING:-0.15}"
IMPUTE_METHOD="${IMPUTE_METHOD:-best_guess}"

INPUT_DIR="${WORK}/input/${STUDY}"
OUTPUT_DIR="${WORK}/output/${STUDY}"
GENO_PREFIX="${INPUT_DIR}/geno"
nullModelPrefix="${OUTPUT_DIR}/step1_allpheno"
FEATURE_INFO_FILE="${OUTPUT_DIR}/info_feature.txt"
STEP2_OUTPUT_PREFIX="${OUTPUT_DIR}/step2"

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
if [[ ! -s "${FEATURE_INFO_FILE}" ]]; then
  echo "step2_1: missing or empty feature info file: ${FEATURE_INFO_FILE}" >&2
  exit 1
fi

TOTAL_TASKS="$(( $(wc -l < "${FEATURE_INFO_FILE}" | tr -d '[:space:]') - 1 ))"
if (( TASK_ID < 1 || TASK_ID > TOTAL_TASKS )); then
  echo "step2_1: task ${TASK_ID} out of range for ${TOTAL_TASKS} features" >&2
  exit 1
fi

FEATURE="$(awk -v n="${TASK_ID}" 'NR == n + 1 {print $1}' "${FEATURE_INFO_FILE}")"
if [[ -z "${FEATURE}" ]]; then
  echo "step2_1: empty feature at task ${TASK_ID} in ${FEATURE_INFO_FILE}" >&2
  exit 1
fi

echo "step2_1: task=${TASK_ID}/${TOTAL_TASKS}"
echo "step2_1: WORK=${WORK}"
echo "step2_1: SIF=${SIF}"
echo "step2_1: step2.1 mode=feature_all_chromosomes"
echo "step2_1: chrom=NULL (all chromosomes)"
echo "step2_1: feature=${FEATURE}"
echo "step2_1: feature info=${FEATURE_INFO_FILE}"
echo "step2_1: outputPrefix=${STEP2_OUTPUT_PREFIX}"

apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}" \
  step2_1_summary.R \
  --genoFile="${GENO_PREFIX}" \
  --nullModelPrefix="${nullModelPrefix}" \
  --outputPrefix="${STEP2_OUTPUT_PREFIX}" \
  --chrom=NULL \
  --featureList="${FEATURE}" \
  --minMAF="${MIN_MAF}" \
  --minMAC="${MIN_MAC}" \
  --maxMissing="${MAX_MISSING}" \
  --impute_method="${IMPUTE_METHOD}" \
  --snpInfoFile=NULL
