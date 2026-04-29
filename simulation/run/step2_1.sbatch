#!/bin/bash
#SBATCH --job-name=step2_1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/step2_1_%A_%a.out
#SBATCH --error=logs/step2_1_%A_%a.err
# Step2.1 in the fixed simulation workflow:
# one Slurm array task per chromosome, fixed as task 1..22 -> chr1..22.

set -euo pipefail
WORK="${WORK:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIF="${SIF:-${WORK}/PALMmGWAS.sif}"
cd "${WORK}"

STUDY="${STUDY:-study1}"
MIN_MAF="${MIN_MAF:-0.05}"
MIN_MAC="${MIN_MAC:-5}"
USE_CLUSTER="${USE_CLUSTER:-FALSE}"
CLUSTER_FILE="${CLUSTER_FILE:-NULL}"

INPUT_DIR="${WORK}/input/${STUDY}"
OUTPUT_DIR="${WORK}/output/${STUDY}"
GENO_PREFIX="${INPUT_DIR}/geno"
NULLMODEL_FILE="${OUTPUT_DIR}/step1_allpheno.rda"
FEATURE_FILE="${OUTPUT_DIR}/feature_list.txt"
TOTAL_TASKS=22

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
if (( TASK_ID < 1 || TASK_ID > TOTAL_TASKS )); then
  echo "step2_1_run_array: task ${TASK_ID} out of range for ${TOTAL_TASKS} chromosomes" >&2
  exit 1
fi

STEP2_SUMMARY_PREFIX="${OUTPUT_DIR}/step2"
CHROM="${TASK_ID}"
CHROM_ARG="${CHROM}"

if [[ ! -s "${FEATURE_FILE}" ]]; then
  echo "step2_1_run_array: missing or empty feature list: ${FEATURE_FILE}" >&2
  exit 1
fi

echo "step2_1_run_array: task=${TASK_ID}/${TOTAL_TASKS}"
echo "step2_1_run_array: WORK=${WORK}"
echo "step2_1_run_array: SIF=${SIF}"
echo "step2_1_run_array: step2.1 mode=chrom_all_features"
echo "step2_1_run_array: chrom=${CHROM}"
echo "step2_1_run_array: feature list=${FEATURE_FILE}"
echo "step2_1_run_array: SummaryPrefix=${STEP2_SUMMARY_PREFIX}"
echo "step2_1_run_array: cluster file=${CLUSTER_FILE}"

apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}" \
  step2_1_summary.R \
  --genoFile="${GENO_PREFIX}" \
  --NULLObjPrefix="${NULLMODEL_FILE%.rda}" \
  --SummaryPrefix="${STEP2_SUMMARY_PREFIX}" \
  --chrom="${CHROM_ARG}" \
  --featureColList=NULL \
  --minMAF="${MIN_MAF}" \
  --minMAC="${MIN_MAC}" \
  --SnpInfoFile=NULL \
  --useCluster="${USE_CLUSTER}" \
  --clusterFile="${CLUSTER_FILE}"
