#!/bin/bash
#SBATCH --job-name=simu_step2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/simu_step2_%A_%a.out
#SBATCH --error=logs/simu_step2_%A_%a.err
#
# Usage:
#   1. Run simulation/generate_simulation_data.R, Step0, and Step1 once first.
#   2. Generate feature list from the fitted null model, e.g.
#      pixi run --manifest-path=pixi.toml Rscript -e '
#        load("simulation/output/study1/step1_allpheno.rda");
#        writeLines(colnames(modglmm[[1]]$Y_I), "simulation/output/study1/feature_list.txt")
#      '
#   3. Set FEATURE_BLOCK and ARRAY_RANGE, then submit:
#      sbatch --array=1-ARRAY_RANGE --export=ALL,ENV_FILE=/path/to/your.env simulation/step2_simu_array.sbatch
#
# Simulation root layout:
#   /mnt/scratch/group/ztang2/pli297/simulation/
#     ├── input/
#     ├── output/
#     ├── feature_list.txt (generated under output/study1/)
#     └── PALMmGWAS.sif
#
# Optional environment overrides:
#   ENV_FILE=/path/to/your.env
#   SIMU_ROOT=/path/to/simulation
#   SIF=/path/to/PALMmGWAS.sif
#   STUDY=study1
#   FEATURE_BLOCK=10
#   CHROMS=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22
#   MIN_MAF=0.05
#   MIN_MAC=5
#   CORRECT=NULL
#   USE_CLUSTER=FALSE

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi
DEFAULT_SIMU_ROOT="/mnt/scratch/group/ztang2/pli297/simulation"
SIMU_ROOT="${SIMU_ROOT:-${PWD:-${DEFAULT_SIMU_ROOT}}}"
SIF="${SIF:-${SIMU_ROOT}/PALMmGWAS.sif}"
cd "${SIMU_ROOT}"

STUDY="${STUDY:-study1}"
FEATURE_BLOCK="${FEATURE_BLOCK:-10}"
CHROMS_CSV="${CHROMS:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"
MIN_MAF="${MIN_MAF:-0.05}"
MIN_MAC="${MIN_MAC:-5}"
CORRECT="${CORRECT:-NULL}"
USE_CLUSTER="${USE_CLUSTER:-FALSE}"

INPUT_DIR="${SIMU_ROOT}/input/${STUDY}"
OUTPUT_DIR="${SIMU_ROOT}/output/${STUDY}"
GENO_PREFIX="${INPUT_DIR}/geno"
NULLMODEL_FILE="${OUTPUT_DIR}/step1_allpheno.rda"
FEATURE_FILE="${OUTPUT_DIR}/feature_list.txt"

mkdir -p "${SIMU_ROOT}/logs" "${OUTPUT_DIR}"

if [[ ! -f "${NULLMODEL_FILE}" ]]; then
  echo "Missing null model: ${NULLMODEL_FILE}" >&2
  exit 1
fi
if [[ ! -f "${FEATURE_FILE}" ]]; then
  echo "Missing feature list: ${FEATURE_FILE}" >&2
  echo "Create it first from Step1 null model." >&2
  exit 1
fi
if ! [[ "${FEATURE_BLOCK}" =~ ^[1-9][0-9]*$ ]]; then
  echo "FEATURE_BLOCK must be a positive integer. Current value: ${FEATURE_BLOCK}" >&2
  exit 1
fi

IFS=',' read -r -a CHROM_LIST <<< "${CHROMS_CSV}"
N_CHROM=${#CHROM_LIST[@]}
if [[ "${N_CHROM}" -lt 1 ]]; then
  echo "No chromosomes provided in CHROMS." >&2
  exit 1
fi

N_FEATURE=$(wc -l < "${FEATURE_FILE}")
if [[ "${N_FEATURE}" -lt 1 ]]; then
  echo "Feature list is empty: ${FEATURE_FILE}" >&2
  exit 1
fi

N_BLOCK=$(( (N_FEATURE + FEATURE_BLOCK - 1) / FEATURE_BLOCK ))
TOTAL_TASKS=$(( N_CHROM * N_BLOCK ))

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
if ! [[ "${TASK_ID}" =~ ^[1-9][0-9]*$ ]]; then
  echo "SLURM_ARRAY_TASK_ID must be a positive integer. Current value: ${TASK_ID}" >&2
  exit 1
fi
if (( TASK_ID > TOTAL_TASKS )); then
  echo "SLURM_ARRAY_TASK_ID=${TASK_ID} exceeds total tasks=${TOTAL_TASKS}." >&2
  exit 1
fi

ZERO_BASED=$(( TASK_ID - 1 ))
CHROM_INDEX=$(( ZERO_BASED / N_BLOCK ))
BLOCK_INDEX=$(( ZERO_BASED % N_BLOCK ))

CHROM="${CHROM_LIST[$CHROM_INDEX]}"
START_LINE=$(( BLOCK_INDEX * FEATURE_BLOCK + 1 ))
END_LINE=$(( START_LINE + FEATURE_BLOCK - 1 ))
if (( END_LINE > N_FEATURE )); then
  END_LINE=${N_FEATURE}
fi

FEATURES=$(sed -n "${START_LINE},${END_LINE}p" "${FEATURE_FILE}" | paste -sd, -)
if [[ -z "${FEATURES}" ]]; then
  echo "Resolved feature block is empty for task ${TASK_ID}." >&2
  exit 1
fi

STEP2_PREFIX="${OUTPUT_DIR}/step2_chr${CHROM}"

echo "step2_simu_array: task=${TASK_ID}/${TOTAL_TASKS}"
echo "step2_simu_array: SIMU_ROOT=${SIMU_ROOT}"
echo "step2_simu_array: SIF=${SIF}"
echo "step2_simu_array: chrom=${CHROM}, block=${BLOCK_INDEX}+1/${N_BLOCK}, feature lines=${START_LINE}-${END_LINE}"
echo "step2_simu_array: outputPrefix=${STEP2_PREFIX}"
echo "step2_simu_array: featureColList=${FEATURES}"

apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  step2_summary.R \
  --genoFile="${GENO_PREFIX}" \
  --NULLmodelFile="${NULLMODEL_FILE}" \
  --PALMOutputFile="${STEP2_PREFIX}" \
  --chrom="${CHROM}" \
  --featureColList="${FEATURES}" \
  --minMAF="${MIN_MAF}" \
  --minMAC="${MIN_MAC}" \
  --outputSnpFile=NULL \
  --correct="${CORRECT}" \
  --useCluster="${USE_CLUSTER}"
