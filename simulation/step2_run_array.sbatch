#!/bin/bash
#SBATCH --job-name=simu_step2_1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/simu_step2_1_%A_%a.out
#SBATCH --error=logs/simu_step2_1_%A_%a.err
#
# Usage:
#   1. Run simulation/generate_simulation_data.R, Step0, and Step1 once first.
#   2. Generate feature list from the fitted null model, e.g.
#      pixi run --manifest-path=pixi.toml Rscript -e '
#        load("simulation/output/study1/step1_allpheno.rda");
#        writeLines(colnames(modglmm[[1]]$Y_I), "simulation/output/study1/feature_list.txt")
#      '
#   3. Set STEP2_MODE and ARRAY_RANGE, then submit Step2.1:
#      sbatch --array=1-ARRAY_RANGE --export=ALL,ENV_FILE=/path/to/your.env simulation/step2_run_array.sbatch
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
#   STEP2_MODE=feature_chrom or feature_allchr
#   CHROMS=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22
#   MIN_MAF=0.05
#   MIN_MAC=5
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
STEP2_MODE="${STEP2_MODE:-feature_chrom}"
CHROMS_CSV="${CHROMS:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"
MIN_MAF="${MIN_MAF:-0.05}"
MIN_MAC="${MIN_MAC:-5}"
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

N_FEATURE=$(wc -l < "${FEATURE_FILE}")
if [[ "${N_FEATURE}" -lt 1 ]]; then
  echo "Feature list is empty: ${FEATURE_FILE}" >&2
  exit 1
fi

IFS=',' read -r -a CHROM_LIST <<< "${CHROMS_CSV}"
N_CHROM=${#CHROM_LIST[@]}
if [[ "${N_CHROM}" -lt 1 ]]; then
  echo "No chromosomes provided in CHROMS." >&2
  exit 1
fi

case "${STEP2_MODE}" in
  feature_allchr)
    TOTAL_TASKS=${N_FEATURE}
    ;;
  feature_chrom)
    TOTAL_TASKS=$(( N_FEATURE * N_CHROM ))
    ;;
  *)
    echo "Invalid STEP2_MODE=${STEP2_MODE}. Use feature_allchr or feature_chrom." >&2
    exit 1
    ;;
esac

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

STEP2_SUMMARY_PREFIX="${OUTPUT_DIR}/step2"

if [[ "${STEP2_MODE}" == "feature_allchr" ]]; then
  FEATURE_LINE=${TASK_ID}
  FEATURE=$(sed -n "${FEATURE_LINE}p" "${FEATURE_FILE}")
  if [[ -z "${FEATURE}" ]]; then
    echo "Resolved feature is empty for task ${TASK_ID}." >&2
    exit 1
  fi
  CHROM_ARG="NULL"
else
  FEATURE_INDEX=$(( ZERO_BASED / N_CHROM + 1 ))
  CHROM_INDEX=$(( ZERO_BASED % N_CHROM ))
  FEATURE=$(sed -n "${FEATURE_INDEX}p" "${FEATURE_FILE}")
  CHROM="${CHROM_LIST[$CHROM_INDEX]}"
  if [[ -z "${FEATURE}" ]]; then
    echo "Resolved feature is empty for task ${TASK_ID}." >&2
    exit 1
  fi
  CHROM_ARG="${CHROM}"
fi

if [[ -z "${FEATURE}" ]]; then
  echo "Resolved feature is empty for task ${TASK_ID}." >&2
  exit 1
fi

echo "step2_1_run_array: task=${TASK_ID}/${TOTAL_TASKS}"
echo "step2_1_run_array: SIMU_ROOT=${SIMU_ROOT}"
echo "step2_1_run_array: SIF=${SIF}"
echo "step2_1_run_array: step2.1 mode=${STEP2_MODE}"
echo "step2_1_run_array: feature=${FEATURE}"
if [[ "${STEP2_MODE}" == "feature_allchr" ]]; then
  echo "step2_1_run_array: running all chromosomes"
else
  echo "step2_1_run_array: chrom=${CHROM}"
fi
echo "step2_1_run_array: SummaryPrefix=${STEP2_SUMMARY_PREFIX}"

apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  step2_1_summary.R \
  --genoFile="${GENO_PREFIX}" \
  --NULLObjPrefix="${NULLMODEL_FILE%.rda}" \
  --SummaryPrefix="${STEP2_SUMMARY_PREFIX}" \
  --chrom="${CHROM_ARG}" \
  --featureColList="${FEATURE}" \
  --minMAF="${MIN_MAF}" \
  --minMAC="${MIN_MAC}" \
  --SnpInfoFile=NULL \
  --useCluster="${USE_CLUSTER}"
