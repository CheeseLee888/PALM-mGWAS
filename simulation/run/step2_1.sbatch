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
#   3. Compute ARRAY_RANGE for chromosome-level tasks, then submit Step2.1:
#      sbatch --array=1-ARRAY_RANGE --export=ALL,ENV_FILE=/path/to/your.env simulation/run/step2_1.sbatch
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
#   CHROMS=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22
#   MIN_MAF=0.05
#   MIN_MAC=5
#   USE_CLUSTER=FALSE

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_SIMU_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SIMU_ROOT="${SIMU_ROOT:-${DEFAULT_SIMU_ROOT}}"
SIF="${SIF:-${SIMU_ROOT}/PALMmGWAS.sif}"
cd "${SIMU_ROOT}"

STUDY="${STUDY:-study1}"
CHROMS_CSV="${CHROMS:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"
MIN_MAF="${MIN_MAF:-0.05}"
MIN_MAC="${MIN_MAC:-5}"
USE_CLUSTER="${USE_CLUSTER:-FALSE}"

INPUT_DIR="${SIMU_ROOT}/input/${STUDY}"
OUTPUT_DIR="${SIMU_ROOT}/output/${STUDY}"
GENO_PREFIX="${INPUT_DIR}/geno"
NULLMODEL_FILE="${OUTPUT_DIR}/step1_allpheno.rda"
FEATURE_FILE="${OUTPUT_DIR}/feature_list.txt"

IFS=',' read -r -a CHROM_LIST <<< "${CHROMS_CSV}"
N_CHROM=${#CHROM_LIST[@]}
TOTAL_TASKS=${N_CHROM}

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

ZERO_BASED=$(( TASK_ID - 1 ))

STEP2_SUMMARY_PREFIX="${OUTPUT_DIR}/step2"
CHROM_INDEX=${ZERO_BASED}
CHROM="${CHROM_LIST[$CHROM_INDEX]}"
CHROM_ARG="${CHROM}"

echo "step2_1_run_array: task=${TASK_ID}/${TOTAL_TASKS}"
echo "step2_1_run_array: SIMU_ROOT=${SIMU_ROOT}"
echo "step2_1_run_array: SIF=${SIF}"
echo "step2_1_run_array: step2.1 mode=chrom_all_features"
echo "step2_1_run_array: chrom=${CHROM}"
echo "step2_1_run_array: SummaryPrefix=${STEP2_SUMMARY_PREFIX}"

apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  step2_1_summary.R \
  --genoFile="${GENO_PREFIX}" \
  --NULLObjPrefix="${NULLMODEL_FILE%.rda}" \
  --SummaryPrefix="${STEP2_SUMMARY_PREFIX}" \
  --chrom="${CHROM_ARG}" \
  --featureColList=NULL \
  --minMAF="${MIN_MAF}" \
  --minMAC="${MIN_MAC}" \
  --SnpInfoFile=NULL \
  --useCluster="${USE_CLUSTER}"
