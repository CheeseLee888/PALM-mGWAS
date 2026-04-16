#!/bin/bash
#SBATCH --job-name=simu_merge_meta
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/simu_merge_meta_%A_%a.out
#SBATCH --error=logs/simu_merge_meta_%A_%a.err
#
# Usage:
#   sbatch --export=ALL,ENV_FILE=/path/to/your.env simulation/merge_step2_meta.sbatch
#
# Optional environment overrides:
#   ENV_FILE=/path/to/your.env
#   SIMU_ROOT=/path/to/simulation
#   SIF=/path/to/PALMmGWAS.sif
#   STUDY=study1
#   STUDIES=study1,study2,study3
#   META_DIR=/path/to/output/meta
#   PLOT_DIR=/path/to/output/plot
#   META_PREFIX=step3_meta_

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
cd "${SIMU_ROOT}"
mkdir -p "${SIMU_ROOT}/logs"

echo "merge_step2_meta.sbatch: simulation root = ${SIMU_ROOT}"
if [[ -n "${STUDIES:-}" ]]; then
  echo "merge_step2_meta.sbatch: studies = ${STUDIES}"
elif [[ -n "${STUDY:-}" ]]; then
  echo "merge_step2_meta.sbatch: study = ${STUDY}"
else
  echo "merge_step2_meta.sbatch: studies = study1,study2,study3"
fi

bash "${SIMU_ROOT}/merge_step2_meta.sh"
