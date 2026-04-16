#!/bin/bash
#SBATCH --job-name=simu_merge_step2
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/simu_merge_step2_%A_%a.out
#SBATCH --error=logs/simu_merge_step2_%A_%a.err

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
cd "${SIMU_ROOT}"
mkdir -p "${SIMU_ROOT}/logs"

echo "merge_step2_outputs.sbatch: simulation root = ${SIMU_ROOT}"
bash "${SIMU_ROOT}/merge_step2_outputs.sh"
