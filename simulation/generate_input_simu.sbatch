#!/bin/bash
#SBATCH --job-name=simu_input
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/simu_input_%j.out
#SBATCH --error=logs/simu_input_%j.err
#
# Usage:
#   sbatch simulation/generate_input_simu.sbatch
#
# Simulation root layout:
#   /mnt/scratch/group/ztang2/pli297/simulation/
#     ├── generate_simulation_data.R
#     ├── input/
#     ├── output/
#     ├── provided/
#     └── PALMmGWAS.sif
#
# Optional environment overrides:
#   ENV_FILE=/path/to/your.env
#   SIMU_ROOT=/path/to/simulation
#   SIF=/path/to/PALMmGWAS.sif

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi
DEFAULT_SIMU_ROOT="/mnt/scratch/group/ztang2/pli297/simulation"
SIMU_ROOT="${SIMU_ROOT:-${PWD:-${DEFAULT_SIMU_ROOT}}}"
SIF="${SIF:-${SIMU_ROOT}/PALMmGWAS.sif}"
cd "${SIMU_ROOT}"

mkdir -p "${SIMU_ROOT}/logs" "${SIMU_ROOT}/input" "${SIMU_ROOT}/output"

echo "generate_input_simu: SIMU_ROOT=${SIMU_ROOT}"
echo "generate_input_simu: SIF=${SIF}"
echo "generate_input_simu: generating reproducible simulation input"

apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  "${SIMU_ROOT}/generate_simulation_data.R"

echo "generate_input_simu: finished"
