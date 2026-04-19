#!/bin/bash
#SBATCH --job-name=simu_step2_2
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/simu_step2_2_%A_%a.out
#SBATCH --error=logs/simu_step2_2_%A_%a.err

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

OUTPUT_DIR="${SIMU_ROOT}/output/${STUDY}"

echo "run_step2_2: simulation root = ${SIMU_ROOT}"
echo "run_step2_2: study = ${STUDY}"

apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  step2_2_merge.R \
  --inputPrefix="${OUTPUT_DIR}/step2" \
  --outputPrefix="${OUTPUT_DIR}/step2"

echo "run_step2_2: finished"
