#!/bin/bash
#SBATCH --job-name=simu_step2_3
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/simu_step2_3_%A_%a.out
#SBATCH --error=logs/simu_step2_3_%A_%a.err

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
if [[ -n "${STUDIES:-}" ]]; then
  STUDIES_CSV="${STUDIES}"
elif [[ -n "${STUDY:-}" ]]; then
  STUDIES_CSV="${STUDY}"
else
  STUDIES_CSV="study1,study2,study3"
fi
SIF="${SIF:-${SIMU_ROOT}/PALMmGWAS.sif}"
OVERWRITE_OUTPUT="${OVERWRITE_OUTPUT:-TRUE}"

cd "${SIMU_ROOT}"
mkdir -p "${SIMU_ROOT}/logs"

IFS=',' read -r -a STUDY_LIST <<< "${STUDIES_CSV}"
if [[ "${#STUDY_LIST[@]}" -eq 0 ]]; then
  echo "run_step2_3_correction: no studies resolved from STUDIES/STUDY" >&2
  exit 1
fi

echo "run_step2_3_correction: simulation root = ${SIMU_ROOT}"
echo "run_step2_3_correction: studies = ${STUDIES_CSV}"
echo "run_step2_3_correction: overwriteOutput = ${OVERWRITE_OUTPUT}"

for study_name in "${STUDY_LIST[@]}"; do
  OUTPUT_DIR="${SIMU_ROOT}/output/${study_name}"
  INPUT_PREFIX="${OUTPUT_DIR}/step2_allchr"

  if ! find "${OUTPUT_DIR}" -maxdepth 1 -type f -name 'step2_allchr_*.txt' | grep -q .; then
    echo "run_step2_3_correction: no merged Step2.2 files found in ${OUTPUT_DIR}" >&2
    exit 1
  fi

  echo "run_step2_3_correction: study=${study_name}"
  echo "run_step2_3_correction: inputPrefix=${INPUT_PREFIX}"

  apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
    step2_3_correct.R \
    --inputPrefix="${INPUT_PREFIX}" \
    --overwriteOutput="${OVERWRITE_OUTPUT}"
done
