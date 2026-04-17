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

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
if [[ -n "${STUDIES:-}" ]]; then
  STUDIES_CSV="${STUDIES}"
elif [[ -n "${STUDY:-}" ]]; then
  STUDIES_CSV="${STUDY}"
else
  STUDIES_CSV="study1,study2,study3"
fi
META_DIR="${META_DIR:-/mnt/scratch/group/ztang2/pli297/simulation/output/meta}"
SIF="${SIF:-${SIMU_ROOT}/PALMmGWAS.sif}"

cd "${SIMU_ROOT}"
mkdir -p "${SIMU_ROOT}/logs" "${META_DIR}"

IFS=',' read -r -a STUDY_LIST <<< "${STUDIES_CSV}"
if [[ "${#STUDY_LIST[@]}" -eq 0 ]]; then
  echo "run_step2_2: no studies resolved from STUDIES/STUDY" >&2
  exit 1
fi

echo "run_step2_2: simulation root = ${SIMU_ROOT}"
echo "run_step2_2: studies = ${STUDIES_CSV}"

STUDY_FILE="${META_DIR}/study_dirs.tsv"
: > "${STUDY_FILE}"

for study_name in "${STUDY_LIST[@]}"; do
  OUTPUT_DIR="${SIMU_ROOT}/output/${study_name}"
  if [[ ! -d "${OUTPUT_DIR}" ]]; then
    echo "run_step2_2: missing output directory ${OUTPUT_DIR}" >&2
    exit 1
  fi

  apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
    step2_2_merge.R \
    --inputPrefix="${OUTPUT_DIR}/step2" \
    --outputPrefix="${OUTPUT_DIR}/step2"

  printf '%s\t%s\n' "${study_name}" "${OUTPUT_DIR}" >> "${STUDY_FILE}"
done

echo "run_step2_2: wrote ${STUDY_FILE}"
