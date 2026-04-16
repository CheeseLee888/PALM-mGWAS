#!/bin/bash
set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
SIF="${SIF:-${SIMU_ROOT}/PALMmGWAS.sif}"
META_DIR="${META_DIR:-/mnt/scratch/group/ztang2/pli297/simulation/output/meta}"
META_PREFIX="${META_PREFIX:-step3_meta_}"
STUDY_FILE="${STUDY_FILE:-${META_DIR}/study_dirs.tsv}"
FEATURE="${FEATURE:-}"

cd "${SIMU_ROOT}"
mkdir -p "${SIMU_ROOT}/logs" "${META_DIR}"

if [[ ! -f "${STUDY_FILE}" ]]; then
  echo "run_step3_meta: missing study dir file ${STUDY_FILE}" >&2
  echo "run_step3_meta: run merge_step2_outputs.sh first." >&2
  exit 1
fi

echo "run_step3_meta: simulation root = ${SIMU_ROOT}"
echo "run_step3_meta: study dir file = ${STUDY_FILE}"
echo "run_step3_meta: meta dir = ${META_DIR}"
if [[ -n "${FEATURE}" ]]; then
  echo "run_step3_meta: feature = ${FEATURE}"
fi

CMD=(
  apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}"
  step3_meta.R
  --studyDirFile="${STUDY_FILE}"
  --pattern='step2_allchr_.*\.txt$'
  --metaDir="${META_DIR}"
  --metaPrefix="${META_PREFIX}"
)

if [[ -n "${FEATURE}" ]]; then
  CMD+=(--features="${FEATURE}")
fi

"${CMD[@]}"

echo "run_step3_meta: finished"
