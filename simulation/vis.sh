#!/bin/bash
set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${WORK:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
SIMU_ROOT="${SIMU_ROOT:-${SCRIPT_DIR}}"
SIF="${SIF:-${SIMU_ROOT}/PALMmGWAS.sif}"
META_DIR="${META_DIR:-/mnt/scratch/group/ztang2/pli297/testsimu/output/meta}"
PLOT_DIR="${PLOT_DIR:-/mnt/scratch/group/ztang2/pli297/testsimu/output/plot}"
VIS_PATTERN="${VIS_PATTERN:-step3_meta_.*\\.txt$}"
PHENO="${PHENO:-g_Acinetobacter}"
SNP="${SNP:-chr4:1682869:G:C}"
PCUT1="${PCUT1:-1e-20}"
PCUT2="${PCUT2:-5e-8}"
SHOW_META="${SHOW_META:-TRUE}"
SHOW_HET="${SHOW_HET:-TRUE}"
PLOT_MIN_P="${PLOT_MIN_P:-1e-10}"
PLOT_WIDTH="${PLOT_WIDTH:-NA}"
PLOT_HEIGHT="${PLOT_HEIGHT:-NA}"

WORK="$(cd "${WORK}" && pwd)"
SIMU_ROOT="$(cd "${SIMU_ROOT}" && pwd)"
SIF="$(cd "$(dirname "${SIF}")" && pwd)/$(basename "${SIF}")"
META_DIR="$(mkdir -p "${META_DIR}" && cd "${META_DIR}" && pwd)"
PLOT_DIR="$(mkdir -p "${PLOT_DIR}" && cd "${PLOT_DIR}" && pwd)"

if [[ ! -d "${META_DIR}" ]]; then
  echo "vis: missing meta directory ${META_DIR}" >&2
  exit 1
fi
if [[ ! -f "${SIF}" ]]; then
  echo "vis: missing sif image ${SIF}" >&2
  exit 1
fi

META_FILE_COUNT=$(
  find "${META_DIR}" -maxdepth 1 -type f | grep -E "/$(printf '%s' "${VIS_PATTERN}")" | wc -l | tr -d ' '
)
if [[ "${META_FILE_COUNT}" -eq 0 ]]; then
  echo "vis: no meta result files found in ${META_DIR} matching ${VIS_PATTERN}" >&2
  exit 1
fi

echo "vis: work = ${WORK}"
echo "vis: simulation root = ${SIMU_ROOT}"
echo "vis: sif = ${SIF}"
echo "vis: meta dir = ${META_DIR}"
echo "vis: plot dir = ${PLOT_DIR}"
echo "vis: pattern = ${VIS_PATTERN}"
echo "vis: pheno = ${PHENO}"
echo "vis: snp = ${SNP}"

echo "vis: running visualization mode 1 (pheno + snp)"
apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  Visualization.R \
  --metaDir="${META_DIR}" \
  --pattern="${VIS_PATTERN}" \
  --plotDir="${PLOT_DIR}" \
  --pheno="${PHENO}" \
  --snp="${SNP}" \
  --showMeta="${SHOW_META}" \
  --showHet="${SHOW_HET}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "vis: running visualization mode 2 (pheno only)"
apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  Visualization.R \
  --metaDir="${META_DIR}" \
  --pattern="${VIS_PATTERN}" \
  --plotDir="${PLOT_DIR}" \
  --pheno="${PHENO}" \
  --plotMinP="${PLOT_MIN_P}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "vis: running visualization mode 3 (snp only)"
apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  Visualization.R \
  --metaDir="${META_DIR}" \
  --pattern="${VIS_PATTERN}" \
  --plotDir="${PLOT_DIR}" \
  --snp="${SNP}" \
  --pCut="${PCUT1}" \
  --showMeta="${SHOW_META}" \
  --showHet="${SHOW_HET}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "vis: running visualization mode 4 (combined)"
apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  Visualization.R \
  --metaDir="${META_DIR}" \
  --pattern="${VIS_PATTERN}" \
  --plotDir="${PLOT_DIR}" \
  --pCut="${PCUT2}" \
  --plotMinP="${PLOT_MIN_P}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "vis: finished"
