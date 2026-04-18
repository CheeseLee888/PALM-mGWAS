#!/bin/bash
#SBATCH --job-name=simu_step4
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/simu_step4_%A_%a.out
#SBATCH --error=logs/simu_step4_%A_%a.err

set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_SIMU_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORK="${WORK:-$(cd "${DEFAULT_SIMU_ROOT}/.." && pwd)}"
SIMU_ROOT="${SIMU_ROOT:-${DEFAULT_SIMU_ROOT}}"
SIF="${SIF:-${SIMU_ROOT}/PALMmGWAS.sif}"
META_DIR="${META_DIR:-/mnt/scratch/group/ztang2/pli297/simulation/output/meta}"
PLOT_DIR="${PLOT_DIR:-/mnt/scratch/group/ztang2/pli297/simulation/output/plot}"
STEP4_PATTERN="${STEP4_PATTERN:-${VIS_PATTERN:-step3_meta_.*\\.txt$}}"
PHENO="${PHENO:-g_Acinetobacter}"
SNP="${SNP:-chr4:1682869:G:C}"
PCUT1="${PCUT1:-5e-8}"
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
  echo "step4: missing meta directory ${META_DIR}" >&2
  exit 1
fi
if [[ ! -f "${SIF}" ]]; then
  echo "step4: missing sif image ${SIF}" >&2
  exit 1
fi

META_FILE_COUNT=$(
  find "${META_DIR}" -maxdepth 1 -type f | grep -E "/$(printf '%s' "${STEP4_PATTERN}")" | wc -l | tr -d ' '
)
if [[ "${META_FILE_COUNT}" -eq 0 ]]; then
  echo "step4: no meta result files found in ${META_DIR} matching ${STEP4_PATTERN}" >&2
  exit 1
fi

echo "step4: work = ${WORK}"
echo "step4: simulation root = ${SIMU_ROOT}"
echo "step4: sif = ${SIF}"
echo "step4: meta dir = ${META_DIR}"
echo "step4: plot dir = ${PLOT_DIR}"
echo "step4: pattern = ${STEP4_PATTERN}"
echo "step4: pheno = ${PHENO}"
echo "step4: snp = ${SNP}"

echo "step4: running visualization mode 1 (pheno + snp)"
apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  step4_visualization.R \
  --metaDir="${META_DIR}" \
  --pattern="${STEP4_PATTERN}" \
  --plotDir="${PLOT_DIR}" \
  --pheno="${PHENO}" \
  --snp="${SNP}" \
  --showMeta="${SHOW_META}" \
  --showHet="${SHOW_HET}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "step4: running visualization mode 2 (pheno only)"
apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  step4_visualization.R \
  --metaDir="${META_DIR}" \
  --pattern="${STEP4_PATTERN}" \
  --plotDir="${PLOT_DIR}" \
  --pheno="${PHENO}" \
  --plotMinP="${PLOT_MIN_P}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "step4: running visualization mode 3 (snp only)"
apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  step4_visualization.R \
  --metaDir="${META_DIR}" \
  --pattern="${STEP4_PATTERN}" \
  --plotDir="${PLOT_DIR}" \
  --snp="${SNP}" \
  --pCut="${PCUT1}" \
  --showMeta="${SHOW_META}" \
  --showHet="${SHOW_HET}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "step4: running visualization mode 4 (combined)"
apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  step4_visualization.R \
  --metaDir="${META_DIR}" \
  --pattern="${STEP4_PATTERN}" \
  --plotDir="${PLOT_DIR}" \
  --pCut="${PCUT2}" \
  --plotMinP="${PLOT_MIN_P}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "step4: finished"
