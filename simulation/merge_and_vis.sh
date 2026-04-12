#!/bin/bash
set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
SIF="${SIF:-${SIMU_ROOT}/PALMmGWAS.sif}"
if [[ -n "${STUDIES:-}" ]]; then
  STUDIES_CSV="${STUDIES}"
elif [[ -n "${STUDY:-}" ]]; then
  STUDIES_CSV="${STUDY}"
else
  STUDIES_CSV="study1,study2,study3"
fi
PLOT_DIR="${PLOT_DIR:-${SIMU_ROOT}/output/plot}"
META_DIR="${META_DIR:-${SIMU_ROOT}/output/meta}"
META_PREFIX="${META_PREFIX:-step3_meta_}"
VIS_PATTERN="${VIS_PATTERN:-step3_meta_.*\\.txt$}"
PHENO="${PHENO:-g_Acinetobacter}"
SNP="${SNP:-chr4:1682869:G:C}"
PCUT1="${PCUT1:-5e-8}"
PCUT2="${PCUT2:-5e-8}"
SHOW_META="${SHOW_META:-TRUE}"
SHOW_HET="${SHOW_HET:-TRUE}"
PLOT_MIN_P="${PLOT_MIN_P:-1e-10}"
PLOT_WIDTH="${PLOT_WIDTH:-NA}"
PLOT_HEIGHT="${PLOT_HEIGHT:-NA}"

cd "${SIMU_ROOT}"
mkdir -p "${SIMU_ROOT}/logs" "${PLOT_DIR}" "${META_DIR}"

IFS=',' read -r -a STUDY_LIST <<< "${STUDIES_CSV}"
if [[ "${#STUDY_LIST[@]}" -eq 0 ]]; then
  echo "merge_and_vis: no studies resolved from STUDIES/STUDY" >&2
  exit 1
fi

echo "merge_and_vis: simulation root = ${SIMU_ROOT}"
echo "merge_and_vis: studies = ${STUDIES_CSV}"
echo "merge_and_vis: meta dir = ${META_DIR}"
echo "merge_and_vis: plot dir = ${PLOT_DIR}"

STUDY_FILE="${META_DIR}/study_dirs.tsv"
: > "${STUDY_FILE}"

for study_name in "${STUDY_LIST[@]}"; do
  OUTPUT_DIR="${SIMU_ROOT}/output/${study_name}"
  if [[ ! -d "${OUTPUT_DIR}" ]]; then
    echo "merge_and_vis: missing output directory ${OUTPUT_DIR}" >&2
    exit 1
  fi

  mapfile -t STEP2_FILES < <(find "${OUTPUT_DIR}" -maxdepth 1 -type f -name 'step2_chr*_*.txt' | sort -V)
  if [[ "${#STEP2_FILES[@]}" -eq 0 ]]; then
    echo "merge_and_vis: no per-chromosome Step2 files found in ${OUTPUT_DIR}" >&2
    exit 1
  fi

  echo "merge_and_vis: ${study_name} found ${#STEP2_FILES[@]} per-chromosome Step2 file(s)"

  mapfile -t PHENOS < <(
    printf '%s\n' "${STEP2_FILES[@]}" |
      xargs -n1 basename |
      sed -E 's/^step2_chr[^_]+_//; s/[.]txt$//' |
      sort -u
  )

  if [[ "${#PHENOS[@]}" -eq 0 ]]; then
    echo "merge_and_vis: could not resolve phenotype names from Step2 filenames in ${OUTPUT_DIR}" >&2
    exit 1
  fi

  for pheno in "${PHENOS[@]}"; do
    out_file="${OUTPUT_DIR}/step2_allchr_${pheno}.txt"
    rm -f "${out_file}"

    mapfile -t pheno_files < <(find "${OUTPUT_DIR}" -maxdepth 1 -type f -name "step2_chr*_${pheno}.txt" | sort -V)
    if [[ "${#pheno_files[@]}" -eq 0 ]]; then
      continue
    fi

    awk 'FNR == 1 && NR != 1 { next } { print }' "${pheno_files[@]}" > "${out_file}"
    echo "merge_and_vis: ${study_name} merged ${#pheno_files[@]} chromosome file(s) -> ${out_file}"
  done

  printf '%s\t%s\n' "${study_name}" "${OUTPUT_DIR}" >> "${STUDY_FILE}"
done

echo "merge_and_vis: running meta-analysis"
apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  step3_meta.R \
  --studyFile="${STUDY_FILE}" \
  --pattern='step2_allchr_.*\.txt$' \
  --metaDir="${META_DIR}" \
  --metaPrefix="${META_PREFIX}"

echo "merge_and_vis: running visualization mode 1 (pheno + snp)"
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

echo "merge_and_vis: running visualization mode 2 (pheno only)"
apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  Visualization.R \
  --metaDir="${META_DIR}" \
  --pattern="${VIS_PATTERN}" \
  --plotDir="${PLOT_DIR}" \
  --pheno="${PHENO}" \
  --plotMinP="${PLOT_MIN_P}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "merge_and_vis: running visualization mode 3 (snp only)"
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

echo "merge_and_vis: running visualization mode 4 (combined)"
apptainer exec --cleanenv --bind "${SIMU_ROOT}:${SIMU_ROOT}" "${SIF}" \
  Visualization.R \
  --metaDir="${META_DIR}" \
  --pattern="${VIS_PATTERN}" \
  --plotDir="${PLOT_DIR}" \
  --pCut="${PCUT2}" \
  --plotMinP="${PLOT_MIN_P}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "merge_and_vis: finished"
