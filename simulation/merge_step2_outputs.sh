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
if [[ -n "${STUDIES:-}" ]]; then
  STUDIES_CSV="${STUDIES}"
elif [[ -n "${STUDY:-}" ]]; then
  STUDIES_CSV="${STUDY}"
else
  STUDIES_CSV="study1,study2,study3"
fi
PLOT_DIR="${PLOT_DIR:-${SIMU_ROOT}/output/plot}"
META_DIR="${META_DIR:-/mnt/scratch/group/ztang2/pli297/simulation/output/meta}"

cd "${SIMU_ROOT}"
mkdir -p "${SIMU_ROOT}/logs" "${PLOT_DIR}" "${META_DIR}"

IFS=',' read -r -a STUDY_LIST <<< "${STUDIES_CSV}"
if [[ "${#STUDY_LIST[@]}" -eq 0 ]]; then
  echo "merge_step2_outputs: no studies resolved from STUDIES/STUDY" >&2
  exit 1
fi

echo "merge_step2_outputs: simulation root = ${SIMU_ROOT}"
echo "merge_step2_outputs: studies = ${STUDIES_CSV}"

STUDY_FILE="${META_DIR}/study_dirs.tsv"
: > "${STUDY_FILE}"

for study_name in "${STUDY_LIST[@]}"; do
  OUTPUT_DIR="${SIMU_ROOT}/output/${study_name}"
  if [[ ! -d "${OUTPUT_DIR}" ]]; then
    echo "merge_step2_outputs: missing output directory ${OUTPUT_DIR}" >&2
    exit 1
  fi

  mapfile -t STEP2_ALLCHR_FILES < <(find "${OUTPUT_DIR}" -maxdepth 1 -type f -name 'step2_allchr_*.txt' | sort -V)
  if [[ "${#STEP2_ALLCHR_FILES[@]}" -gt 0 ]]; then
    echo "merge_step2_outputs: ${study_name} found ${#STEP2_ALLCHR_FILES[@]} all-chromosome Step2.1 file(s)"
  else
    mapfile -t STEP2_FILES < <(find "${OUTPUT_DIR}" -maxdepth 1 -type f -name 'step2_chr*_*.txt' | sort -V)
    if [[ "${#STEP2_FILES[@]}" -eq 0 ]]; then
      echo "merge_step2_outputs: no Step2.1 files found in ${OUTPUT_DIR}" >&2
      exit 1
    fi

    echo "merge_step2_outputs: ${study_name} found ${#STEP2_FILES[@]} per-chromosome Step2.1 file(s)"

    mapfile -t PHENOS < <(
      printf '%s\n' "${STEP2_FILES[@]}" |
        xargs -n1 basename |
        sed -E 's/^step2_chr[^_]+_//; s/[.]txt$//' |
        sort -u
    )

    if [[ "${#PHENOS[@]}" -eq 0 ]]; then
      echo "merge_step2_outputs: could not resolve phenotype names from Step2.1 filenames in ${OUTPUT_DIR}" >&2
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
      echo "merge_step2_outputs: ${study_name} merged ${#pheno_files[@]} chromosome file(s) -> ${out_file}"
    done
  fi

  printf '%s\t%s\n' "${study_name}" "${OUTPUT_DIR}" >> "${STUDY_FILE}"
done

echo "merge_step2_outputs: wrote ${STUDY_FILE}"
