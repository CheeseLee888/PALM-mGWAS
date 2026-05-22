#!/bin/bash
#SBATCH --job-name=step4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/step4_%A_%a.out
#SBATCH --error=logs/step4_%A_%a.err

set -euo pipefail
WORK="${WORK:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "${WORK}"

SIF="${SIF:-${WORK}/PALMGWAS.sif}"
INPUT_PREFIX="${INPUT_PREFIX:-${WORK}/output/meta/step3_meta}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-${WORK}/output/plot/step4}"
FEATURE="${FEATURE:-g_Acinetobacter}"
SNP="${SNP:-chr4:1682869:G:C}"
PCUT1="${PCUT1:-5e-8}"
PCUT2="${PCUT2:-5e-8}"
PLOT_MIN_P="${PLOT_MIN_P:-1e-10}"
PLOT_WIDTH="${PLOT_WIDTH:-NA}"
PLOT_HEIGHT="${PLOT_HEIGHT:-NA}"

echo "step4: work = ${WORK}"
echo "step4: sif = ${SIF}"
echo "step4: input prefix = ${INPUT_PREFIX}"
echo "step4: output prefix = ${OUTPUT_PREFIX}"
echo "step4: feature = ${FEATURE}"
echo "step4: snp = ${SNP}"

echo "step4: running reporting mode 1 (feature + snp)"
apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}" \
  step4_reporting.R \
  --inputPrefix="${INPUT_PREFIX}" \
  --outputPrefix="${OUTPUT_PREFIX}" \
  --feature="${FEATURE}" \
  --snp="${SNP}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "step4: running reporting mode 2 (feature only)"
apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}" \
  step4_reporting.R \
  --inputPrefix="${INPUT_PREFIX}" \
  --outputPrefix="${OUTPUT_PREFIX}" \
  --feature="${FEATURE}" \
  --plotMinP="${PLOT_MIN_P}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "step4: running reporting mode 3 (snp only)"
apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}" \
  step4_reporting.R \
  --inputPrefix="${INPUT_PREFIX}" \
  --outputPrefix="${OUTPUT_PREFIX}" \
  --snp="${SNP}" \
  --pCut="${PCUT1}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "step4: running reporting mode 4 (combined)"
apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}" \
  step4_reporting.R \
  --inputPrefix="${INPUT_PREFIX}" \
  --outputPrefix="${OUTPUT_PREFIX}" \
  --pCut="${PCUT2}" \
  --plotMinP="${PLOT_MIN_P}" \
  --width="${PLOT_WIDTH}" \
  --height="${PLOT_HEIGHT}"

echo "step4: finished"
