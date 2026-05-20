#!/bin/bash
#SBATCH --job-name=step2_2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/step2_2_%j.out
#SBATCH --error=logs/step2_2_%j.err
# Step2.2 in the fixed simulation workflow:
# one job per study, correcting all-chromosome Step2.1 files.

set -euo pipefail
WORK="${WORK:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIF="${SIF:-${WORK}/PALMGWAS.sif}"

cd "${WORK}"

OUTPUT_DIR="${WORK}/output/${STUDY}"
INPUT_PREFIX="${OUTPUT_DIR}/step2"
OVERWRITE_OUTPUT="${OVERWRITE_OUTPUT:-TRUE}"
STEP2_CORRECT_MODE="${STEP2_CORRECT_MODE:-median}"
CHROM="NULL"

echo "run_step2_2: WORK = ${WORK}"
echo "run_step2_2: study = ${STUDY}"
echo "run_step2_2: mode = all chromosomes"
echo "run_step2_2: chrom = ${CHROM}"
echo "run_step2_2: overwriteOutput = ${OVERWRITE_OUTPUT}"
echo "run_step2_2: correct = ${STEP2_CORRECT_MODE}"

echo "run_step2_2: inputPrefix=${INPUT_PREFIX}"

apptainer exec --cleanenv --bind "${WORK}:${WORK}" "${SIF}" \
  step2_2_correction.R \
  --inputPrefix="${INPUT_PREFIX}" \
  --chrom="${CHROM}" \
  --overwriteOutput="${OVERWRITE_OUTPUT}" \
  --correct="${STEP2_CORRECT_MODE}"
