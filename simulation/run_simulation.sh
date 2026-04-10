#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PIXIMANIFEST="${ROOT}/pixi.toml"

INPUT_DIR="${ROOT}/simulation/input/study1"
OUTPUT_DIR="${ROOT}/simulation/output/study1"
PLOT_DIR="${ROOT}/simulation/output/plot"

ABD_FILE="${INPUT_DIR}/abd.txt"
COV_FILE="${INPUT_DIR}/cov.txt"
GENO_PREFIX="${INPUT_DIR}/geno"

ABD_ALIGNED_FILE="${INPUT_DIR}/abd_aligned.txt"
COV_ALIGNED_FILE="${INPUT_DIR}/cov_aligned.txt"
NULLMODEL_FILE="${OUTPUT_DIR}/step1_allpheno.rda"
STEP2_PREFIX="${OUTPUT_DIR}/step2_allchr"
FEATURE_INFO_FILE="${OUTPUT_DIR}/info_feature.txt"
SNP_INFO_FILE="${OUTPUT_DIR}/info_snp.txt"
SEQDEPTH_INFO_FILE="${OUTPUT_DIR}/info_seqdepth.txt"

echo "Generating reproducible simulation input..."
pixi run --manifest-path="${PIXIMANIFEST}" Rscript "${ROOT}/simulation/generate_simulation_data.R"

mkdir -p "${OUTPUT_DIR}" "${PLOT_DIR}"

echo "Step0: aligning abundance/covariate/genotype inputs..."
pixi run --manifest-path="${PIXIMANIFEST}" Rscript "${ROOT}/extdata/step0_checkInput.R" \
  --abdFile="${ABD_FILE}" \
  --covFile="${COV_FILE}" \
  --abdAlignedFile="${ABD_ALIGNED_FILE}" \
  --covAlignedFile="${COV_ALIGNED_FILE}" \
  --covarColList=AGE,SEX,PC01,PC02,PC03,PC04,PC05 \
  --depthCol=NULL \
  --depth.filter=0 \
  --genoFile="${GENO_PREFIX}" \
  --outputSeqDepthFile="${SEQDEPTH_INFO_FILE}"

echo "Step1: fitting PALM null model..."
pixi run --manifest-path="${PIXIMANIFEST}" Rscript "${ROOT}/extdata/step1_null.R" \
  --abdFile="${ABD_ALIGNED_FILE}" \
  --covFile="${COV_ALIGNED_FILE}" \
  --covarColList=AGE,SEX,PC01,PC02,PC03,PC04,PC05 \
  --depthCol=NULL \
  --prev.filter=0.1 \
  --outputFeatureFile="${FEATURE_INFO_FILE}" \
  --NULLmodelFile="${NULLMODEL_FILE}"

echo "Step2: generating per-feature summary statistics..."
pixi run --manifest-path="${PIXIMANIFEST}" Rscript "${ROOT}/extdata/step2_summary.R" \
  --genoFile="${GENO_PREFIX}" \
  --NULLmodelFile="${NULLMODEL_FILE}" \
  --PALMOutputFile="${STEP2_PREFIX}" \
  --chrom=NULL \
  --minMAF=0.05 \
  --minMAC=5 \
  --outputSnpFile="${SNP_INFO_FILE}" \
  --correct=NULL \
  --useCluster=FALSE

echo "Visualization: writing Manhattan/QQ plot for g_Acinetobacter..."
pixi run --manifest-path="${PIXIMANIFEST}" Rscript "${ROOT}/extdata/Visualization.R" \
  --metaDir="${OUTPUT_DIR}" \
  --pattern='step2_allchr_.*\.txt$' \
  --plotDir="${PLOT_DIR}" \
  --pheno=g_Acinetobacter \
  --plotMinP=1e-12

echo "Simulation workflow finished."
echo "Key outputs:"
echo "  ${OUTPUT_DIR}/step2_allchr_g_Acinetobacter.txt"
echo "  ${PLOT_DIR}/manhattan_g_Acinetobacter.png"
