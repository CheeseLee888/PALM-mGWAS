#!/bin/bash
set -euo pipefail

# One-command Slurm submission for simulation analysis only:
# Step0 -> Step1 -> Step2 array.
# Step2 is submitted by a lightweight bridge job after Step1 finishes.

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
STUDY="${STUDY:-study1}"
FEATURE_BLOCK="${FEATURE_BLOCK:-10}"
CHROMS="${CHROMS:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"

cd "${SIMU_ROOT}"

echo "submit_all: simulation root = ${SIMU_ROOT}"
echo "submit_all: study = ${STUDY}"
echo "submit_all: feature block = ${FEATURE_BLOCK}"
echo "submit_all: chrom list = ${CHROMS}"

jid0=$(sbatch --parsable step0_simu.sbatch)
echo "submit_all: submitted Step0 job ${jid0}"

jid1=$(sbatch --parsable --dependency=afterok:${jid0} step1_simu.sbatch)
echo "submit_all: submitted Step1 job ${jid1}"

jid2=$(sbatch --parsable \
  --dependency=afterok:${jid1} \
  --export=ALL,SIMU_ROOT=${SIMU_ROOT},STUDY=${STUDY},FEATURE_BLOCK=${FEATURE_BLOCK} \
  submit_step2_bridge.sbatch)
echo "submit_all: submitted Step2 bridge job ${jid2}"
