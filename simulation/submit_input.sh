#!/bin/bash
set -euo pipefail

# Submit only the simulation input generation job.

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
cd "${SIMU_ROOT}"

jid=$(sbatch --parsable generate_input_simu.sbatch)
echo "submit_input: submitted input generation job ${jid}"
