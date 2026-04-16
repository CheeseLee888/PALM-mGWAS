#!/bin/bash
set -euo pipefail

if [[ -n "${ENV_FILE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

SIMU_ROOT="${SIMU_ROOT:-/mnt/scratch/group/ztang2/pli297/simulation}"
cd "${SIMU_ROOT}"

echo "merge_step2_meta: running merge step"
bash "${SIMU_ROOT}/merge_step2_outputs.sh"

echo "merge_step2_meta: running meta step"
bash "${SIMU_ROOT}/run_step3_meta.sh"

echo "merge_step2_meta: finished merge + meta"
