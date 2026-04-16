#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################

# step2.2: merge per-chromosome Step2.1 outputs into per-feature all-chromosome files
input_prefix="${step2MergeInputPrefix:-${palm1_step2_prefix}}"
output_prefix="${step2MergeOutputPrefix:-${step2InputPrefix:-NULL}}"

if [[ -z "${input_prefix}" || "${input_prefix}" == "NULL" ]]; then
  echo "step2.2 merge: step2MergeInputPrefix (or palm1_step2_prefix) is required." >&2
  exit 1
fi

echo "Start: Perform step2.2 merge."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step2_2_merge.R \
    --inputPrefix=${input_prefix} \
    --outputPrefix=${output_prefix}

echo "Finish: Perform step2.2 merge."
