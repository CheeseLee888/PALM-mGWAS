#!/usr/bin/env bash
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################

# merge: combine per-chromosome Step2.1 outputs into per-feature all-chromosome files
input_prefix="${mergeInputPrefix:-NULL}"
output_prefix="${mergeOutputPrefix:-${step2InputPrefix:-NULL}}"

if [[ -z "${input_prefix}" || "${input_prefix}" == "NULL" ]]; then
  input_prefix="${SummaryPrefix}"
fi

echo "Start: Perform merge."
pixi run --manifest-path=${WORK}/pixi.toml Rscript ${WORK}/extdata/step_merge.R \
    --inputPrefix=${input_prefix} \
    --featureColList=${mergeFeatureColList:-NULL} \
    --outputPrefix=${output_prefix}

echo "Finish: Perform merge."
