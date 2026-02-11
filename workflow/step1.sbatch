#!/bin/bash
#SBATCH --job-name=palm-step1
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/scratch/group/ztang2/pli297/PALMmbQTL/logs/%x_%A_%a.out
#SBATCH --error=/mnt/scratch/group/ztang2/pli297/PALMmbQTL/logs/%x_%A_%a.err
set -euo pipefail
source "$WORK/config.sh"


################################# workflow below (do not modify) #################################
if [[ "${PALMmethod}" != 1 && "${PALMmethod}" != 2 ]]; then
    echo "PALMmethod must be specified as 1 (PALM) or 2 (PALM-mbQTL)."
    exit 1
fi

if [[ "${PALMmethod}" == 1 ]]; then
    echo "Running PALM method..."
else
    echo "Running PALM-mbQTL method..."
fi

mkdir -p "${outputFolder}"

# step0: check input files
echo "Start: Check input files."
apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
    step0_checkInput.R \
    --abdFile=${abdFile} \
    --covFile=${covFile} \
    --genoFile=${genoFile}
echo "Finish: Check input files."    

# # step0: generate GRM from genotype data (only for PALM-mbQTL)
# if [[ "${PALMmethod}" == 2 ]]; then
#     echo "Start: Generate GRM."
#     if [[ -f "${grmFile}" ]]; then
#         echo "GRM already exists at: ${grmFile}"
#         echo "Skip generating GRM and reuse the existing GRM."
#     else
#         apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
#             step0_generateGRM.R \
#             --genoFile=${genoFile} \
#             --grmFile=${grmFile}
#     fi
#     echo "Finish: Generate GRM."
# fi

# step1: fit null model for all phenotypes
echo "Start: Fit null model."
if [[ "${PALMmethod}" == 1 ]]; then
    apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
        step1_palm.R \
        --abdFile=${abdFile} \
        --covFile=${covFile} \
        --outputPrefix=${palm1_step1_prefix}
else
    # Progress file: records "which batch_idx to run next"
    RESET_STEP1=1
    PROGRESS="${palm2_step1_prefix}_progress.txt"

    if [ "${RESET_STEP1}" -eq 1 ]; then
        echo "[step1] RESET enabled: removing progress and restarting from 1"
        rm -f "$PROGRESS"
        echo 1 > "$PROGRESS"
    else
        if [ ! -f "$PROGRESS" ]; then
            echo 1 > "$PROGRESS"
        fi
    fi

    attempt=1
    set +e   # Allow commands to fail inside the loop
    
    while true; do
        batch_idx=$(cat "$PROGRESS")

        echo "=== attempt=${attempt} resume batch_idx=${batch_idx} ==="

        apptainer exec --cleanenv --bind $WORK:$WORK $SIF \
            step1_fitNULL.R \
            --grmFile=${grmFile} \
            --abdFile=${abdFile} \
            --covFile=${covFile} \
            --outputPrefix=${palm2_step1_prefix} \
            --batch_idx=${batch_idx}

        rc=$?

        # exit after all phenotypes are processed
        if [ $rc -eq 0 ]; then
            echo "All phenotypes processed. Stop."
            rm -f "$PROGRESS"
            break
        fi

        if [ $rc -eq 2 ]; then
            echo "Fatal preflight error detected. Abort." >&2
            exit 1
        fi

        echo "Run crashed. Restart; will resume from progress=$(cat "$PROGRESS")" >&2
        attempt=$((attempt + 1))
    done
    set -e  # Re-enable exit on error
fi
echo "Finish: Fit null model."