# Simulation Dataset

This directory provides a reproducible single-study simulation intended to show a larger sample size than `example/` and give users a more realistic sense of `step2` runtime and memory usage.

The working directory is:

- `/mnt/scratch/group/ztang2/pli297/simulation`

That directory directly contains:

- `provided/`
- `input/`
- `output/`
- `PALMmGWAS.sif`
- the submission scripts in this folder

Core design:

- The sample size is fixed at `502`, defined in `simulation/provided/sample_config.tsv`
- The SNP list and all `CHR/POS/A1/A2/AF` values are taken from `simulation/provided/snp_reference.tsv`
- Feature names and covariate names are stored in `simulation/provided/feature_names.txt` and `simulation/provided/covariate_names.txt`
- The generated `abd` and `cov` tables follow the same overall column structure as `example/input/study1`
- The random seed is fixed at `20260409`, so any user running the same script will generate identical inputs and outputs
- A strong genetic signal is planted in `g_Acinetobacter` so the Manhattan plot contains clear high peaks

## Validation Status

The simulation workflow has been run successfully through `Step0 -> Step2` on the cluster layout described below.

In particular, `step2` has been verified under Slurm with a `chrom x feature block` decomposition:

- chromosome subsetting is controlled by `--chrom`
- feature subsetting is controlled by `--featureColList`
- parallelism is handled at the scheduler level by submitting multiple jobs, rather than inside one R process

This means the current implementation supports user-managed `chrom x feature` parallel execution for `step2`.

## One-command run

From the simulation working directory, run:

```bash
bash run_simulation.sh
```

The script will automatically:

1. Generate `input/study1/abd.txt`
2. Generate `input/study1/cov.txt`
3. Generate `input/study1/geno.{bed,bim,fam}`
4. Run `Step0 -> Step1 -> Step2`
5. Write the Manhattan/QQ plot for `g_Acinetobacter`

## Slurm submission order

You can split the workflow into two phases.

Phase 1: generate simulation input

```bash
cd /mnt/scratch/group/ztang2/pli297/simulation
bash submit_input.sh
```

This submits [generate_input_simu.sbatch](/Users/peterli/Desktop/PALM-mGWAS/simulation/generate_input_simu.sbatch), which creates the files under `input/study1/`.

Phase 2: run analysis

```bash
cd /mnt/scratch/group/ztang2/pli297/simulation
bash submit_all.sh
```

This submits `Step0 -> Step1 -> Step2`.

You can also run the analysis chain manually in three scheduler steps:

```bash
ENV_FILE=/path/to/your.env
jid0=$(sbatch --parsable --export=ALL,ENV_FILE=${ENV_FILE} step0_simu.sbatch)
jid1=$(sbatch --parsable --dependency=afterok:${jid0} --export=ALL,ENV_FILE=${ENV_FILE} step1_simu.sbatch)
```

The default simulation root used by these sbatch scripts is:

```bash
/mnt/scratch/group/ztang2/pli297/simulation
```

That directory is expected to contain at least `generate_simulation_data.R`, `input/`, `output/`, `provided/`, and `PALMmGWAS.sif`. You can still override `SIMU_ROOT` or `SIF` from your env file if needed.

Default simulation parameters are written directly inside the sbatch scripts. If needed, you can override them at submit time with `sbatch --export=...`.

Then compute the Step2 array size and submit the chromosome x feature-block jobs:

```bash
FEATURE_BLOCK=10
N_FEATURE=$(wc -l < output/study1/feature_list.txt)
N_BLOCK=$(( (N_FEATURE + FEATURE_BLOCK - 1) / FEATURE_BLOCK ))
N_CHROM=22
TOTAL_TASKS=$(( N_BLOCK * N_CHROM ))

sbatch \
  --dependency=afterok:${jid1} \
  --array=1-${TOTAL_TASKS} \
  --export=ALL,ENV_FILE=${ENV_FILE} \
  step2_simu_array.sbatch
```

`step1_simu.sbatch` writes `output/study1/feature_list.txt` automatically, so Step2 can use it directly.

`submit_all.sh` submits all three analysis stages in one command: Step0, Step1, and then a bridge job which computes `TOTAL_TASKS` after Step1 finishes and submits the Step2 array automatically.

## Main outputs

- `output/study1/step2_allchr_g_Acinetobacter.txt`
- `output/study1/info_snp.txt`
- `output/study1/info_feature.txt`
- `output/plot/manhattan_g_Acinetobacter.png`

## Slurm array for Step2

If you want to parallelize simulated `step2` runs across both chromosome and feature subsets, run `Step0` and `Step1` once first, then submit `step2_simu_array.sbatch`.

This has already been validated in the simulation workflow: `step2` array jobs were launched across multiple chromosomes and feature blocks, and the expected per-chromosome outputs were produced.

Choose a feature block size, for example `10` features per task:

```bash
ENV_FILE=/path/to/your.env
FEATURE_BLOCK=10
N_FEATURE=$(wc -l < output/study1/feature_list.txt)
N_BLOCK=$(( (N_FEATURE + FEATURE_BLOCK - 1) / FEATURE_BLOCK ))
N_CHROM=22
TOTAL_TASKS=$(( N_BLOCK * N_CHROM ))

sbatch --array=1-${TOTAL_TASKS} --export=ALL,ENV_FILE=${ENV_FILE} step2_simu_array.sbatch
```

You can override settings at submit time:

```bash
sbatch \
  --array=1-${TOTAL_TASKS} \
  --export=ALL,ENV_FILE=${ENV_FILE},FEATURE_BLOCK=10,STUDY=study1 \
  step2_simu_array.sbatch
```

If you need a custom chromosome list, edit the `CHROMS` default directly in the sbatch script instead of passing a comma-separated list through `sbatch --export`, because Slurm uses commas as variable separators in `--export`.

This script uses 1-based array indexing:

- task `1` is chromosome `1`, feature block `1`
- task `2` is chromosome `1`, feature block `2`
- ...
- after the last block on chromosome `1`, the next task moves to chromosome `2`

## Note

This simulation currently uses a single-study workflow, so the visualization step reads `step2` outputs directly instead of `step3` meta-analysis outputs.

The simulation working directory keeps the lightweight metadata under `provided/`. Users regenerate `input/` and `output/` locally by running the scripts above.
