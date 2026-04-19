# Simulation Dataset

This directory provides a reproducible multi-study simulation intended to show a larger sample size than `example/`, give users a more realistic sense of `step2.1` runtime and memory usage, and provide compatible inputs for later meta-analysis testing.

The working directory is:

- `/mnt/scratch/group/ztang2/pli297/simulation`

That directory directly contains:

- `provided/`
- `input/`
- `output/`
- `run/`
- `PALMmGWAS.sif`
- `generate_simulation_data.R`

Directory layout:

- `simulation/` root: submit scripts that only orchestrate and call `sbatch`
- `run/`: scripts that do the actual computation or visualization

Core design:

- The generator now writes three studies: `study1`, `study2`, and `study3`
- The baseline sample size is defined in `simulation/provided/sample_config.tsv`, and `study2` and `study3` use different sample sizes to mimic heterogeneous cohorts
- The SNP list and all `CHR/POS/A1/A2/AF` values are taken from `simulation/provided/snp_reference.tsv`
- Feature names and covariate names are stored in `simulation/provided/feature_names.txt` and `simulation/provided/covariate_names.txt`
- The three studies intentionally share the same SNP and feature backbone, but `study2` and `study3` each drop a small subset of SNPs and features so later meta-analysis can test imperfect overlap
- The generated `abd` and `cov` tables follow the same overall column structure as `example/input/study1`
- The random seed is fixed at `20260409`, so any user running the same script will generate identical inputs and outputs
- A strong genetic signal is planted in `g_Acinetobacter` so the Manhattan plot contains clear high peaks

## Validation Status

The simulation workflow has been run successfully through `Step0 -> Step2.3` on the cluster layout described below.

In particular, `step2.1` has been verified under Slurm with a `one chromosome per job` decomposition:

- chromosome subsetting is controlled by `--chrom`
- all features on that chromosome are modeled together with `--featureColList=NULL`
- parallelism is handled at the scheduler level by submitting multiple jobs, rather than inside one R process

This means the simulation workflow follows one fixed route by default: `step2.1` runs as `chrom only`, `step2.2` merges chromosome-split files into `step2_allchr_*`, `step2.3` runs with `--chrom=NULL` and overwrites those merged files in place, and `step3` reads `step2_allchr_*.txt`.

## Slurm submission order

You can split the workflow into two phases.

Phase 1: generate simulation input

```bash
cd /mnt/scratch/group/ztang2/pli297/simulation
Rscript generate_simulation_data.R
```

This writes the files under `input/study1/`, `input/study2/`, and `input/study3/`.

Phase 2: run analysis

```bash
cd /mnt/scratch/group/ztang2/pli297/simulation
bash submit_all.sh
```

This is the only submit entrypoint. It submits:

- `run/step0.sbatch`
- `run/step1.sbatch`
- `run/step2_1.sbatch` as a one-chromosome-per-task array for each study
- `run/step2_2.sbatch` for each study
- `run/step2_3.sbatch` for each study
- `run/step3.sbatch` as a feature-level array
- `run/step4.sbatch`

- merge per-chromosome Step2.1 outputs into per-feature `step2_allchr_*` files within each study
- run Step2.3 median correction on the merged per-feature files within each study
- run Step3 meta-analysis across studies as one array task per feature
- run the same four visualization modes used in the main workflow on the meta-analysis results

- phenotype + SNP forest
- phenotype-only Manhattan/QQ
- SNP-only forest across phenotypes
- combined overview across phenotypes

`run/step4.sbatch` reads the meta-analysis results from `/mnt/scratch/group/ztang2/pli297/simulation/output/meta` by default and writes plots to `/mnt/scratch/group/ztang2/pli297/simulation/output` by default:

```bash
sbatch --export=ALL,PHENO=g_Acinetobacter,SNP=chr4:1682869:G:C run/step4.sbatch
```

You can also run the first two scheduler steps manually:

```bash
ENV_FILE=/path/to/your.env
jid0=$(sbatch --parsable --export=ALL,ENV_FILE=${ENV_FILE} run/step0.sbatch)
jid1=$(sbatch --parsable --dependency=afterok:${jid0} --export=ALL,ENV_FILE=${ENV_FILE} run/step1.sbatch)
```

The default simulation root used by these sbatch scripts is:

```bash
/mnt/scratch/group/ztang2/pli297/simulation
```

That directory is expected to contain at least `generate_simulation_data.R`, `input/`, `output/`, `provided/`, and `PALMmGWAS.sif`. You can still override `SIMU_ROOT` or `SIF` from your env file if needed.

Default simulation parameters are written directly inside the sbatch scripts. If needed, you can override them at submit time with `sbatch --export=...`.

Then compute the Step2.1 array size and submit the per-chromosome jobs:

```bash
N_CHROM=22
TOTAL_TASKS=${N_CHROM}

sbatch \
  --dependency=afterok:${jid1} \
  --array=1-${TOTAL_TASKS} \
  --export=ALL,ENV_FILE=${ENV_FILE} \
  run/step2_1.sbatch
```

`run/step1.sbatch` writes `output/<study>/feature_list.txt` automatically for the selected study, so Step2.1 can verify the feature backbone before running the per-chromosome jobs.

`submit_all.sh` is the only submission entrypoint for this simulation workflow. The three studies run independently through `Step0 -> Step2.3` in parallel. Step3 and Step4 run once after all studies complete Step2.3. Step2.1 is always submitted as one chromosome per task, Step2.2 always merges them to `step2_allchr_*`, Step2.3 always overwrites those merged files in place, Step3 always reads `step2_allchr_*.txt`, and Step4 always runs after Step3.

## Main outputs

- `input/study1/`, `input/study2/`, `input/study3/`
- `output/study1/step2_allchr_g_Acinetobacter.txt`
- `output/study1/info_snp.txt`
- `output/study1/info_feature.txt`
- `output/plot/manhattan_g_Acinetobacter.png`

## Slurm array for Step2.1

If you want to parallelize simulated `step2.1` runs, run `Step0` and `Step1` once first, then submit `run/step2_1.sbatch`.

This has already been validated in the simulation workflow: `step2.1` array jobs were launched as one chromosome per task, and the expected per-chromosome outputs were produced.

After the Step2.1 array completes, the simulation workflow merges files sharing the same phenotype across chromosomes. For example:

- `step2_chr1_g_Acinetobacter.txt`
- `step2_chr2_g_Acinetobacter.txt`
- ...

are merged into:

- `step2_allchr_g_Acinetobacter.txt`

After Step2.2 merge, Step2.3 applies median correction to the merged `step2_allchr_*.txt` files and overwrites them in place. The script now also supports direct chromosome-specific correction through `--chrom=1..22`, which targets `step2_chrN_*.txt` for one chromosome without requiring merge first. Step3 still reads the corrected `step2_allchr_*.txt` files and writes `step3_meta_*.txt` through a feature-level array. The visualization step reads those meta files through `run/step4.sbatch`.

Choose the chromosome count for the array size:

```bash
ENV_FILE=/path/to/your.env
N_CHROM=22
TOTAL_TASKS=${N_CHROM}

sbatch --array=1-${TOTAL_TASKS} --export=ALL,ENV_FILE=${ENV_FILE} run/step2_1.sbatch
```

You can override settings at submit time:

```bash
sbatch \
  --array=1-${TOTAL_TASKS} \
  --export=ALL,ENV_FILE=${ENV_FILE},STUDY=study1 \
  run/step2_1.sbatch
```

If you need a custom chromosome list, edit the `CHROMS` default directly in the sbatch script instead of passing a comma-separated list through `sbatch --export`, because Slurm uses commas as variable separators in `--export`.

This script uses 1-based array indexing:

- task `1` is the first chromosome in `CHROMS`
- task `2` is the second chromosome in `CHROMS`
- ...
- the last task is the last chromosome in `CHROMS`

## Note

The input generator now creates three related studies for later meta-analysis work. The current single-study sbatch scripts still default to `study1`, so if you want to run `study2` or `study3` manually, override `STUDY` at submit time or edit the default inside the sbatch files. The default reproducible route is still `bash submit_all.sh`.

The simulation working directory keeps the lightweight metadata under `provided/`. Users regenerate `input/` and `output/` locally by running the scripts above.
