# Simulation Dataset

This directory provides a reproducible multi-study simulation intended to show a larger sample size than `example/`, give users a more realistic sense of `step2.1` runtime and memory usage, and provide compatible inputs for later meta-analysis testing.

The working directory is:

- `/mnt/scratch/group/ztang2/pli297/simulation`

That directory directly contains:

- `provided/`
- `input/`
- `output/`
- `PALMmGWAS.sif`
- the submission scripts in this folder

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

In particular, `step2.1` has been verified under Slurm with a `chrom x feature block` decomposition:

- chromosome subsetting is controlled by `--chrom`
- feature subsetting is controlled by `--featureColList`
- parallelism is handled at the scheduler level by submitting multiple jobs, rather than inside one R process

This means the current implementation supports user-managed `chrom x feature` parallel execution for `step2.1`, followed by a separate `step2.2` merge stage and `step2.3` correction stage.

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
bash submit_pipeline.sh
```

This submits `Step0 -> Step1 -> Step2.1`.

If you want one command that submits the full workflow from `Step0` through visualization, run:

```bash
cd /mnt/scratch/group/ztang2/pli297/simulation
bash submit_all.sh
```

This submits:

- `step0_align_inputs.sbatch`
- `step1_fit_null_model.sbatch`
- `submit_step2_array.sbatch`
- `run_step2_2.sbatch`
- `run_step2_3_correction.sbatch`
- `run_step3_meta.sh` as a feature-level array
- `visualize_meta_results.sh`

After all Step2.1 array jobs finish, run:

```bash
cd /mnt/scratch/group/ztang2/pli297/simulation
sbatch run_step2_2.sbatch
sbatch run_step2_3_correction.sbatch
sbatch submit_step3_meta_array.sbatch
bash visualize_meta_results.sh
```

This now performs four stages across four commands:

- merge per-chromosome Step2.1 outputs into per-phenotype `step2_allchr_*` files within each study
- run Step2.3 median correction on the merged per-feature files within each study
- run Step3 meta-analysis across studies as one array task per feature
- run the same four visualization modes used in the main workflow on the meta-analysis results

- phenotype + SNP forest
- phenotype-only Manhattan/QQ
- SNP-only forest across phenotypes
- combined overview across phenotypes

By default, `run_step2_2.sbatch` uses `study1,study2,study3`. If needed, you can restrict the Step2.2 stage to one study or a subset:

```bash
STUDY=study2 bash run_step2_2.sbatch
STUDIES=study1,study3 bash run_step2_2.sbatch
```

Submit the split merge and meta steps as Slurm jobs so runtime and memory are recorded separately:

```bash
sbatch run_step2_2.sbatch
sbatch --export=ALL,STUDY=study2 run_step2_2.sbatch
sbatch submit_step3_meta_array.sbatch
```

`visualize_meta_results.sh` reads the meta-analysis results from `/mnt/scratch/group/ztang2/pli297/simulation/output/meta` by default and writes plots to `/mnt/scratch/group/ztang2/pli297/simulation/output` by default:

```bash
PHENO=g_Acinetobacter SNP=chr4:1682869:G:C bash visualize_meta_results.sh
```

You can also run the analysis chain manually in three scheduler steps:

```bash
ENV_FILE=/path/to/your.env
jid0=$(sbatch --parsable --export=ALL,ENV_FILE=${ENV_FILE} step0_align_inputs.sbatch)
jid1=$(sbatch --parsable --dependency=afterok:${jid0} --export=ALL,ENV_FILE=${ENV_FILE} step1_fit_null_model.sbatch)
```

The default simulation root used by these sbatch scripts is:

```bash
/mnt/scratch/group/ztang2/pli297/simulation
```

That directory is expected to contain at least `generate_simulation_data.R`, `input/`, `output/`, `provided/`, and `PALMmGWAS.sif`. You can still override `SIMU_ROOT` or `SIF` from your env file if needed.

Default simulation parameters are written directly inside the sbatch scripts. If needed, you can override them at submit time with `sbatch --export=...`.

Then compute the Step2.1 array size and submit the chromosome x feature-block jobs:

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
  step2_run_array.sbatch
```

`step1_fit_null_model.sbatch` writes `output/<study>/feature_list.txt` automatically for the selected study, so Step2.1 can use it directly.

`submit_pipeline.sh` submits the analysis stages through Step2.1. The Step2.2 merge, Step2.3 correction, meta stage, and visualization stage are intentionally run later as separate manual commands after Step2.1 finishes.

## Main outputs

- `input/study1/`, `input/study2/`, `input/study3/`
- `output/study1/step2_allchr_g_Acinetobacter.txt`
- `output/study1/info_snp.txt`
- `output/study1/info_feature.txt`
- `output/plot/manhattan_g_Acinetobacter.png`

## Slurm array for Step2.1

If you want to parallelize simulated `step2.1` runs across both chromosome and feature subsets, run `Step0` and `Step1` once first, then submit `step2_run_array.sbatch`.

This has already been validated in the simulation workflow: `step2.1` array jobs were launched across multiple chromosomes and feature blocks, and the expected per-chromosome outputs were produced.

After the Step2.1 array completes, the simulation workflow merges files sharing the same phenotype across chromosomes. For example:

- `step2_chr1_g_Acinetobacter.txt`
- `step2_chr2_g_Acinetobacter.txt`
- ...

are merged into:

- `step2_allchr_g_Acinetobacter.txt`

After Step2.2 merge, Step2.3 updates the merged `step2_allchr_*.txt` files in place using median correction. The meta-analysis step then reads those corrected files and writes `step3_meta_*.txt` through a feature-level array. The visualization step reads those meta files.

Choose a feature block size, for example `10` features per task:

```bash
ENV_FILE=/path/to/your.env
FEATURE_BLOCK=10
N_FEATURE=$(wc -l < output/study1/feature_list.txt)
N_BLOCK=$(( (N_FEATURE + FEATURE_BLOCK - 1) / FEATURE_BLOCK ))
N_CHROM=22
TOTAL_TASKS=$(( N_BLOCK * N_CHROM ))

sbatch --array=1-${TOTAL_TASKS} --export=ALL,ENV_FILE=${ENV_FILE} step2_run_array.sbatch
```

You can override settings at submit time:

```bash
sbatch \
  --array=1-${TOTAL_TASKS} \
  --export=ALL,ENV_FILE=${ENV_FILE},FEATURE_BLOCK=10,STUDY=study1 \
  step2_run_array.sbatch
```

If you need a custom chromosome list, edit the `CHROMS` default directly in the sbatch script instead of passing a comma-separated list through `sbatch --export`, because Slurm uses commas as variable separators in `--export`.

This script uses 1-based array indexing:

- task `1` is chromosome `1`, feature block `1`
- task `2` is chromosome `1`, feature block `2`
- ...
- after the last block on chromosome `1`, the next task moves to chromosome `2`

## Note

The input generator now creates three related studies for later meta-analysis work. The current Slurm analysis scripts still default to `study1`, so if you want to run `study2` or `study3`, override `STUDY` at submit time or edit the default inside the sbatch files. By contrast, `run_step2_2.sbatch` defaults to combining all three studies for later meta-analysis unless you restrict `STUDY` or `STUDIES`.

The simulation working directory keeps the lightweight metadata under `provided/`. Users regenerate `input/` and `output/` locally by running the scripts above.
