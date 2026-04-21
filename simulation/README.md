# Simulation Dataset

This directory provides a reproducible multi-study simulation intended to show a larger sample size than `example/` and give users a more realistic sense of runtime and memory usage.

The working directory is:

- `/path/to/simulation`

That directory should contain:

- `provided/`
- `run/`
- `generate_simulation_data.R`

Users generate `input/` and `output/` locally by running the simulation workflow. `PALMmGWAS.sif` is not shipped in this directory and should be built by following the container instructions in the guide.

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

The simulation workflow has been run successfully through `Step0 -> Step2.2` on the cluster layout described below.

In particular, `step2.1` has been verified under Slurm with a `one chromosome per job` decomposition where each job processes all features for that chromosome:

- chromosome subsetting is controlled by `--chrom`
- all features on that chromosome are modeled together with `--featureColList=NULL`
- parallelism is handled at the scheduler level by submitting multiple jobs, rather than inside one R process

This means the simulation workflow now follows one fixed route by default: `step2.1` runs as `chrom only` with all features, `step2.2` runs as `chrom only` and overwrites those chromosome-split files in place, `merge` then combines the corrected chromosome-split files into final `step2_allchr_*`, and `step3` reads those merged `step2_allchr_*.txt` files.

## Reproducible Route

The simulation workflow has one supported route.

Phase 1: generate simulation input

```bash
cd /path/to/simulation
Rscript generate_simulation_data.R
```

This writes the files under `input/study1/`, `input/study2/`, and `input/study3/`.

Phase 2: submit jobs step by step

```bash
cd /path/to/simulation
source .env
mkdir -p logs
```

Before submitting the next stage, check that the current stage has really finished cleanly:

```bash
squeue -u "$USER"
sacct -j <jobid> --format=JobID,JobName%20,State,ExitCode,Elapsed
```

Use `squeue` to see whether there are still jobs running or pending. Use `sacct -j <jobid>` to inspect whether a job finished as `COMPLETED` or ended as `FAILED` or `CANCELLED`. Make sure the previous step completed correctly before moving on to the next one.

Submit `Step0 -> Step1` for each study. Wait for `step0` to finish before submitting `step1`.

```bash
STUDY=study1 sbatch run/step0.sbatch
STUDY=study2 sbatch run/step0.sbatch
STUDY=study3 sbatch run/step0.sbatch

STUDY=study1 sbatch run/step1.sbatch
STUDY=study2 sbatch run/step1.sbatch
STUDY=study3 sbatch run/step1.sbatch
```

For each study, submit `Step2.1 -> Step2.2 -> merge` only after the previous stage has finished.
If you want to reproduce the same results currently stored under `output/`, do not run `step2.2`.

```bash
STUDY=study1 sbatch --array=1-22 run/step2_1.sbatch
STUDY=study2 sbatch --array=1-22 run/step2_1.sbatch
STUDY=study3 sbatch --array=1-22 run/step2_1.sbatch

STUDY=study1 sbatch --array=1-22 run/step2_2.sbatch
STUDY=study2 sbatch --array=1-22 run/step2_2.sbatch
STUDY=study3 sbatch --array=1-22 run/step2_2.sbatch

n_feature_study1=$(wc -l < output/study1/feature_list.txt)
n_feature_study2=$(wc -l < output/study2/feature_list.txt)
n_feature_study3=$(wc -l < output/study3/feature_list.txt)

STUDY=study1 sbatch --array=1-${n_feature_study1} run/merge.sbatch
STUDY=study2 sbatch --array=1-${n_feature_study2} run/merge.sbatch
STUDY=study3 sbatch --array=1-${n_feature_study3} run/merge.sbatch
```

This means:

- `run/step2_1.sbatch` runs as a one-chromosome-per-task array
- `run/step2_2.sbatch` runs as a one-chromosome-per-task array
- `run/merge.sbatch` runs as a one-feature-per-task array

The default simulation root used by these sbatch scripts is:

```bash
/path/to/simulation
```

That directory is expected to contain at least `generate_simulation_data.R`, `input/`, `output/`, `provided/`, and `PALMmGWAS.sif`.

After all three study-level merge stages have completed, write `study_dirs.tsv`, then submit `step3`. Wait for `step3` to finish before submitting `step4`.

```bash
mkdir -p output/meta

cat > output/meta/study_dirs.tsv <<EOF
study1	${PWD}/output/study1
study2	${PWD}/output/study2
study3	${PWD}/output/study3
EOF

n_feature_meta=$(wc -l < output/study1/feature_list.txt)
sbatch --array=1-${n_feature_meta} run/step3.sbatch

sbatch run/step4.sbatch
```

## Main outputs

- `input/study1/`, `input/study2/`, `input/study3/`
- `output/study1/step2_allchr_g_Acinetobacter.txt`
- `output/study1/info_snp.txt`
- `output/study1/info_feature.txt`
- `output/plot/manhattan_g_Acinetobacter.png`

## Step2 Structure

`step2.1` is always launched as a fixed chromosome-level array, and this has already been validated in the simulation workflow.

After the Step2.1 array completes, the simulation workflow runs Step2.2 median correction chromosome by chromosome, overwriting files such as:

- `step2_chr1_g_Acinetobacter.txt`
- `step2_chr2_g_Acinetobacter.txt`
- ...

After the Step2.2 array completes, the simulation workflow merges corrected files sharing the same phenotype across chromosomes. These chromosome shards:

- `step2_chr1_g_Acinetobacter.txt`
- `step2_chr2_g_Acinetobacter.txt`
- `...`

are merged into:

- `step2_allchr_g_Acinetobacter.txt`

Step2.2 applies median correction to the selected Step2 scope. In the simulation workflow it runs on chromosome-split files with `--overwriteOutput=TRUE`, so merge reads `step2_chr1` through `step2_chr22` after correction and writes final `step2_allchr_*.txt`. If the chromosome-split files do not exist, merge will not work. The visualization step reads the downstream meta files through `run/step4.sbatch`.

Step2.1 and Step2.2 always use a fixed `1-22` Slurm array in this simulation workflow.

This script uses 1-based array indexing:

- task `1` is chromosome 1
- task `2` is chromosome 2
- ...
- task `22` is chromosome 22

## Note

The input generator creates three related studies for later meta-analysis work. The intended reproducible route is the step-by-step submission sequence above.

The simulation working directory keeps the lightweight metadata under `provided/`. Users regenerate `input/` and `output/` locally by running the scripts above.
