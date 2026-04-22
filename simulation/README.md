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

- `provided/`: reference inputs used by the simulation generator, including feature names, covariate names, SNP definitions, and study sample-size settings
- `run/`: Slurm submission scripts for each simulation stage
- `generate_simulation_data.R`: the R script that creates the simulation inputs for `study1`, `study2`, and `study3`

Core design:

- The generator now writes three studies: `study1`, `study2`, and `study3`
- The baseline sample size is defined in `simulation/provided/sample_config.tsv`, and `study2` and `study3` use different sample sizes to mimic heterogeneous cohorts
- The SNP list and all `CHR/POS/A1/A2/AF` values are taken from `simulation/provided/snp_reference.tsv`
- Feature names and covariate names are stored in `simulation/provided/feature_names.txt` and `simulation/provided/covariate_names.txt`
- The three studies intentionally share the same SNP and feature backbone, but `study2` and `study3` each drop a small subset of SNPs and features so later meta-analysis can test imperfect overlap
- The generated `abd` and `cov` tables follow the same overall column structure as `example/input/study1`
- The random seed is fixed at `20260409`, so any user running the same script will generate identical inputs and outputs
- A strong genetic signal is planted in `g_Acinetobacter` so the Manhattan plot contains clear high peaks


## Reproducible Route

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
If you want to reproduce the same results currently stored under `output/`, skip `step2.2`.

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

## Outputs

For routine inspection, users only need to look at `output/plot/`, which contains the downstream reporting outputs from the simulation workflow.
