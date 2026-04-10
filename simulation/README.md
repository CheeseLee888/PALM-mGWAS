# Simulation Dataset

This directory provides a reproducible single-study simulation intended to show a larger sample size than `example/` and give users a more realistic sense of `step2` runtime and memory usage.

The directory is split into two parts:

- `simulation/provided/` stores the small source files that should be kept in the repository
- `simulation/input/` and `simulation/output/` are generated locally and should not be pushed to GitHub

Core design:

- The sample size is fixed at `502`, defined in `simulation/provided/sample_config.tsv`
- The SNP list and all `CHR/POS/A1/A2/AF` values are taken from `simulation/provided/snp_reference.tsv`
- Feature names and covariate names are stored in `simulation/provided/feature_names.txt` and `simulation/provided/covariate_names.txt`
- The generated `abd` and `cov` tables follow the same overall column structure as `example/input/study1`
- The random seed is fixed at `20260409`, so any user running the same script will generate identical inputs and outputs
- A strong genetic signal is planted in `g_Acinetobacter` so the Manhattan plot contains clear high peaks

## One-command run

From the repository root, run:

```bash
bash simulation/run_simulation.sh
```

The script will automatically:

1. Generate `simulation/input/study1/abd.txt`
2. Generate `simulation/input/study1/cov.txt`
3. Generate `simulation/input/study1/geno.{bed,bim,fam}`
4. Run `Step0 -> Step1 -> Step2`
5. Write the Manhattan/QQ plot for `g_Acinetobacter`

## Main outputs

- `simulation/output/study1/step2_allchr_g_Acinetobacter.txt`
- `simulation/output/study1/info_snp.txt`
- `simulation/output/study1/info_feature.txt`
- `simulation/output/plot/manhattan_g_Acinetobacter.png`

## Note

This simulation currently uses a single-study workflow, so the visualization step reads `step2` outputs directly instead of `step3` meta-analysis outputs.

The repository is intended to keep only the lightweight metadata under `simulation/provided/`. Users regenerate `simulation/input/` and `simulation/output/` locally by running the script above.
