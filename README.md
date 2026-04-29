# PALM-mGWAS

`PALM-mGWAS` is an R package and command-line workflow for microbiome
genome-wide association studies based on the `PALM` framework. It supports
input checking, PALM null model fitting, per-study association summary
generation, multi-study meta-analysis, and downstream reporting.

Current genotype input support includes PLINK (`.bed/.bim/.fam`) and VCF
(`.vcf`, `.vcf.gz`, `.vcf.bgz`).

For installation, command-line usage, input formats, examples, Docker, and HPC
instructions, see the full guide:

<https://cheeselee888.github.io/PALM-mGWAS/>

## Package Information

- R package name: `PALMmGWAS`
- Version: `1.0`
- License: `GPL (>= 2)`
- Required R version: `R >= 4.1.0`

## Repository Map

This repository contains both the R package source and reproducible workflow
materials.

| Path | Purpose |
| --- | --- |
| `DESCRIPTION` | R package metadata: package name, version, authors, license, R version, dependencies, project URL, and issue tracker. |
| `NAMESPACE` | Exported functions and imported symbols used by R when loading the package. This file is generated from roxygen comments. |
| `LICENSE` | GPL license text for the GitHub repository. |
| `README.md` | This repository overview. Detailed usage is intentionally kept in the online guide. |
| `.Rbuildignore` | Files and directories excluded from the R source package tarball, such as examples, docs, local environments, and large runtime artifacts. |
| `.gitignore` | Local files ignored by Git, including `.env`, `.pixi/`, R check outputs, Docker/SIF archives, and temporary files. |
| `.dockerignore` | Files excluded from Docker build context. |
| `pixi.toml` / `pixi.lock` | Reproducible pixi environment definition and lockfile used for local development and source runs. |

## Source Code

| Path | Purpose |
| --- | --- |
| `R/fitNULL.R` | Fits the PALM null model from abundance and covariate inputs. This is the package function behind Step1. |
| `R/getSummary.R` | Reads PLINK or VCF genotype input and generates per-feature association summary statistics. This is the package function behind Step2.1. |
| `R/correctSummary.R` | Applies median- or tune-based compositional correction to Step2 result files. This is the package function behind Step2.2. |
| `R/metaSummary.R` | Runs multi-study meta-analysis across per-study summary files. This is the package function behind Step3. |
| `R/reporting.R` | Generates reporting outputs, including Manhattan plots, QQ plots, combined-hit summaries, and forest plots. This is the package code behind Step4. |
| `R/generateInfo.R` | Writes auxiliary information tables, such as feature, SNP, sequencing-depth, and run summaries. |
| `R/genoInput.R` | Handles genotype input details, including PLINK and VCF sample/SNP parsing. |
| `R/utils.R` | Shared helper functions used across the package. |

## Command-Line Scripts

The `extdata/` directory contains user-facing command-line wrappers around the
package functions. These scripts are the recommended workflow entry points.

| Path | Purpose |
| --- | --- |
| `extdata/step0_checkInput.R` | Checks sample IDs, filters/aligns abundance and covariate tables, and verifies genotype sample compatibility. |
| `extdata/step1_null.R` | Fits the null model and writes the Step1 `.rda` object. |
| `extdata/step2_1_summary.R` | Generates Step2.1 association summaries for one chromosome or all chromosomes. |
| `extdata/step2_2_correction.R` | Applies compositional correction to Step2 outputs. |
| `extdata/step3_meta.R` | Runs feature-wise meta-analysis across studies. |
| `extdata/step4_reporting.R` | Generates downstream plots and reporting tables. |

## Documentation

| Path | Purpose |
| --- | --- |
| `man/` | R help files (`.Rd`) generated from roxygen comments in `R/`. These power `?function_name` after package installation. |
| `docs/index.html` | GitHub Pages redirect to the main guide. |
| `docs/guide.html` | Full user guide with installation routes, workflow commands, parameter tables, Docker/SIF instructions, and examples. |
| `docs/assets/plot/` | Plot images used by the online guide. |

## Example Data

| Path | Purpose |
| --- | --- |
| `example/input/` | Small example abundance, covariate, genotype, and study-list inputs. |
| `example/output/` | Example outputs from the workflow, including Step2 summaries, meta-analysis results, and plots. |

These files are included so users can inspect expected input and output shapes
without running a full simulation.

## Simulation Workflow

| Path | Purpose |
| --- | --- |
| `simulation/README.md` | Detailed instructions for the larger reproducible simulation. |
| `simulation/generate_simulation_data.R` | Generates multi-study simulation input data. |
| `simulation/provided/` | Reference files used by the simulation generator, including feature names, covariate names, sample settings, and SNP definitions. |
| `simulation/run/` | Slurm scripts for Step0, Step1, Step2.1, Step2.2, Step3, and Step4 on an HPC cluster. |

Generated simulation `input/`, `output/`, logs, and local `.env` files are not
intended to be committed.

## Container and Third-Party Materials

| Path | Purpose |
| --- | --- |
| `docker/Dockerfile` | Docker recipe for building a PALM-mGWAS runtime with pixi, PALM, PALMmGWAS, and R dependencies installed. |
| `thirdParty/PALM_0.1.0.tar.gz` | Bundled PALM source package used by the local and container installation workflows. |
| `thirdParty/README_PALM.md` | Notes about the bundled PALM materials. |
| `thirdParty/*.pdf` | PALM package and manuscript reference materials. |

Large local runtime artifacts such as `PALMmGWAS.sif` and Docker archive tarballs
are ignored and should not be committed.

## Development Notes

- Use roxygen comments in `R/` as the source of truth for `NAMESPACE` and
  `man/*.Rd`.
- Run `roxygen2::roxygenise()` or `devtools::document()` after changing exported
  functions or documentation.
- Build and check the package from a source tarball before release.
- Keep local runtime state out of Git: `.env`, `.pixi/`, `.Rcheck/`, SIF images,
  Docker archives, and generated simulation outputs are ignored.
