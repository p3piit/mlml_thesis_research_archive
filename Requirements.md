# Requirements

## R Version

R ≥ 4.6.0

The analysis was developed and run on R 4.6.0.
Minor numerical differences may arise on other operating systems due to
platform-specific floating-point behaviour in `lme4` and `ranger`.

## Key Packages

The full list of package versions is recorded in `renv.lock`. The most
important packages are listed below for quick reference.

| Package    | Version  | Role                                              |
|------------|----------|---------------------------------------------------|
| mlml       | 0.3.2    | GMERT and GMERF model fitting and prediction      |
| lme4       | 2.0-1    | GLMM fitting (binary simulation)                  |
| nnet       | 7.3-19   | Multinomial logit fitting                         |
| mclogit    | 0.9.15   | Mixed-effects multinomial logit fitting           |
| rpart      | 4.1.27   | Classification tree fitting                       |
| ranger     | 0.18.0   | Random forest fitting                             |
| tidyverse  | 2.0.0    | Data manipulation and visualization               |
| xtable     | 1.8-8    | LaTeX table output                                |
| gridExtra  | 2.3      | PNG table output                                  |
| here       | 1.0.2    | Project-relative path construction                |
| MASS       | 7.3-60   | Multivariate normal simulation                    |

## Hardware

The analysis was run on a machine with the following specifications:

- **CPU**: 20 cores
- **RAM**: 32 GB

At least **11 cores** are recommended — 10 for the parallel CV workers and 1
for the main R session. The analysis can run on fewer cores but computation
time will increase substantially (see computation time notes in the README).

## Restoring the Environment

All package versions are managed via `renv`. To restore the exact environment
used for the analysis, open `mlml_thesis_research_archive.Rproj` in RStudio
and run:

```r
renv::restore()
```

This will install all packages at the versions recorded in `renv.lock`,
including `mlml` from GitHub. An internet connection is required.
