# MLML: Multilevel Machine Learning
## Research Archive — Master's Thesis

**Title:** MLML: Multilevel Machine Learning — Multilevel machine learning for categorical classification  
**Author:** Paolo Colussi  
**Student number:** 2390191  
**Email:** p.colussi@uu.nl  
**ORCID:** 0009-0009-7048-7559  
**Supervisors:** Joep Burger (Statistics Netherlands, CBS) — Jonas Klingwort (Statistics Netherlands, CBS)  
**Ethics approval:** FERB protocol 25-1979  
**Programme:** Research Master in Methodology and Statistics for Behavioural, Biomedical and Social Sciences — Utrecht University  
**Candidate journal:** Journal of Official Statistics  

---

## Introduction

This repository contains the data and R scripts to reproduce the simulation results for the master's thesis *MLML: Multilevel Machine Learning — Multilevel machine learning for categorical classification*, which investigates whether explicitly modelling the hierarchical structure of clustered data improves predictive performance in classification tasks.

The thesis focuses on GMERT and GMERF, two tree-based extensions that combine flexible fixed-effect estimation with random effects for clustered data. These methods are compared against standard generalized linear models, mixed-effects models, classification trees, and random forests across four controlled simulation scenarios (linear, nonlinear, interactions, multicollinearity) and a real-data application using smartphone-based travel survey data from Statistics Netherlands (CBS).

The simulation study considers both binary and multinomial outcomes. The binary simulation uses a random intercept and random slope on x1, while the multinomial simulation uses a random intercept only with five outcome classes. Both studies use 10-fold cluster-aware cross-validation and evaluate models on accuracy, F1 scores, and relative bias.

The real-data application uses CBS travel survey data and cannot be shared publicly. Further details on data access are provided in the [Permissions and Access](#permissions-and-access) section.

---

## The mlml Package

All GMERT and GMERF models in this archive are fitted using the `mlml` R package, developed as part of this thesis and available at [https://github.com/p3piit/mlml](https://github.com/p3piit/mlml).

The package implements the full family of Generalised Mixed-Effects Regression Tree (GMERT) and Generalised Mixed-Effects Random Forest (GMERF) algorithms for non-Gaussian clustered outcomes, including the multinomial extension derived in Appendix B of the thesis. It provides a unified interface for fitting and predicting these models across binary and multinomial settings.

The core functions used in this archive are:

- `fit_gmert_small()` — fits a binary GMERT model via the PQL-EM algorithm (Algorithm 1 of the thesis).
- `fit_gmerf_small()` — fits a binary GMERF model, replacing the tree step with a random forest.
- `fit_gmert_cat()` — fits a multinomial GMERT model using the categorical extension (Appendix B).
- `fit_gmerf_cat()` — fits a multinomial GMERF model using the categorical extension (Appendix B).
- `predict_gmert()` / `predict_gmerf()` — generate predictions for new observations; population-level predictions (random effects set to zero) are used throughout the cross-validation.
- `predict_gmert_cat()` / `predict_gmerf_cat()` — the multinomial equivalents, with an additional `prob_saved = TRUE` argument to return the full N × K probability matrix.

The package can be installed directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("p3piit/mlml")
```

The package is also managed via `renv` in this archive, so `renv::restore()` will install the correct version automatically without requiring a manual GitHub installation.

---

## Data

This repository contains data in the following structure:

``` plaintext
This repository contains data in the following structure:
├── Binary_analysis
│   ├── R
│   │   ├── analysis.R
│   │   ├── data_simulation.R
│   │   ├── image.R
│   │   ├── main_binary.R
│   │   └── table.R
│   ├── configuration
│   │   └── configuration.txt
│   ├── data
│   │   ├── simulated_data1.csv
│   │   ├── simulated_data2.csv
│   │   ├── simulated_data3.csv
│   │   └── simulated_data4.csv
│   ├── images
│   │   └── graph.png
│   ├── results
│   │   └── results_20260212-182435
│   │       ├── acc_table.png
│   │       ├── bias_table.png
│   │       ├── configuration.txt
│   │       ├── f1_mag_table.png
│   │       ├── f1_min_table.png
│   │       └── graph.png
│   └── tables
│       ├── acc_table.png
│       ├── bias_table.png
│       ├── f1_mag_table.png
│       └── f1_min_table.png
├── LICENSE.md
├── Multinominal_analysis
│   ├── R
│   │   ├── analysis.R
│   │   ├── data_simulation.R
│   │   ├── image.R
│   │   ├── main_multinominal.R
│   │   └── table.R
│   ├── configuration
│   │   └── configuration.txt
│   ├── data
│   │   ├── simulated_data1.csv
│   │   ├── simulated_data2.csv
│   │   ├── simulated_data3.csv
│   │   └── simulated_data4.csv
│   ├── images
│   │   └── graph.png
│   ├── results
│   │   └── results_20260412-113837
│   │       ├── acc_table.png
│   │       ├── configuration.txt
│   │       ├── f1_table.png
│   │       ├── graph.png
│   │       └── wf1_table.png
│   └── tables
│       ├── acc_table.png
│       ├── f1_table.png
│       └── wf1_table.png
├── README.md
├── mlml_thesis_research_archive.Rproj
└── renv.lock
```

### `Binary_analysis/`

Contains all scripts, data, and outputs for the binary simulation study (Section 3.2 and Section 4.1.1 of the thesis).

#### `Binary_analysis/data/`

Four CSV files containing the simulated clustered binary datasets, one per scenario. Each dataset contains 30 clusters of 60 observations each (N = 1800), with columns `id` (cluster identifier), `y` (binary outcome), and predictors `x1`–`x4` depending on the scenario.

- `simulated_data1.csv` — Scenario 1: linear / additive fixed component, mildly correlated predictors (rho = 0.2).
- `simulated_data2.csv` — Scenario 2: nonlinear fixed component; x2 enters via log(x2).
- `simulated_data3.csv` — Scenario 3: interaction fixed component; x2 × x4 product term.
- `simulated_data4.csv` — Scenario 4: linear fixed component with strongly correlated predictors (rho = 0.9).

All datasets are generated using `data_simulation.R` with `seed = 123`. The random-effects structure is identical across scenarios: a random intercept and a random slope on x1 with covariance matrix D = [[0.64, 0.08], [0.08, 0.25]].

#### `Binary_analysis/R/`

- `data_simulation.R` — Generates the four binary simulation datasets and writes them to `data/`.
- `analysis.R` — Defines helper functions, the single-fold CV worker, the summary function, and runs the parallelised 10-fold cross-validation across all four scenarios. Produces the `results` object.
- `image.R` — Takes `results` and produces the faceted dot-and-errorbar figure (Figure 1 of the thesis). Writes to `images/graph.png`.
- `table.R` — Takes `results` and produces formatted summary tables for accuracy, majority-class F1, minority-class F1, and relative bias. Writes PNG and LaTeX output to `tables/`.
- `main_binary.R` — Master script. Sources all of the above in order, writes the configuration record, and archives outputs in a timestamped subfolder under `results/`.

#### `Binary_analysis/configuration/`

- `configuration.txt` — Plain-text snapshot of all analysis parameters (number of folds, seeds, model formulas, hyperparameters) written automatically by `main_binary.R` at each run.

#### `Binary_analysis/results/`

Each run of `main_binary.R` creates a timestamped subfolder (e.g. `results_20260412-134604/`) containing a frozen copy of all outputs and the configuration file. This ensures previous results are never overwritten.

#### `Binary_analysis/tables/` and `Binary_analysis/images/`

Contain the most recent outputs of `table.R` and `image.R` respectively. These are overwritten on each run; archived copies are preserved in `results/`.

---

### `Multinominal_analysis/`

Contains all scripts, data, and outputs for the multinomial simulation study (Section 3.2 and Section 4.1.2 of the thesis).

#### `Multinominal_analysis/data/`

Four CSV files containing the simulated clustered multinomial datasets, one per scenario. Each dataset contains 100 clusters of 18 observations each (N = 1800), with columns `id` (cluster identifier), `y` (five-class outcome, class 5 as reference), and predictors `x1`–`x10`.

- `simulated_data1.csv` — Scenario 1: linear / additive fixed component (rho = 0.2).
- `simulated_data2.csv` — Scenario 2: nonlinear fixed component (log, square, sin, cos, tanh transformations).
- `simulated_data3.csv` — Scenario 3: interaction fixed component (two-way products of predictors).
- `simulated_data4.csv` — Scenario 4: linear fixed component with strongly correlated predictors (rho = 0.85).

All datasets are generated using `data_simulation.R` with `seed = 123`. The random-effects structure is identical across scenarios: a class-specific random intercept with SD = 5 for each of the four non-reference log-odds equations.

#### `Multinominal_analysis/R/`

- `data_simulation.R` — Generates the four multinomial simulation datasets and writes them to `data/`.
- `analysis.R` — Defines helper functions, the single-fold CV worker, the summary function, and runs the parallelised 10-fold cross-validation across all four scenarios. Produces the `results` object.
- `image.R` — Takes `results` and produces the faceted dot-and-errorbar figure (Figure 2 of the thesis). Writes to `images/graph.png`.
- `table.R` — Takes `results` and produces formatted summary tables for accuracy, macro F1, and weighted F1. Writes PNG and LaTeX output to `tables/`.
- `main_multinominal.R` — Master script. Sources all of the above in order, writes the configuration record, and archives outputs in a timestamped subfolder under `results/`.

#### `Multinominal_analysis/configuration/`

- `configuration.txt` — Plain-text snapshot of all analysis parameters written automatically by `main_multinominal.R` at each run.

#### `Multinominal_analysis/results/`

Each run of `main_multinominal.R` creates a timestamped subfolder containing a frozen copy of all outputs and the configuration file.

#### `Multinominal_analysis/tables/` and `Multinominal_analysis/images/`

Contain the most recent outputs of `table.R` and `image.R` respectively.

---

### `mlml_thesis_research_archive.Rproj`

The RStudio project file. Opening this file in RStudio sets the working directory to the project root automatically. All scripts use `here::here()` to construct paths relative to this root, so they run correctly regardless of the platform or the directory from which they are called.

### `renv.lock`

Records the exact versions of all R packages used in the analysis. Used by `renv::restore()` to recreate the R environment (see [Prerequisites](#prerequisites) below).

---

## Reproducing Results

### Prerequisites

To reproduce the simulation results, the following are needed:

- **R** (≥ 4.6.0)
- **RStudio** (recommended; ensures the working directory is set correctly by the `.Rproj` file)
- The **`renv`** package for restoring the R environment
- The **`mlml`** package (see [The mlml Package](#the-mlml-package) above)

Install the required packages by running the following in the R console after opening the `.Rproj` file:

```r
renv::restore()
```

> **Hardware recommendation:** The cross-validation runs in parallel using a PSOCK cluster with one worker per fold (10 folds). Ideally, the machine running the analysis should have **at least 11 cores** — 10 for the parallel workers and 1 for the main R session. The analysis can in principle run on fewer cores, but each missing core forces workers to share CPU time, which increases computation time substantially. Thread-level parallelism is disabled within each worker to prevent over-subscription across processes.

### Running the Scripts

All results can be reproduced by running the two master scripts. Open `mlml_thesis_research_archive.Rproj` in RStudio first, then run from the R console:

```r
# Binary simulation study (Section 4.1.1)
source(here::here("Binary_analysis", "R", "main_binary.R"))

# Multinomial simulation study (Section 4.1.2)
source(here::here("Multinominal_analysis", "R", "main_multinominal.R"))
```

Each master script executes the following steps automatically:

1. Reads the simulated datasets from `data/`.
2. Runs the parallelised 10-fold cross-validation across all four scenarios.
3. Generates the summary figure and saves it to `images/`.
4. Generates the summary tables and saves them to `tables/`.
5. Writes the configuration record to `configuration/`.
6. Archives all outputs in a new timestamped subfolder under `results/`.

If you wish to regenerate the simulated datasets from scratch before running the analysis, source `data_simulation.R` first:

```r
source(here::here("Binary_analysis",       "R", "data_simulation.R"))
source(here::here("Multinominal_analysis", "R", "data_simulation.R"))
```

> **Note on computation time:** The binary simulation takes approximately **2–3 hours** on a machine with 11 or more cores. The multinomial simulation is substantially more demanding due to the iterative PQL-EM estimation of K − 1 = 4 separate tree or forest components at each iteration, and takes approximately **16–18 hours** under the same conditions. Both analyses were run on a machine with 11 cores and 32 GB of RAM. Computation time will increase on machines with fewer cores.

> **Note on reproducibility:** All random seeds are fixed in the master scripts (`seed_folds`, `seed_cluster`) and inside the data-generating functions (`seed = 123`). Results should be exactly reproducible on the same platform and R version. Minor numerical differences may arise across operating systems due to platform-specific floating-point behaviour in `lme4` and `ranger`.

---

## Real-Data Application

The real-data application (Section 3.3 and Section 4.2 of the thesis) uses smartphone-based travel survey data collected by Statistics Netherlands (CBS) as part of the CBS mobility research programme. The dataset contains 5,020 tracks nested within 369 respondents, with approximately 300 predictors derived from GPS and accelerometer signals recorded by the CBS app. For each track, the app records the transport mode used (car, bike, walk, bus, metro, tram, or train), which serves as the outcome variable in the analysis.

Three prediction tasks are considered: a binary task distinguishing car from all other modes, a binary task distinguishing urban public transport (bus, metro, tram) from all other modes, and a full multinomial task predicting all seven transport-mode categories jointly.

This data cannot be shared publicly due to CBS data-access restrictions and privacy considerations applying to survey respondents. Researchers wishing to reproduce the real-data results or obtain access to the dataset should contact Statistics Netherlands directly:

- **Joep Burger** — j.burger@cbs.nl
- **Jonas Klingwort** — j.klingwort@cbs.nl

The analysis code for the real-data application is structurally identical to the simulation scripts provided in this archive. The only additions are a data-cleaning and preprocessing step that transforms the raw GPS and accelerometer signals into the predictor matrix used for modelling, and minor adjustments to the cross-validation design to accommodate the strongly unbalanced cluster sizes and class distribution present in the real data. This code can also be shared upon request to the supervisors above, subject to CBS approval.

---

## Ethics and Privacy

Ethics approval was granted by the Ethics Review Board of the Faculty of Social and Behavioural Sciences at Utrecht University (FERB protocol 25-1979).

Due to the synthetic nature of the data used in both simulation studies, no privacy concerns are relevant. The simulated datasets are generated from known parametric models and do not contain any real-world personal data.

The real-data application uses travel survey data collected by CBS. Respondents provided informed consent for the use of their mobility data for research purposes. The data are stored and processed in accordance with the CBS data-governance framework and the privacy regulations of the Faculty of Social and Behavioural Sciences of Utrecht University. No personally identifiable information is included in the analysis.

---

## Permissions and Access

This archive is publicly available at [https://github.com/p3piit/mlml_thesis_research_archive](https://github.com/p3piit/mlml_thesis_research_archive).

For any questions regarding this archive, do not hesitate to contact:

**Paolo Colussi**  
Email: p.colussi@uu.nl

---

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE.md) file for details.
