# ===============================================================
# main_binary.R
#
# Master script for the binary simulation study.
# This is the single entry point to reproduce all results
# reported in Section 4.1.1 and Appendix Table 9 of the
# thesis. Running this script in order will:
#
#   1. Load and prepare the four simulated datasets.
#   2. Define model formulas and hyperparameters for each
#      simulation scenario.
#   3. Execute the parallelised 10-fold cross-validation
#      (sourced from R/analysis_binary.R).
#   4. Generate the summary figure (sourced from R/image_binary.R).
#   5. Generate the summary tables (sourced from R/tables_binary.R).
#   6. Write the full configuration to a plain-text record.
#   7. Archive all outputs in a timestamped results folder.
#
# Working directory must be the project root (i.e. the folder
# containing the R/, data/, images/, tables/, configuration/,
# and results/ subdirectories). In RStudio, opening the
# .Rproj file sets this automatically.
#
# Dependencies: mlml, tidyverse, xtable, gridExtra
#   All other packages (lme4, rpart, ranger, MASS, parallel)
#   are loaded inside R/analysis_binary.R.
# ===============================================================


# ---------------------------------------------------------------
# Housekeeping
# ---------------------------------------------------------------
rm(list = ls())

library(mlml)
library(tidyverse)
library(xtable)
library(gridExtra)


# ===============================================================
# SECTION 1: Load and prepare data
#
# The four CSV files are the output of data_simulation.R
# (stored in data/). Each corresponds to one simulation scenario:
#   s1 — linear / additive   (Scenario 1)
#   s2 — nonlinear           (Scenario 2)
#   s3 — interactions        (Scenario 3)
#   s4 — multicollinearity   (Scenario 4)
#
# id must be a factor so that cluster-aware fold construction
# and mixed-effects model fitting treat it as a grouping
# variable rather than a numeric covariate.
# y is kept as integer (0/1) because the binary models expect
# a numeric response, unlike the multinomial case.
# ===============================================================
df_list <- list(
  s1 = read.csv("data/simulated_data_1.csv"),
  s2 = read.csv("data/simulated_data_2.csv"),
  s3 = read.csv("data/simulated_data_3.csv"),
  s4 = read.csv("data/simulated_data_4.csv")
)

df_list <- lapply(df_list, function(d) {
  d$id <- as.factor(d$id)
  d
})


# ===============================================================
# SECTION 2: Cross-validation and scenario configuration
#
# Global CV settings shared across all four scenarios.
# ===============================================================
K_folds      <- 10
seed_folds   <- 42
seed_cluster <- 42

# ---------------------------------------------------------------
# scenario_cfg: one list entry per scenario, each containing
# model formulas and hyperparameter lists for GMERT and GMERF.
#
# The GLMM formula includes a random intercept and a random
# slope on x1, matching the data-generating random-effects
# structure (Table 2). The logit, CART, and RF formulas use
# only fixed predictors.
#
# Scenario 3 (interactions) includes x4 in the GLMM and logit
# formulas because x4 is part of the DGP; CART and RF use
# only x1, x2, x3 as in the original notebook.
#
# GMERT hyperparameters (passed to mlml::fit_gmert_small()):
#   max_iter_out  Maximum outer PQL iterations.
#   tol           Convergence tolerance.
#
# GMERF hyperparameters (passed to mlml::fit_gmerf_small()):
#   max_iter_out  Maximum outer PQL iterations.
#   tol           Convergence tolerance.
# ---------------------------------------------------------------
scenario_cfg <- list(
  s1 = list(
    glmm_formula  = y ~ x1 + x2 + x3 + (1 + x1 | id),
    logit_formula = y ~ x1 + x2 + x3,
    cart_formula  = y ~ x1 + x2 + x3,
    rf_formula    = y ~ x1 + x2 + x3,
    gmert_args = list(max_iter_out = 100, tol = 1e-4),
    gmerf_args = list(max_iter_out = 100, tol = 1e-4)
  ),
  s2 = list(
    glmm_formula  = y ~ x1 + x2 + x3 + (1 + x1 | id),
    logit_formula = y ~ x1 + x2 + x3,
    cart_formula  = y ~ x1 + x2 + x3,
    rf_formula    = y ~ x1 + x2 + x3,
    gmert_args = list(max_iter_out = 100, tol = 1e-2),
    gmerf_args = list(max_iter_out = 100, tol = 1e-2)
  ),
  s3 = list(
    glmm_formula  = y ~ x1 + x2 + x3 + x4 + (1 + x1 | id),
    logit_formula = y ~ x1 + x2 + x3 + x4,
    cart_formula  = y ~ x1 + x2 + x3,
    rf_formula    = y ~ x1 + x2 + x3,
    gmert_args = list(max_iter_out = 50, tol = 1e-4),
    gmerf_args = list(max_iter_out = 50, tol = 1e-4)
  ),
  s4 = list(
    glmm_formula  = y ~ x1 + x2 + x3 + (1 + x1 | id),
    logit_formula = y ~ x1 + x2 + x3,
    cart_formula  = y ~ x1 + x2 + x3,
    rf_formula    = y ~ x1 + x2 + x3,
    gmert_args = list(max_iter_out = 100, tol = 1e-4),
    gmerf_args = list(max_iter_out = 100, tol = 1e-4)
  )
)


# ===============================================================
# SECTION 3: Run analysis, visualization, and tables
#
# Each source() call expects the objects defined above to be
# present in the global environment:
#   R/analysis.R         needs df_list, scenario_cfg, K_folds,
#                        seed_folds, seed_cluster;
#                        produces `results`
#   R/image.R            needs results;
#                        produces images/graph.png
#   R/table.R            needs results;
#                        produces tables/*.png and
#                        prints LaTeX to the console
# ===============================================================
source("R/analysis.R")
source("R/image.R")
source("R/table.R")


# ===============================================================
# SECTION 4: Write configuration record
#
# Saves a plain-text snapshot of all analysis parameters to
# configuration/configuration.txt for inclusion in the
# timestamped results archive and the research archive.
# ===============================================================
config_params <- list(
  K_folds      = K_folds,
  seed_folds   = seed_folds,
  seed_cluster = seed_cluster,
  scenario_cfg = scenario_cfg
)
writeLines(
  capture.output(str(config_params)),
  con = "configuration/configuration.txt"
)


# ===============================================================
# SECTION 5: Archive outputs with timestamp
#
# Creates a uniquely named subfolder under results/ and copies
# all outputs into it so that re-running the script never
# overwrites previous results.
# ===============================================================
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
dir.create(paste0("results/results_", timestamp))

file.copy("tables/acc_table_binary.png",
          paste0("results/results_", timestamp, "/acc_table.png"))
file.copy("tables/f1_mag_table.png",
          paste0("results/results_", timestamp, "/f1_mag_table.png"))
file.copy("tables/f1_min_table.png",
          paste0("results/results_", timestamp, "/f1_min_table.png"))
file.copy("tables/bias_table.png",
          paste0("results/results_", timestamp, "/bias_table.png"))
file.copy("images/graph.png",
          paste0("results/results_", timestamp, "/graph.png"))
file.copy("configuration/configuration.txt",
          paste0("results/results_", timestamp, "/configuration.txt"))

