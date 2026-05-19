# ===============================================================
# main_multinominal.R
#
# Master script for the multinomial simulation study.
# This is the single entry point to reproduce all results
# reported in Section 4.1.2 and Appendix Table 10 of the
# thesis. Running this script in order will:
#
#   1. Load and prepare the four simulated datasets.
#   2. Define model formulas and hyperparameters for each
#      simulation scenario.
#   3. Execute the parallelised 10-fold cross-validation
#      (sourced from Multinominal_analysis/R/analysis.R).
#   4. Generate the summary figure (sourced from Multinominal_analysis/R/image.R).
#   5. Generate the summary tables (sourced from Multinominal_analysis/R/table.R).
#   6. Write the full configuration to a plain-text record.
#   7. Archive all outputs in a timestamped results folder.
#
# Working directory is always the project root thanks to here().
# In RStudio, opening mlml_thesis_research_archive.Rproj sets
# the project root automatically; here() anchors all paths to
# that root regardless of which subfolder the script lives in.
#
# Dependencies: mlml, tidyverse, xtable, gridExtra, here
#   All other packages (nnet, mclogit, rpart, ranger, parallel)
#   are loaded inside Multinominal_analysis/R/analysis.R.
# ===============================================================


# ---------------------------------------------------------------
# Housekeeping
# ---------------------------------------------------------------
rm(list = ls())

library(mlml)
library(tidyverse)
library(xtable)
library(gridExtra)
library(here)


# ===============================================================
# SECTION 1: Load and prepare data
#
# The four CSV files are the output of data_simulation.R
# (stored in Multinominal_analysis/data/). Each corresponds
# to one simulation scenario:
#   s1 — linear / additive   (Scenario 1)
#   s2 — nonlinear           (Scenario 2)
#   s3 — interactions        (Scenario 3)
#   s4 — multicollinearity   (Scenario 4)
# ===============================================================
df_list <- list(
  s1 = read.csv(here("Multinominal_analysis", "data", "simulated_data1.csv")),
  s2 = read.csv(here("Multinominal_analysis", "data", "simulated_data2.csv")),
  s3 = read.csv(here("Multinominal_analysis", "data", "simulated_data3.csv")),
  s4 = read.csv(here("Multinominal_analysis", "data", "simulated_data4.csv"))
)

df_list <- lapply(df_list, function(d) {
  d$id <- as.factor(d$id)
  d$y  <- as.factor(d$y)
  d
})


# ===============================================================
# SECTION 2: Cross-validation and scenario configuration
# ===============================================================
K_folds      <- 10
seed_folds   <- 30
seed_cluster <- 30

x_terms <- paste0("x", 1:10, collapse = " + ")

scenario_cfg <- list(
  s1 = list(
    multinom_formula = as.formula(paste("y ~", x_terms)),
    mblogit_formula  = as.formula(paste("y ~", x_terms)),
    cart_formula     = as.formula(paste("y ~", x_terms)),
    rf_formula       = as.formula(paste("y ~", x_terms)),
    gmert_args = list(
      max_iter_out = 5, max_iter_inn = 50, tol = 1e-4,
      xval = 0, maxdepth = 3, minsplit = 30, minbucket = 10
    ),
    gmerf_args = list(
      max_iter_out = 5, max_iter_inn = 50, tol = 1e-4,
      ntrees = 200, min_node_size = 15
    )
  ),
  s2 = list(
    multinom_formula = as.formula(paste("y ~", x_terms)),
    mblogit_formula  = as.formula(paste("y ~", x_terms)),
    cart_formula     = as.formula(paste("y ~", x_terms)),
    rf_formula       = as.formula(paste("y ~", x_terms)),
    gmert_args = list(
      max_iter_out = 5, max_iter_inn = 50, tol = 1e-4,
      xval = 0, maxdepth = 4, minsplit = 30, minbucket = 10
    ),
    gmerf_args = list(
      max_iter_out = 5, max_iter_inn = 50, tol = 1e-4,
      ntrees = 200, min_node_size = 15
    )
  ),
  s3 = list(
    multinom_formula = as.formula(paste("y ~", x_terms)),
    mblogit_formula  = as.formula(paste("y ~", x_terms)),
    cart_formula     = as.formula(paste("y ~", x_terms)),
    rf_formula       = as.formula(paste("y ~", x_terms)),
    gmert_args = list(
      max_iter_out = 5, max_iter_inn = 50, tol = 1e-4,
      xval = 0, maxdepth = 2, minsplit = 30, minbucket = 10
    ),
    gmerf_args = list(
      max_iter_out = 5, max_iter_inn = 50, tol = 1e-4,
      ntrees = 200, min_node_size = 15
    )
  ),
  s4 = list(
    multinom_formula = as.formula(paste("y ~", x_terms)),
    mblogit_formula  = as.formula(paste("y ~", x_terms)),
    cart_formula     = as.formula(paste("y ~", x_terms)),
    rf_formula       = as.formula(paste("y ~", x_terms)),
    gmert_args = list(
      max_iter_out = 5, max_iter_inn = 50, tol = 1e-4,
      xval = 0, maxdepth = 3, minsplit = 30, minbucket = 10
    ),
    gmerf_args = list(
      max_iter_out = 5, max_iter_inn = 50, tol = 1e-4,
      ntrees = 200, min_node_size = 15
    )
  )
)


# ===============================================================
# SECTION 3: Run analysis, visualization, and tables
# ===============================================================
source(here("Multinominal_analysis", "R", "analysis.R"))
source(here("Multinominal_analysis", "R", "image.R"))
source(here("Multinominal_analysis", "R", "table.R"))


# ===============================================================
# SECTION 4: Write configuration record
# ===============================================================
config_params <- list(
  K_folds      = K_folds,
  seed_folds   = seed_folds,
  seed_cluster = seed_cluster,
  scenario_cfg = scenario_cfg
)
writeLines(
  capture.output(str(config_params)),
  con = here("Multinominal_analysis", "configuration", "configuration.txt")
)


# ===============================================================
# SECTION 5: Archive outputs with timestamp
# ===============================================================
timestamp   <- format(Sys.time(), "%Y%m%d-%H%M%S")
results_dir <- here("Multinominal_analysis", "results", paste0("results_", timestamp))
dir.create(results_dir)

file.copy(here("Multinominal_analysis", "tables", "acc_table.png"),
          file.path(results_dir, "acc_table.png"))
file.copy(here("Multinominal_analysis", "tables", "f1_table.png"),
          file.path(results_dir, "f1_table.png"))
file.copy(here("Multinominal_analysis", "tables", "wf1_table.png"),
          file.path(results_dir, "wf1_table.png"))
file.copy(here("Multinominal_analysis", "images", "graph.png"),
          file.path(results_dir, "graph.png"))
file.copy(here("Multinominal_analysis", "configuration", "configuration.txt"),
          file.path(results_dir, "configuration.txt"))