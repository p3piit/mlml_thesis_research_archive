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
#      (sourced from R/analysis.R).
#   4. Generate the summary figure (sourced from R/image.R).
#   5. Generate the summary tables (sourced from R/table.R).
#   6. Write the full configuration to a plain-text record.
#   7. Archive all outputs in a timestamped results folder.
#
# Working directory must be the project root (i.e. the folder
# containing the R/, data/, images/, tables/, configuration/,
# and results/ subdirectories). In RStudio, opening the
# .Rproj file sets this automatically.
#
# Dependencies: mlml, tidyverse, xtable, gridExtra
#   All other packages (nnet, mclogit, rpart, ranger,
#   parallel) are loaded inside R/analysis.R.
# ===============================================================


# ---------------------------------------------------------------
# Housekeeping
# Clear the workspace before running to ensure no objects from
# a previous session interfere with the current analysis.
# ---------------------------------------------------------------
rm(list = ls())

library(mlml)
library(tidyverse)
library(xtable)
library(gridExtra)


# ===============================================================
# SECTION 1: Load and prepare data
#
# The four CSV files are the output of generate_categorical_data.R
# (stored in data/). Each corresponds to one simulation scenario:
#   s1 — linear / additive   (Scenario 1)
#   s2 — nonlinear           (Scenario 2)
#   s3 — interactions        (Scenario 3)
#   s4 — multicollinearity   (Scenario 4)
#
# id must be a factor so that cluster-aware fold construction
# and mixed-effects model fitting treat it as a grouping
# variable rather than a numeric covariate.
# y must be a factor so that all classifiers recognise the
# task as categorical and return probability matrices with
# consistently named columns across folds.
# ===============================================================
df_list <- list(
  s1 = read.csv("data/simulated_data1.csv"),
  s2 = read.csv("data/simulated_data2.csv"),
  s3 = read.csv("data/simulated_data3.csv"),
  s4 = read.csv("data/simulated_data4.csv")
)

df_list <- lapply(df_list, function(d) {
  d$id <- as.factor(d$id)
  d$y  <- as.factor(d$y)
  d
})


# ===============================================================
# SECTION 2: Cross-validation and scenario configuration
#
# Global CV settings shared across all four scenarios.
# Seeds are fixed here so that fold construction and the
# parallel RNG stream are fully reproducible; the same values
# must be used if results are to be exactly replicated.
# ===============================================================
K_folds      <- 10   # number of cross-validation folds
seed_folds   <- 30   # seed for make_cluster_folds()
seed_cluster <- 30   # seed for clusterSetRNGStream()

# Formula shared by all parametric and tree-based models:
# y regressed on all ten predictors x1 … x10
x_terms <- paste0("x", 1:10, collapse = " + ")

# ---------------------------------------------------------------
# scenario_cfg: one list entry per scenario, each containing
# model formulas and hyperparameter lists for GMERT and GMERF.
#
# All four scenarios use the same formula (all ten predictors)
# and the same GMERF hyperparameters. GMERT maxdepth differs
# across scenarios to reflect the varying complexity of the
# data-generating fixed component:
#   s1 (linear)         maxdepth = 3  — shallow tree sufficient
#   s2 (nonlinear)      maxdepth = 4  — extra depth for smooth
#                                       nonlinear surfaces
#   s3 (interactions)   maxdepth = 2  — interactions are
#                                       captured at depth 2;
#                                       deeper trees overfit
#   s4 (multicollinearity) maxdepth = 3 — same as s1 since the
#                                       fixed component is linear
#
# GMERT hyperparameters (passed to mlml::fit_gmert_cat()):
#   max_iter_out  Maximum outer PQL iterations (convergence of
#                 the linear predictor; Algorithm 1 outer loop).
#   max_iter_inn  Maximum inner EM iterations per outer step
#                 (Algorithm 1 inner loop).
#   tol           Convergence tolerance for the GLL criterion.
#   xval          Number of cross-validations for rpart pruning;
#                 0 disables pruning (cp = 0 in rpart).
#   maxdepth      Maximum tree depth for the fixed-effects tree.
#   minsplit      Minimum observations to attempt a node split.
#   minbucket     Minimum observations in any terminal node.
#
# GMERF hyperparameters (passed to mlml::fit_gmerf_cat()):
#   max_iter_out, max_iter_inn, tol — as above.
#   ntrees        Number of trees in the random forest fixed-
#                 effects component.
#   min_node_size Minimum terminal node size for the forest.
# ---------------------------------------------------------------
scenario_cfg <- list(
  s1 = list(
    multinom_formula = as.formula(paste("y ~", x_terms)),
    mblogit_formula  = as.formula(paste("y ~", x_terms)),
    cart_formula     = as.formula(paste("y ~", x_terms)),
    rf_formula       = as.formula(paste("y ~", x_terms)),
    gmert_args = list(
      max_iter_out = 5,
      max_iter_inn = 50,
      tol          = 1e-4,
      xval         = 0,
      maxdepth     = 3,   # linear DGP: moderate depth
      minsplit     = 30,
      minbucket    = 10
    ),
    gmerf_args = list(
      max_iter_out  = 5,
      max_iter_inn  = 50,
      tol           = 1e-4,
      ntrees        = 200,
      min_node_size = 15
    )
  ),
  s2 = list(
    multinom_formula = as.formula(paste("y ~", x_terms)),
    mblogit_formula  = as.formula(paste("y ~", x_terms)),
    cart_formula     = as.formula(paste("y ~", x_terms)),
    rf_formula       = as.formula(paste("y ~", x_terms)),
    gmert_args = list(
      max_iter_out = 5,
      max_iter_inn = 50,
      tol          = 1e-4,
      xval         = 0,
      maxdepth     = 4,   # nonlinear DGP: extra depth needed
      minsplit     = 30,
      minbucket    = 10
    ),
    gmerf_args = list(
      max_iter_out  = 5,
      max_iter_inn  = 50,
      tol           = 1e-4,
      ntrees        = 200,
      min_node_size = 15
    )
  ),
  s3 = list(
    multinom_formula = as.formula(paste("y ~", x_terms)),
    mblogit_formula  = as.formula(paste("y ~", x_terms)),
    cart_formula     = as.formula(paste("y ~", x_terms)),
    rf_formula       = as.formula(paste("y ~", x_terms)),
    gmert_args = list(
      max_iter_out = 5,
      max_iter_inn = 50,
      tol          = 1e-4,
      xval         = 0,
      maxdepth     = 2,   # interaction DGP: shallow to avoid
      # overfitting pure interaction structure
      minsplit     = 30,
      minbucket    = 10
    ),
    gmerf_args = list(
      max_iter_out  = 5,
      max_iter_inn  = 50,
      tol           = 1e-4,
      ntrees        = 200,
      min_node_size = 15
    )
  ),
  s4 = list(
    multinom_formula = as.formula(paste("y ~", x_terms)),
    mblogit_formula  = as.formula(paste("y ~", x_terms)),
    cart_formula     = as.formula(paste("y ~", x_terms)),
    rf_formula       = as.formula(paste("y ~", x_terms)),
    gmert_args = list(
      max_iter_out = 5,
      max_iter_inn = 50,
      tol          = 1e-4,
      xval         = 0,
      maxdepth     = 3,   # multicollinearity DGP: same linear
      # structure as s1, same depth
      minsplit     = 30,
      minbucket    = 10
    ),
    gmerf_args = list(
      max_iter_out  = 5,
      max_iter_inn  = 50,
      tol           = 1e-4,
      ntrees        = 200,
      min_node_size = 15
    )
  )
)


# ===============================================================
# SECTION 3: Run analysis, visualization, and tables
#
# Each source() call expects the objects defined above to be
# present in the global environment:
#   R/analysis.R  needs df_list, scenario_cfg, K_folds,
#                 seed_folds, seed_cluster; produces `results`
#   R/image.R     needs results; produces images/graph.png
#   R/table.R     needs results; produces tables/*.png and
#                 prints LaTeX to the console
# ===============================================================
source("R/analysis.R")
source("R/image.R")
source("R/table.R")


# ===============================================================
# SECTION 4: Write configuration record
#
# Saves a plain-text snapshot of all analysis parameters to
# configuration/configuration.txt. This file is included in
# the timestamped results archive (Section 5) and in the
# research archive, ensuring that the exact settings used to
# produce any given set of outputs are always recoverable.
#
# str() produces a compact, human-readable representation of
# nested lists, which is more legible than dput() for this
# purpose.
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
# all outputs into it. This means re-running the script never
# overwrites previous results — each run produces its own
# dated folder, making it easy to track how outputs changed
# across analysis iterations.
#
# ===============================================================
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
dir.create(paste0("results/results_", timestamp))

file.copy("tables/acc_table.png",
          paste0("results/results_", timestamp, "/acc_table.png"))
file.copy("tables/f1_table.png",
          paste0("results/results_", timestamp, "/f1_table.png"))
file.copy("tables/wf1_table.png",
          paste0("results/results_", timestamp, "/wf1_table.png"))
file.copy("images/graph.png",
          paste0("results/results_", timestamp, "/graph.png"))
file.copy("configuration/configuration.txt",
          paste0("results/results_", timestamp, "/configuration.txt"))