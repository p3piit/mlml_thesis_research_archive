# ===============================================================
# Multiclass simulation study: parallelized cross-validation
# ===============================================================

rm(list = ls())

library(mlml)
library(tidyverse)
library(xtable)
library(gridExtra)

# ===============================================================
# Read data
df_list <- list(
  s1 = read.csv("data/simulated_data1.csv"),
  s2 = read.csv("data/simulated_data2.csv"),
  s3 = read.csv("data/simulated_data3.csv"),
  s4 = read.csv("data/simulated_data4.csv")
)

# Make sure id is factor and y is factor
df_list <- lapply(df_list, function(d) {
  d$id <- as.factor(d$id)
  d$y  <- as.factor(d$y)
  d
})

# ===============================================================
# Scenario configurations
K_folds <- 10
seed_folds <- 30
seed_cluster <- 30

x_terms <- paste0("x", 1:10, collapse = " + ")

scenario_cfg <- list(
  s1 = list(
    multinom_formula = as.formula(paste("y ~", x_terms)),
    mblogit_formula  = as.formula(paste("y ~", x_terms)),
    cart_formula     = as.formula(paste("y ~", x_terms)),
    rf_formula       = as.formula(paste("y ~", x_terms)),
    gmert_args = list(
      max_iter_out = 5,
      max_iter_inn = 50,
      tol = 1e-4,
      xval = 0,
      maxdepth = 3,
      minsplit = 30,
      minbucket = 10
    ),
    gmerf_args = list(
      max_iter_out = 5,
      max_iter_inn = 50,
      tol = 1e-4,
      ntrees = 200,
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
      tol = 1e-4,
      xval = 0,
      maxdepth = 4,
      minsplit = 30,
      minbucket = 10
    ),
    gmerf_args = list(
      max_iter_out = 5,
      max_iter_inn = 50,
      tol = 1e-4,
      ntrees = 200,
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
      tol = 1e-4,
      xval = 0,
      maxdepth = 2,
      minsplit = 30,
      minbucket = 10
    ),
    gmerf_args = list(
      max_iter_out = 5,
      max_iter_inn = 50,
      tol = 1e-4,
      ntrees = 200,
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
      tol = 1e-4,
      xval = 0,
      maxdepth = 3,
      minsplit = 30,
      minbucket = 10
    ),
    gmerf_args = list(
      max_iter_out = 5,
      max_iter_inn = 50,
      tol = 1e-4,
      ntrees = 200,
      min_node_size = 15
    )
  )
)

# ===============================================================
# Run parallelized cross-validation for all scenarios
source("R/analysis.R")

# ==================================================
# Run image script to create plots
source("R/image.R")

# ==================================================
# Run table script to create tables
source("R/table.R")

# ==================================================
# Write all the parameters setted in this configuration file to a text file for record-keeping
config_params <- list(
  K_folds = K_folds,
  seed_folds = seed_folds,
  seed_cluster = seed_cluster,
  scenario_cfg = scenario_cfg
)
writeLines(capture.output(str(config_params)), con = "configuration/configuration.txt")

# ==================================================
# Save tables, images and the configuration in a folder with timestamp
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
dir.create(paste0("results/results_", timestamp))
file.copy("tables/acc_table.png", paste0("results/results_", timestamp, "/acc_table.png"))
file.copy("tables/f1_table.png", paste0("results/results_", timestamp, "/f1_table.png"))
file.copy("tables/wf1_table.png", paste0("results/results_", timestamp, "/wf1_table.png"))
file.copy("tables/brier_table.png", paste0("results/results_", timestamp, "/brier_table.png"))
file.copy("tables/ce_table.png", paste0("results/results_", timestamp, "/ce_table.png"))
file.copy("images/graph.png", paste0("results/results_", timestamp, "/graph.png"))
file.copy("configuration/configuration.txt", paste0("results/results_", timestamp, "/configuration.txt"))

