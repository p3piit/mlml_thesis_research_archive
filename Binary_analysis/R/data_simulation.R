# ===============================================================
# data_simulation.R
#
# Generates simulated clustered datasets with a binary outcome
# for the binary simulation study (Section 3.2 of the thesis).
#
# The simulation mirrors the design described in Table 2:
# 30 clusters, 60 observations each, a random intercept and
# a random slope on x1, and a binary outcome drawn from a
# Bernoulli logistic model.
#
# Four scenarios are implemented, each varying one aspect of
# the data-generating fixed component while keeping the
# cluster structure and random-effects covariance constant:
#   1. Linear / additive      — correctly specified baseline
#   2. Nonlinear              — log transformation of x2
#   3. Interactions           — x2 × x4 interaction term
#   4. Multicollinearity      — strongly correlated predictors
#
# The fixed-effects specifications match Appendix Table 5
# of the thesis.
#
# Output: four CSV files saved to data/
#   simulated_data.csv  (scenario 1)
#   simulated_data.csv  (scenario 2)
#   simulated_data.csv  (scenario 3)
#   simulated_data.csv  (scenario 4)
#
# Dependencies: MASS (mvrnorm)
# ===============================================================


# ---------------------------------------------------------------
# Scenario 1: linear / additive  (Appendix Table 5, "Baseline")
#
# The fixed component is a standard linear logistic predictor.
# This is the correctly specified scenario for GLM and GLMM,
# and serves as the benchmark against which the other scenarios
# are evaluated.
#
# Predictors:
#   x1  Uniform(-2, 2)   — also carries the random slope
#   x2  Normal(0, 1)     — standard normal
#   x3  Bernoulli(0.5)   — binary covariate
#
# Random-effects structure (Table 2):
#   (b0_j, b1_j) ~ N(0, D),  D = [[0.64, 0.08], [0.08, 0.25]]
#   which corresponds to sigma_b0 = 0.8, sigma_b1 = 0.5,
#   rho = 0.2 (correlation between intercept and slope).
#
# Fixed component (Appendix Table 5):
#   eta = beta0 + beta1*x1 + beta2*x2 + beta3*x3
#       = 0.5   + 1.2*x1   - 0.8*x2   + 0.6*x3
# ---------------------------------------------------------------
sim_data1 <- function(G = 50,           # number of clusters
                      n_i = 40,         # observations per cluster
                      beta0 = 0.5,      # fixed intercept
                      beta1 = 1.2,      # fixed effect for x1
                      beta2 = -0.8,     # fixed effect for x2
                      beta3 = 0.6,      # fixed effect for x3
                      sigma_b0 = 0.8,   # SD of random intercept
                      sigma_b1 = 0.5,   # SD of random slope for x1
                      rho = 0.2,        # correlation between b0 and b1
                      seed = 123) {     # random seed for reproducibility
  
  set.seed(seed)
  
  # Covariance matrix of random effects (b0_j, b1_j); Table 2
  D <- matrix(c(sigma_b0^2,             rho * sigma_b0 * sigma_b1,
                rho * sigma_b0 * sigma_b1, sigma_b1^2), 2, 2)
  
  # Draw one (b0_j, b1_j) pair per cluster
  b <- MASS::mvrnorm(G, mu = c(0, 0), Sigma = D)
  
  # Cluster membership vector
  id <- rep(1:G, each = n_i)
  
  # Generate covariates
  x1 <- runif(G * n_i, -2, 2)      # uniform; also the random-slope variable
  x2 <- rnorm(G * n_i, 0, 1)       # standard normal
  x3 <- rbinom(G * n_i, 1, 0.5)    # binary
  
  # Linear predictor: fixed part + cluster-specific random intercept
  # and random slope on x1
  eta <- numeric(G * n_i)
  for (g in seq_len(G)) {
    idx <- which(id == g)
    eta[idx] <- beta0 + beta1 * x1[idx] + beta2 * x2[idx] + beta3 * x3[idx] +
      b[g, 1] + b[g, 2] * x1[idx]
  }
  
  # Convert linear predictor to probability via logistic link
  p <- 1 / (1 + exp(-eta))
  
  # Draw binary outcome from Bernoulli(p)
  y <- rbinom(G * n_i, 1, p)
  
  data.frame(
    id = factor(id),
    y  = y,
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
}


# ---------------------------------------------------------------
# Scenario 2: nonlinear  (Appendix Table 5, "Nonlinear")
#
# The fixed component replaces the linear x2 term with
# beta2 * log(x2). To ensure log() is well-defined, x2 is
# drawn from Normal(5, 1) so that all values are strictly
# positive with overwhelming probability.
#
# This tests whether tree-based methods that can approximate
# nonlinear functions outperform the misspecified linear GLM
# and GLMM (which still receive x2 untransformed).
#
# Predictors:
#   x1  Uniform(-2, 2)
#   x2  Normal(5, 1)    — positive mean ensures log(x2) is safe
#   x3  Bernoulli(0.5)
#
# Fixed component (Appendix Table 5):
#   eta = 0.5 + 1.2*x1 - 1.0*log(x2) + 0.6*x3
# ---------------------------------------------------------------
sim_data2 <- function(G = 50,
                      n_i = 40,
                      beta0 = 0.5,
                      beta1 = 1.2,
                      beta2 = -1,       # coefficient on log(x2)
                      beta3 = 0.6,
                      sigma_b0 = 0.8,
                      sigma_b1 = 0.5,
                      rho = 0.2,
                      seed = 123) {
  
  set.seed(seed)
  
  D <- matrix(c(sigma_b0^2,             rho * sigma_b0 * sigma_b1,
                rho * sigma_b0 * sigma_b1, sigma_b1^2), 2, 2)
  
  b <- MASS::mvrnorm(G, mu = c(0, 0), Sigma = D)
  
  id <- rep(1:G, each = n_i)
  
  x1 <- runif(G * n_i, -2, 2)
  x2 <- rnorm(G * n_i, 5, 1)    # shifted mean ensures log(x2) > 0
  x3 <- rbinom(G * n_i, 1, 0.5)
  
  # Fixed component includes log(x2) rather than x2 directly
  eta <- numeric(G * n_i)
  for (g in seq_len(G)) {
    idx <- which(id == g)
    eta[idx] <- beta0 + beta1 * x1[idx] + beta2 * log(x2[idx]) + beta3 * x3[idx] +
      b[g, 1] + b[g, 2] * x1[idx]
  }
  
  p <- 1 / (1 + exp(-eta))
  y <- rbinom(G * n_i, 1, p)
  
  data.frame(
    id  = factor(id),
    y   = y,
    x1  = x1,
    x2  = x2,
    x3  = x3,
    eta = eta   # retained for diagnostic use (e.g. checking DGP)
  )
}


# ---------------------------------------------------------------
# Scenario 3: interactions  (Appendix Table 5, "Interaction")
#
# The x2 main effect is replaced by the product x2 * x4, where
# x4 is an additional binary covariate. This introduces a
# two-way interaction that misspecifies the linear GLM and GLMM
# (which receive only main effects) and tests whether tree-based
# methods recover the interaction structure through recursive
# splitting.
#
# Note: CART and RF in the CV loop use only x1, x2, x3 as
# predictors (x4 is excluded from their formula). This is a
# known feature of the script as written — the interaction term
# is present in the DGP but x4 is not passed to the tree-based
# models, making the scenario a test of implicit vs explicit
# interaction recovery.
#
# Predictors:
#   x1  Uniform(-2, 2)
#   x2  Normal(5, 1)
#   x3  Bernoulli(0.5)
#   x4  Bernoulli(0.5)   — additional binary covariate
#
# Fixed component (Appendix Table 5):
#   eta = 0.5 + 1.2*x1 - 0.8*(x2*x4) + 0.6*x3
# ---------------------------------------------------------------
sim_data3 <- function(G = 50,
                      n_i = 40,
                      beta0 = 0.5,
                      beta1 = 1.2,
                      beta2 = -0.8,     # coefficient on the x2*x4 interaction
                      beta3 = 0.6,
                      sigma_b0 = 0.8,
                      sigma_b1 = 0.5,
                      rho = 0.2,
                      seed = 123) {
  
  set.seed(seed)
  
  D <- matrix(c(sigma_b0^2,             rho * sigma_b0 * sigma_b1,
                rho * sigma_b0 * sigma_b1, sigma_b1^2), 2, 2)
  
  b <- MASS::mvrnorm(G, mu = c(0, 0), Sigma = D)
  
  id <- rep(1:G, each = n_i)
  
  x1 <- runif(G * n_i, -2, 2)
  x2 <- rnorm(G * n_i, 5, 1)
  x3 <- rbinom(G * n_i, 1, 0.5)
  x4 <- rbinom(G * n_i, 1, 0.5)   # extra binary covariate for interaction
  
  # Fixed component: x2 enters only through the x2*x4 product
  eta <- numeric(G * n_i)
  for (g in seq_len(G)) {
    idx <- which(id == g)
    eta[idx] <- beta0 + beta1 * x1[idx] + beta2 * x2[idx] * x4[idx] + beta3 * x3[idx] +
      b[g, 1] + b[g, 2] * x1[idx]
  }
  
  p <- 1 / (1 + exp(-eta))
  y <- rbinom(G * n_i, 1, p)
  
  data.frame(
    id = factor(id),
    y  = y,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    x4 = x4
  )
}


# ---------------------------------------------------------------
# Scenario 4: multicollinearity  (Appendix Table 5, "Multicollinearity")
#
# The fixed component is identical to Scenario 1 (linear), but
# the three predictors are drawn from a multivariate normal
# distribution with strong pairwise correlations:
#   Cor(x1, x2) = 0.9,  Cor(x1, x3*) = 0.7,  Cor(x2, x3*) = 0.6
# where x3* is the latent continuous variable thresholded to
# produce binary x3 = I(x3* > 0).
#
# This tests whether models degrade differently under near-
# collinearity. Because the fixed structure is linear, the
# GLMM is correctly specified but its coefficient estimates
# become unstable; tree-based methods may be more robust.
#
# Predictor covariance (Table 2):
#   Sigma_X = [[1, 0.9, 0.7], [0.9, 1, 0.6], [0.7, 0.6, 1]]
#
# Fixed component (Appendix Table 5):
#   eta = 0.5 + 1.2*x1 - 0.8*x2 + 0.6*x3   (same as Scenario 1)
# ---------------------------------------------------------------
sim_data4 <- function(G = 50,
                      n_i = 40,
                      beta0 = 0.5,
                      beta1 = 1.2,
                      beta2 = -0.8,
                      beta3 = 0.6,
                      sigma_b0 = 0.8,
                      sigma_b1 = 0.5,
                      rho = 0.2,
                      seed = 123) {
  
  set.seed(seed)
  
  D <- matrix(c(sigma_b0^2,             rho * sigma_b0 * sigma_b1,
                rho * sigma_b0 * sigma_b1, sigma_b1^2), 2, 2)
  
  b <- MASS::mvrnorm(G, mu = c(0, 0), Sigma = D)
  
  id <- rep(1:G, each = n_i)
  
  N <- G * n_i
  
  # Strongly correlated predictor covariance matrix (Table 2)
  Sigma_X <- matrix(c(1,   0.9, 0.7,
                      0.9, 1,   0.6,
                      0.7, 0.6, 1),   ncol = 3)
  X  <- MASS::mvrnorm(N, mu = c(0, 0, 0), Sigma = Sigma_X)
  x1 <- X[, 1]
  x2 <- X[, 2]
  x3 <- ifelse(X[, 3] > 0, 1, 0)   # threshold latent variable to binary
  
  # Fixed component: same linear structure as Scenario 1
  eta <- numeric(N)
  for (g in seq_len(G)) {
    idx <- which(id == g)
    eta[idx] <- beta0 + beta1 * x1[idx] + beta2 * x2[idx] + beta3 * x3[idx] +
      b[g, 1] + b[g, 2] * x1[idx]
  }
  
  p <- 1 / (1 + exp(-eta))
  y <- rbinom(N, 1, p)
  
  data.frame(
    id = factor(id),
    y  = y,
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
}


# ---------------------------------------------------------------
# Generate and save all four scenario datasets
#
# Default arguments reproduce the design in Table 2:
# G = 30 clusters, n_i = 60 observations each (N = 1800).
# seed = 123 is set inside each function for reproducibility
# of both predictor generation and outcome sampling.
#
# Output files are written to data/ relative to the working
# directory. The directory must exist before running this
# script; it is created by the folder structure of the
# research archive.
# ---------------------------------------------------------------
df_bin1 <- sim_data1(G = 30, n_i = 60)
df_bin2 <- sim_data2(G = 30, n_i = 60)
df_bin3 <- sim_data3(G = 30, n_i = 60)
df_bin4 <- sim_data4(G = 30, n_i = 60)

write.csv(df_bin1, here("Binary_analysis", "data", "simulated_data1.csv"), row.names = FALSE)
write.csv(df_bin2, here("Binary_analysis", "data", "simulated_data2.csv"), row.names = FALSE)
write.csv(df_bin3, here("Binary_analysis", "data", "simulated_data3.csv"), row.names = FALSE)
write.csv(df_bin4, here("Binary_analysis", "data", "simulated_data4.csv"), row.names = FALSE)


