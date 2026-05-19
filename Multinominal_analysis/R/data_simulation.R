# ===============================================================
# data_simulation.R
#
# Generates simulated clustered datasets with a 5-class
# categorical outcome for the multinomial simulation study.
#
# The simulation mirrors the design described in Section 3.2
# of the thesis (Table 3): 100 clusters, 18 observations each,
# 10 multivariate-normal predictors, a random-intercept-only
# structure, and class 5 as the reference category.
#
# Four scenarios are implemented, each varying the functional
# form of the fixed component while keeping the cluster
# structure and predictor distribution constant:
#   1. Linear / additive      (rho = 0.20)
#   2. Nonlinear              (rho = 0.20)
#   3. Interactions           (rho = 0.20)
#   4. Multicollinearity      (rho = 0.85)
#
# The exact fixed-effects specifications match Appendix Table 6
# of the thesis.
#
# Output: four CSV files saved to data/
#   simulated_data1.csv  (scenario 1)
#   simulated_data2.csv  (scenario 2)
#   simulated_data3.csv  (scenario 3)
#   simulated_data4.csv  (scenario 4)
#
# Dependencies: MASS (mvrnorm)
# ===============================================================


# ---------------------------------------------------------------
# sim_cat5_core()
#
# Shared engine called by all four scenario functions.
#
# Given a pre-built predictor matrix X and a matrix of
# fixed-part log-odds eta_fixed (one column per non-reference
# class), it:
#   (1) draws group-specific random intercepts from a
#       multivariate normal distribution (one intercept per
#       log-odds equation, i.e. 4-dimensional);
#   (2) adds those intercepts to the fixed-part log-odds to
#       obtain the full linear predictor eta;
#   (3) applies the softmax transformation with class 5 as
#       reference to obtain class probabilities P;
#   (4) samples the observed outcome y by drawing from the
#       resulting categorical distribution.
#
# Arguments:
#   X            N x 10 predictor matrix (numeric).
#   id           Length-N vector of cluster identifiers.
#   eta_fixed    N x 4 matrix of fixed log-odds
#                (columns = logit 1 … logit 4).
#   intercept_sd Scalar SD for independent random intercepts
#                when D_intercept is not supplied. Default 0.8.
#   D_intercept  Optional 4 x 4 covariance matrix for the
#                random intercepts. When NULL, a diagonal
#                matrix with variance intercept_sd^2 is used,
#                implying independent intercepts across logits.
#   seed         Integer seed passed to set.seed(). Ignored
#                when NULL.
#   return_prob  Logical. If TRUE, the five softmax
#                probabilities are appended to the returned
#                data frame as columns prob_1 … prob_5.
#
# Returns:
#   A data frame with columns:
#     id  – cluster identifier (factor)
#     y   – outcome factor with levels class1 … class5
#     x1 … x10 – predictors
#     prob_1 … prob_5 (only when return_prob = TRUE)
# ---------------------------------------------------------------
sim_cat5_core <- function(
    X,                        # N x 10 predictor matrix
    id,                       # cluster ids
    eta_fixed,                # N x 4 fixed-part logits
    intercept_sd = 0.8,       # SD of random intercepts
    D_intercept = NULL,       # optional 4 x 4 covariance across logits
    seed = NULL,
    return_prob = FALSE
) {
  if (!is.null(seed)) set.seed(seed)
  
  X <- as.matrix(X)
  eta_fixed <- as.matrix(eta_fixed)
  
  N <- nrow(X)
  p <- ncol(X)
  K <- 5    # total number of outcome classes
  K1 <- 4   # number of non-reference classes (= K - 1)
  
  # Convert cluster ids to consecutive integers 1 … G
  groups <- as.integer(as.factor(id))
  G <- length(unique(groups))
  
  colnames(X) <- paste0("x", seq_len(p))
  
  # ------------------------------------------------------------------
  # Random-intercept covariance matrix (4 x 4)
  #
  # Each of the K-1 = 4 log-odds equations gets its own random
  # intercept. When D_intercept is NULL, the intercepts are
  # independent with variance intercept_sd^2 (diagonal D).
  # A user-supplied D_intercept allows for correlated intercepts
  # across logit equations, which would capture, for example,
  # cluster-level tendencies that affect several classes at once.
  # ------------------------------------------------------------------
  if (is.null(D_intercept)) {
    D_intercept <- diag(intercept_sd^2, K1)
  } else {
    D_intercept <- as.matrix(D_intercept)
    if (!all(dim(D_intercept) == c(K1, K1))) {
      stop("D_intercept must be a 4 x 4 matrix.")
    }
  }
  
  # Draw one 4-dimensional random intercept per cluster:
  # b0 is a G x 4 matrix; row g = (b0_g^(1), …, b0_g^(4))
  b0 <- MASS::mvrnorm(n = G, mu = rep(0, K1), Sigma = D_intercept)
  
  # ------------------------------------------------------------------
  # Add cluster-specific random intercepts to the fixed log-odds
  #
  # For observation i in cluster g, the full linear predictor is
  #   eta[i, k] = eta_fixed[i, k] + b0[g, k],  k = 1 … 4
  # sweep() adds b0[g, ] to every row of eta_fixed that belongs
  # to cluster g.
  # ------------------------------------------------------------------
  eta <- eta_fixed
  for (g in seq_len(G)) {
    idx <- which(groups == g)
    eta[idx, ] <- sweep(eta[idx, , drop = FALSE], 2, b0[g, ], FUN = "+")
  }
  
  # ------------------------------------------------------------------
  # Softmax transformation (class 5 = reference category)
  #
  # The multinomial logit model with K = 5 classes and reference K is:
  #   P(y = k) = exp(eta_k) / (1 + sum_{j=1}^{4} exp(eta_j)),  k = 1…4
  #   P(y = 5) = 1           / (1 + sum_{j=1}^{4} exp(eta_j))
  #
  # This matches the link function described in Section 3.1.1 and
  # Appendix B of the thesis.
  # ------------------------------------------------------------------
  exp_eta <- exp(eta)
  denom <- 1 + rowSums(exp_eta)
  
  P <- matrix(0, N, K)
  P[, 1:K1] <- exp_eta / denom       # classes 1–4
  P[, K] <- 1 / denom                # class 5 (reference)
  
  # Draw the observed categorical outcome for each observation
  # by sampling from its individual probability vector
  y_int <- apply(P, 1, function(pv) sample.int(K, size = 1L, prob = pv))
  y <- factor(y_int, levels = 1:5, labels = paste0("class", 1:5))
  
  # Assemble output data frame
  df <- data.frame(
    id = factor(groups),
    y = y,
    X,
    stringsAsFactors = FALSE
  )
  
  # Optionally append the five class probabilities (useful for
  # diagnostics or computing the oracle Bayes error rate)
  if (return_prob) {
    colnames(P) <- paste0("prob_", 1:5)
    df <- cbind(df, as.data.frame(P))
  }
  
  df
}


# ---------------------------------------------------------------
# Scenario 1: linear / additive  (Appendix Table 6, "Baseline")
#
# All four log-odds are strictly linear combinations of the
# predictors — no transformations, interactions, or collinearity.
# This is the correctly specified scenario for GLM and GLMM,
# and serves as the benchmark against which nonlinearity,
# interactions, and multicollinearity are evaluated.
#
# Predictor correlation: rho = 0.20 (mild).
# Random-intercept SD: 5 (large relative to the fixed effects,
#   so the clustering structure is prominent).
#
# Fixed components (from Appendix Table 6):
#   eta1 =  0.8 + 1.0*x1 - 0.8*x2 + 0.6*x3 + 0.5*x4 - 0.4*x5
#   eta2 = -0.4 - 0.7*x1 + 0.9*x2 - 0.5*x3 + 0.4*x6 + 0.3*x7
#   eta3 =  0.2 + 0.6*x2 - 0.6*x4 + 0.7*x8 - 0.5*x9
#   eta4 = -0.6 + 0.5*x1 + 0.4*x5 - 0.8*x7 + 0.6*x10
# ---------------------------------------------------------------
sim_cat5_data1 <- function(
    G = 100,
    n_i = 18,
    rho = 0.2,
    intercept_sd = 5,
    seed = 123,
    return_prob = FALSE
) {
  set.seed(seed)
  
  N <- G * n_i   # total number of observations (1800 by default)
  p <- 10        # number of predictors
  id <- rep(seq_len(G), each = n_i)   # cluster membership vector
  
  # Equicorrelation covariance matrix for the predictors:
  # all off-diagonal entries equal rho, diagonal entries equal 1
  Sigma_X <- matrix(rho, p, p)
  diag(Sigma_X) <- 1
  X <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = Sigma_X)
  
  # Fixed log-odds: linear functions of the predictors
  # Each eta corresponds to one non-reference class logit
  eta1 <-  0.8 + 1.0 * X[,1] - 0.8 * X[,2] + 0.6 * X[,3] + 0.5 * X[,4] - 0.4 * X[,5]
  eta2 <- -0.4 - 0.7 * X[,1] + 0.9 * X[,2] - 0.5 * X[,3] + 0.4 * X[,6] + 0.3 * X[,7]
  eta3 <-  0.2 + 0.6 * X[,2] - 0.6 * X[,4] + 0.7 * X[,8] - 0.5 * X[,9]
  eta4 <- -0.6 + 0.5 * X[,1] + 0.4 * X[,5] - 0.8 * X[,7] + 0.6 * X[,10]
  
  eta_fixed <- cbind(eta1, eta2, eta3, eta4)
  
  sim_cat5_core(
    X = X,
    id = id,
    eta_fixed = eta_fixed,
    intercept_sd = intercept_sd,
    return_prob = return_prob
  )
}


# ---------------------------------------------------------------
# Scenario 2: nonlinear  (Appendix Table 6, "Nonlinear")
#
# The fixed log-odds contain smooth nonlinear transformations
# of the predictors: logarithm, square, square root, sine,
# cosine, and hyperbolic tangent. This tests whether machine-
# learning methods that can capture nonlinearity (CART, RF,
# GMERT, GMERF) outperform parametric models (GLM, GLMM) that
# assume a linear fixed component.
#
# Predictor correlation: rho = 0.20 (same as Scenario 1).
#
# Fixed components (from Appendix Table 6):
#   eta1 =  0.8 + 1.0*x1 - 0.8*log(|x2|+1) + 0.7*x3^2 - 0.5*sin(x4)
#   eta2 = -0.4 - 0.7*x1 + 0.9*sqrt(|x2|+0.5) - 0.6*x5^2 + 0.4*cos(x3)
#   eta3 =  0.2 + 0.6*tanh(x6) - 0.7*log(|x7|+1) + 0.5*x8^2
#   eta4 = -0.6 + 0.8*sin(x1) + 0.5*x9^2 - 0.6*sqrt(|x10|+0.5)
#
# Note: log(|x|+1) and sqrt(|x|+c) are used instead of log(x)
# and sqrt(x) to handle negative predictor values that arise
# from the multivariate-normal distribution.
# ---------------------------------------------------------------
sim_cat5_data2 <- function(
    G = 100,
    n_i = 18,
    rho = 0.2,
    intercept_sd = 5,
    seed = 123,
    return_prob = FALSE
) {
  set.seed(seed)
  
  N <- G * n_i
  p <- 10
  id <- rep(seq_len(G), each = n_i)
  
  Sigma_X <- matrix(rho, p, p)
  diag(Sigma_X) <- 1
  X <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = Sigma_X)
  
  # Fixed log-odds: nonlinear transformations of the predictors
  eta1 <-  0.8 + 1.0 * X[,1] - 0.8 * log(abs(X[,2]) + 1) + 0.7 * X[,3]^2 - 0.5 * sin(X[,4])
  eta2 <- -0.4 - 0.7 * X[,1] + 0.9 * sqrt(abs(X[,2]) + 0.5) - 0.6 * X[,5]^2 + 0.4 * cos(X[,3])
  eta3 <-  0.2 + 0.6 * tanh(X[,6]) - 0.7 * log(abs(X[,7]) + 1) + 0.5 * X[,8]^2
  eta4 <- -0.6 + 0.8 * sin(X[,1]) + 0.5 * X[,9]^2 - 0.6 * sqrt(abs(X[,10]) + 0.5)
  
  eta_fixed <- cbind(eta1, eta2, eta3, eta4)
  
  sim_cat5_core(
    X = X,
    id = id,
    eta_fixed = eta_fixed,
    intercept_sd = intercept_sd,
    return_prob = return_prob
  )
}


# ---------------------------------------------------------------
# Scenario 3: interactions  (Appendix Table 6, "Interaction")
#
# The fixed log-odds are driven primarily by two-way products
# of predictors, with smaller main effects. This tests whether
# tree-based methods that implicitly capture interaction
# structure (through recursive splits) outperform additive
# parametric models on interaction-heavy data.
#
# Predictor correlation: rho = 0.20 (same as Scenarios 1–2).
#
# Fixed components (from Appendix Table 6):
#   eta1 =  0.8 + 0.8*x1 - 0.5*x2 + 1.2*x2*x4 - 0.7*x3*x5
#   eta2 = -0.4 - 0.6*x1 + 0.7*x3 - 1.0*x1*x6 + 0.8*x4*x5
#   eta3 =  0.2 + 0.9*x2*x7 - 0.6*x8 + 0.7*x9*x10
#   eta4 = -0.6 + 0.5*x1 + 0.6*x5*x6 - 0.8*x7*x8
# ---------------------------------------------------------------
sim_cat5_data3 <- function(
    G = 100,
    n_i = 18,
    rho = 0.2,
    intercept_sd = 5,
    seed = 123,
    return_prob = FALSE
) {
  set.seed(seed)
  
  N <- G * n_i
  p <- 10
  id <- rep(seq_len(G), each = n_i)
  
  Sigma_X <- matrix(rho, p, p)
  diag(Sigma_X) <- 1
  X <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = Sigma_X)
  
  # Fixed log-odds: two-way interaction terms between predictors
  eta1 <-  0.8 + 0.8 * X[,1] - 0.5 * X[,2] + 1.2 * X[,2] * X[,4] - 0.7 * X[,3] * X[,5]
  eta2 <- -0.4 - 0.6 * X[,1] + 0.7 * X[,3] - 1.0 * X[,1] * X[,6] + 0.8 * X[,4] * X[,5]
  eta3 <-  0.2 + 0.9 * X[,2] * X[,7] - 0.6 * X[,8] + 0.7 * X[,9] * X[,10]
  eta4 <- -0.6 + 0.5 * X[,1] + 0.6 * X[,5] * X[,6] - 0.8 * X[,7] * X[,8]
  
  eta_fixed <- cbind(eta1, eta2, eta3, eta4)
  
  sim_cat5_core(
    X = X,
    id = id,
    eta_fixed = eta_fixed,
    intercept_sd = intercept_sd,
    return_prob = return_prob
  )
}


# ---------------------------------------------------------------
# Scenario 4: multicollinearity  (Appendix Table 6, "Multicollinearity")
#
# The fixed-effects specification is identical to Scenario 1
# (linear / additive), but the predictor covariance matrix
# uses rho = 0.85 instead of 0.20. This creates strong
# pairwise correlations among all predictors and evaluates
# how regularization implicit in tree-based methods compares
# with parametric approaches under near-collinearity.
#
# Fixed components: same as Scenario 1 (linear).
# Predictor correlation: rho = 0.85 (strong multicollinearity).
# ---------------------------------------------------------------
sim_cat5_data4 <- function(
    G = 100,
    n_i = 18,
    rho = 0.85,          # high equicorrelation among predictors
    intercept_sd = 5,
    seed = 123,
    return_prob = FALSE
) {
  set.seed(seed)
  
  N <- G * n_i
  p <- 10
  id <- rep(seq_len(G), each = n_i)
  
  # High equicorrelation: rho = 0.85 for all predictor pairs
  Sigma_X <- matrix(rho, p, p)
  diag(Sigma_X) <- 1
  X <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = Sigma_X)
  
  # Fixed log-odds: same linear structure as Scenario 1;
  # only the predictor correlation structure differs
  eta1 <-  0.8 + 1.0 * X[,1] - 0.8 * X[,2] + 0.6 * X[,3] + 0.5 * X[,4] - 0.4 * X[,5]
  eta2 <- -0.4 - 0.7 * X[,1] + 0.9 * X[,2] - 0.5 * X[,3] + 0.4 * X[,6] + 0.3 * X[,7]
  eta3 <-  0.2 + 0.6 * X[,2] - 0.6 * X[,4] + 0.7 * X[,8] - 0.5 * X[,9]
  eta4 <- -0.6 + 0.5 * X[,1] + 0.4 * X[,5] - 0.8 * X[,7] + 0.6 * X[,10]
  
  eta_fixed <- cbind(eta1, eta2, eta3, eta4)
  
  sim_cat5_core(
    X = X,
    id = id,
    eta_fixed = eta_fixed,
    intercept_sd = intercept_sd,
    return_prob = return_prob
  )
}


# ---------------------------------------------------------------
# Generate and save all four scenario datasets
#
# Each call uses seed = 123 for reproducibility. The seed is
# set both at the top of each scenario function (which controls
# predictor generation) and inside sim_cat5_core() (which
# controls random intercept draws and outcome sampling),
# so the full pipeline is deterministic.
#
# Output files are written to data/ relative to the working
# directory. The directory must exist before running this script;
# it is created by the folder structure of the research archive.
# ---------------------------------------------------------------
df1 <- sim_cat5_data1()
df2 <- sim_cat5_data2()
df3 <- sim_cat5_data3()
df4 <- sim_cat5_data4()

write.csv(df1, here("Multinominal_analysis", "data", "simulated_data1.csv"), row.names = FALSE)
write.csv(df2, here("Multinominal_analysis", "data", "simulated_data2.csv"), row.names = FALSE)
write.csv(df3, here("Multinominal_analysis", "data", "simulated_data3.csv"), row.names = FALSE)
write.csv(df4, here("Multinominal_analysis", "data", "simulated_data4.csv"), row.names = FALSE)
