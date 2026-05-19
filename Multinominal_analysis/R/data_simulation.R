# ===============================================================
# Simulate clustered categorical data with 5 classes
# Random intercept only, 10 predictors
# ===============================================================

# ---------------------------------------------------------------
# Common core: takes X and a fixed-part eta matrix with 4 logits
# (class 5 is the reference category)
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
  K <- 5
  K1 <- 4

  groups <- as.integer(as.factor(id))
  G <- length(unique(groups))

  colnames(X) <- paste0("x", seq_len(p))

  # Random intercept covariance across the 4 logits
  if (is.null(D_intercept)) {
    D_intercept <- diag(intercept_sd^2, K1)
  } else {
    D_intercept <- as.matrix(D_intercept)
    if (!all(dim(D_intercept) == c(K1, K1))) {
      stop("D_intercept must be a 4 x 4 matrix.")
    }
  }

  # Draw group random intercepts: G x 4
  b0 <- MASS::mvrnorm(n = G, mu = rep(0, K1), Sigma = D_intercept)

  # Add random intercepts
  eta <- eta_fixed
  for (g in seq_len(G)) {
    idx <- which(groups == g)
    eta[idx, ] <- sweep(eta[idx, , drop = FALSE], 2, b0[g, ], FUN = "+")
  }

  # Softmax with class 5 as reference
  exp_eta <- exp(eta)
  denom <- 1 + rowSums(exp_eta)

  P <- matrix(0, N, K)
  P[, 1:K1] <- exp_eta / denom
  P[, K] <- 1 / denom

  # Draw categorical outcome
  y_int <- apply(P, 1, function(pv) sample.int(K, size = 1L, prob = pv))
  y <- factor(y_int, levels = 1:5, labels = paste0("class", 1:5))

  df <- data.frame(
    id = factor(groups),
    y = y,
    X,
    stringsAsFactors = FALSE
  )

  if (return_prob) {
    colnames(P) <- paste0("prob_", 1:5)
    df <- cbind(df, as.data.frame(P))
  }

  df
}


# ---------------------------------------------------------------
# Scenario 1: linear / additive
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

  N <- G * n_i
  p <- 10
  id <- rep(seq_len(G), each = n_i)

  Sigma_X <- matrix(rho, p, p)
  diag(Sigma_X) <- 1
  X <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = Sigma_X)

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
# Scenario 2: nonlinearity
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
# Scenario 3: interactions
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
# Scenario 4: multicollinearity
# ---------------------------------------------------------------
sim_cat5_data4 <- function(
    G = 100,
    n_i = 18,
    rho = 0.85,
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

df1 <- sim_cat5_data1()
df2 <- sim_cat5_data2()
df3 <- sim_cat5_data3()
df4 <- sim_cat5_data4()

write.csv(df1, "data/simulated_data1.csv", row.names = FALSE)
write.csv(df2, "data/simulated_data2.csv", row.names = FALSE)
write.csv(df3, "data/simulated_data3.csv", row.names = FALSE)
write.csv(df4, "data/simulated_data4.csv", row.names = FALSE)

table(df1$y)
table(df2$y)
table(df3$y)
table(df4$y)

table(df1$id)
table(df2$id)
table(df3$id)
table(df4$id)

# Binary real-data results for car
df_car_binary <- data.frame(
  Model = c("GMERT", "GMERF", "CART", "RF"),
  acc_mean = c(0.846, 0.903, 0.889, 0.932),
  acc_sd = c(0.141, 0.012, 0.011, 0.004),
  f1_maj_mean = c(0.871, 0.909, 0.897, 0.936),
  f1_maj_sd = c(0.090, 0.012, 0.011, 0.004),
  f1_min_mean = c(0.793, 0.895, 0.880, 0.927),
  f1_min_sd = c(0.277, 0.012, 0.014, 0.005),
  bias_maj_mean = c(-0.099, -0.014, 0.003, -0.026),
  bias_maj_sd = c(0.307, 0.035, 0.059, 0.015),
  bias_min_mean = c(-0.055, 0.017, -0.003, 0.031),
  bias_min_sd = c(0.169, 0.041, 0.070, 0.018)
)

# Binary real-data results for public transport
df_public_transport_binary <- data.frame(
  Model = c("GMERT", "GMERF", "CART", "RF"),
  acc_mean = c(0.9277727, 0.8950735, 0.9779384, 0.9782031),
  acc_sd = c(0.058470372, 0.136633959, 0.001268718, 0.001410695),
  f1_maj_mean = c(0.9612439, 0.8793850, 0.9888458, 0.9889769),
  f1_maj_sd = c(0.0325691062, 0.2787843865, 0.0006483672, 0.0007178679),
  f1_min_mean = c(0.10853527, 0.16682834, 0.00000000, 0.03279757),
  f1_min_sd = c(0.10981918, 0.18723219, 0.00000000, 0.05266226),
  bias_maj_mean = c(1.7852679, 3.2233698, -1.0000000, -0.9765451),
  bias_maj_sd = c(3.03852142, 6.35746774, 0.00000000, 0.04681103),
  bias_min_mean = c(0.04814595, -0.07536780, 0.02256086, 0.02204511),
  bias_min_sd = c(0.079419514, 0.151795309, 0.001325692, 0.001851763)
)

# Multinomial real-data results
df_multinomial <- data.frame(
  Model = c("GMERT-cat", "GMERF-cat", "CART", "RF"),
  acc_mean = c(0.665, 0.459, 0.863, 0.904),
  acc_sd = c(0.290, 0.296, 0.013, 0.009),
  macro_f1_mean = c(0.584, 0.399, 0.840, 0.859),
  macro_f1_sd = c(0.218, 0.195, 0.029, 0.026),
  weighted_f1_mean = c(0.645, 0.401, 0.852, 0.893),
  weighted_f1_sd = c(0.308, 0.331, 0.014, 0.011),
  br_mean = c(1.476, 1.377, 0.222, 0.191),
  br_sd = c(0.222, 0.225, 0.027, 0.010),
  ce_mean = c(25.147, 25.050, 0.554, 0.417),
  ce_sd = c(4.717, 4.877, 0.069, 0.019)
)