# ===============================================================
# analysis.R
#
# 10-fold cross-validation analysis for the binary simulation
# study (Section 3.2 / Section 4.1.1 of the thesis).
#
# Six models are compared across four data-generating scenarios
# (linear, nonlinear, interactions, multicollinearity):
#
#   GMERT     Generalised Mixed-Effects Regression Tree,
#             binary case; fitted via mlml::fit_gmert_small()
#   GMERF     Generalised Mixed-Effects Random Forest,
#             binary case; fitted via mlml::fit_gmerf_small()
#   GLMM      Mixed-effects logistic regression (GLMM baseline);
#             lme4::glmer() with random intercept and random
#             slope on x1
#   Logit     Logistic regression (GLM baseline); stats::glm()
#   CART      Classification tree; rpart::rpart()
#   RF        Random forest; ranger::ranger()
#
# Evaluation metrics (Section 3.4):
#   Accuracy, majority-class F1, minority-class F1, relative bias
#
# Cross-validation design (Section 3.4):
#   K = 10 folds; observations within each cluster are
#   distributed across folds in a round-robin fashion so that
#   every fold contains data from every cluster. Folds are
#   constructed once per scenario and shared across all workers.
#
# Parallelisation:
#   One fold per worker (PSOCK cluster). Thread-level
#   parallelism inside each worker is disabled via environment
#   variables to prevent over-subscription.
#
# Inputs (expected in the calling environment):
#   df_list       Named list of data frames, one per scenario
#                 (s1, s2, s3, s4); produced by
#                 generate_binary_data.R.
#   scenario_cfg  Named list; each element holds model formulas
#                 and hyperparameter lists for one scenario:
#                   $glmm_formula, $logit_formula,
#                   $cart_formula, $rf_formula,
#                   $gmert_args,   $gmerf_args
#   K_folds       Integer; number of CV folds (10).
#   seed_folds    Integer; seed for fold construction.
#   seed_cluster  Integer; seed for the parallel RNG stream.
#
# Output:
#   results       Named list with one entry per scenario, each
#                 containing:
#                   $metrics_df   per-fold metric data frame
#                   $cv_summary   mean ± SD summary across folds
#                   $fold_table   per-fold convergence diagnostics
# ===============================================================


# ===============================================================
# SECTION 1: Helper functions
# ===============================================================

# ---------------------------------------------------------------
# make_cluster_folds()
#
# Constructs K cross-validation folds that respect the clustered
# structure of the data. Within each cluster, observations are
# shuffled and assigned to folds in round-robin order, ensuring
# that every fold contains observations from every cluster.
# This is required for mixed-effects models: a cluster absent
# from the training set would produce unestimable random effects.
#
# Arguments:
#   df    Data frame with an 'id' column identifying clusters.
#   K     Number of folds. Default 10.
#   seed  RNG seed for reproducibility. Default 42.
#
# Returns:
#   A list of K integer vectors, each containing the row
#   indices assigned to that fold (used as the test set).
# ---------------------------------------------------------------
make_cluster_folds <- function(df, K = 10, seed = 42) {
  set.seed(seed)
  folds <- vector("list", K)
  for (k in seq_len(K)) folds[[k]] <- integer(0)
  
  for (g in unique(df$id)) {
    idx_g      <- which(df$id == g)
    idx_g      <- sample(idx_g)    # shuffle within cluster
    split_g    <- split(idx_g, cut(seq_along(idx_g), K, labels = FALSE))
    for (k in seq_len(K)) folds[[k]] <- c(folds[[k]], split_g[[k]])
  }
  
  folds
}

# ---------------------------------------------------------------
# f1_fun()
#
# Computes the F1 score for either the majority or the minority
# class from a 2 × 2 confusion matrix (binary outcome).
#
# In the binary simulation the outcome is balanced by design
# (Table 2 / Appendix Table 7), but majority and minority F1
# are still reported separately because the real-data application
# is strongly imbalanced and the same function is reused there.
#
# The majority class is identified as the class with the higher
# marginal frequency in the confusion matrix (column sums).
# The minority class is the other one.
#
# Arguments:
#   cm        2 × 2 confusion matrix from table(predicted, actual).
#   majority  Logical; TRUE returns majority-class F1,
#             FALSE returns minority-class F1. Default TRUE.
#
# Returns: scalar F1 in [0, 1], or NA if precision or recall
#          is undefined (zero denominator).
# ---------------------------------------------------------------
f1_fun <- function(cm, majority = TRUE) {
  # Identify which class has more true observations
  class_counts  <- colSums(cm)
  majority_class <- names(which.max(class_counts))
  minority_class <- names(which.min(class_counts))
  
  target <- if (majority) majority_class else minority_class
  
  tp <- cm[target, target]
  fp <- sum(cm[target, ]) - tp
  fn <- sum(cm[, target]) - tp
  
  if ((tp + fp) == 0 || (tp + fn) == 0) return(NA_real_)
  
  prec <- tp / (tp + fp)
  rec  <- tp / (tp + fn)
  
  if ((prec + rec) == 0) return(NA_real_)
  2 * prec * rec / (prec + rec)
}

# ---------------------------------------------------------------
# bias_fun()
#
# Computes the relative prediction bias for the positive class
# (y = 1), defined in Section 3.4 as:
#
#   relative bias = (sum(y_hat) / sum(y_true)) - 1
#
# A value of 0 indicates the model predicts the positive class
# at exactly the correct frequency. Positive values indicate
# overprediction; negative values indicate underprediction.
#
# Arguments:
#   y_hat   Integer vector of predicted class labels (0/1).
#   y_true  Integer vector of true class labels (0/1).
#
# Returns: scalar; unbounded, centred at 0.
# ---------------------------------------------------------------
bias_fun <- function(y_hat, y_true) {
  sum(y_hat) / sum(y_true) - 1
}


# ===============================================================
# SECTION 2: Single-fold cross-validation worker
# ===============================================================

# ---------------------------------------------------------------
# cv_one_fold_bin()
#
# Fits all six models on the training portion of fold k and
# evaluates them on the held-out test portion. This function
# is the unit of work dispatched to each parallel worker.
#
# All six models receive exactly the same training and test
# split to ensure a fair comparison. Hard class predictions
# (0/1) are obtained by thresholding predicted probabilities
# at 0.5 for GLMM and Logit; GMERT, GMERF, CART, and RF
# return hard labels directly.
#
# GMERT and GMERF are fitted using mlml::fit_gmert_small() and
# mlml::fit_gmerf_small() respectively, which implement the
# binary PQL-EM algorithm described in Section 3.1.4 and
# Algorithm 1. Predictions for new observations use population-
# level fitted values (random effects set to zero), consistent
# with the cross-validation design where test observations come
# from partially held-out clusters.
#
# The fold result includes fitted model objects so that
# convergence diagnostics can be extracted in the summary step.
#
# Arguments:
#   k             Fold index (integer in 1:K).
#   df            Full data frame (all folds).
#   folds         List of K integer index vectors.
#   gmert_args    Named list of additional arguments for
#                 mlml::fit_gmert_small().
#   gmerf_args    Named list of additional arguments for
#                 mlml::fit_gmerf_small().
#   glmm_formula  Formula for lme4::glmer().
#   logit_formula Formula for stats::glm().
#   cart_formula  Formula for rpart::rpart().
#   rf_formula    Formula for ranger::ranger().
#   rf_seed       Integer seed for ranger. Default 42.
#
# Returns:
#   A list with elements:
#     fold_id    Integer fold index k.
#     train_idx  Integer vector of training row indices.
#     test_idx   Integer vector of test row indices.
#     gmert_fit  Fitted GMERT object (for convergence diagnostics).
#     gmerf_fit  Fitted GMERF object.
#     glmm_fit   Fitted glmer object.
#     logit_fit  Fitted glm object.
#     ctree_fit  Fitted rpart object.
#     rf_fit     Fitted ranger object.
#     metrics    Named list of scalar performance metrics.
#     confusion  Named list of 2 × 2 confusion matrices.
#     df_train   Training data frame for this fold.
#     df_test    Test data frame for this fold.
# ---------------------------------------------------------------
cv_one_fold_bin <- function(k, df, folds,
                            gmert_args    = list(),
                            gmerf_args    = list(),
                            glmm_formula,
                            logit_formula,
                            cart_formula,
                            rf_formula,
                            rf_seed = 42) {
  
  test_idx  <- folds[[k]]
  train_idx <- setdiff(seq_len(nrow(df)), test_idx)
  df_train  <- df[train_idx, ]
  df_test   <- df[test_idx, ]
  
  # -----------------------------------------------------------
  # Model 1: GMERT (binary)
  # Generalised Mixed-Effects Regression Tree for binary
  # outcomes, implemented in the mlml package (Section 3.1.4).
  # fit_gmert_small() is the binary-specific wrapper around
  # the PQL-EM algorithm (Algorithm 1).
  # -----------------------------------------------------------
  fit_gmert_k <- do.call(
    fit_gmert_small,
    c(list(df = df_train), gmert_args)
  )
  yhat_gmert <- predict_gmert(fit_gmert_k, df_test)
  
  # -----------------------------------------------------------
  # Model 2: GMERF (binary)
  # Generalised Mixed-Effects Random Forest for binary outcomes
  # (Section 3.1.4). Replaces the tree step of GMERT with a
  # random forest for a more stable fixed-effects estimate.
  # -----------------------------------------------------------
  fit_gmerf_k <- do.call(
    fit_gmerf_small,
    c(list(df = df_train), gmerf_args)
  )
  yhat_gmerf <- predict_gmerf(fit_gmerf_k, df_test)
  
  # -----------------------------------------------------------
  # Model 3: GLMM (mixed-effects logistic regression baseline)
  # lme4::glmer() with a random intercept and random slope on
  # x1, matching the data-generating random-effects structure
  # (Table 2). bobyqa optimizer improves convergence stability.
  # Predictions use type = "response" (fitted probabilities),
  # thresholded at 0.5 to produce hard class labels.
  # -----------------------------------------------------------
  fit_glmm_k <- glmer(
    glmm_formula,
    data    = df_train,
    family  = binomial(link = "logit"),
    control = glmerControl(optimizer = "bobyqa")
  )
  p_glmm    <- as.numeric(predict(fit_glmm_k, newdata = df_test, type = "response"))
  yhat_glmm <- as.integer(p_glmm >= 0.5)
  
  # -----------------------------------------------------------
  # Model 4: Logit (logistic regression GLM baseline)
  # Standard logistic regression with no random effects;
  # treats all observations as independent (Section 3.1.1).
  # Serves as the parametric single-level baseline.
  # -----------------------------------------------------------
  fit_logit_k <- glm(logit_formula, data = df_train, family = binomial())
  p_logit     <- as.numeric(predict(fit_logit_k, newdata = df_test, type = "response"))
  yhat_logit  <- as.integer(p_logit >= 0.5)
  
  # -----------------------------------------------------------
  # Model 5: CART (single-level tree baseline)
  # rpart::rpart() with no random effects (Section 3.1.3).
  # y must be a factor for method = "class"; it is converted
  # locally to avoid modifying df_train for subsequent models.
  #
  # Hyperparameters:
  #   cp = 0.0        no cost-complexity pruning
  #   minsplit = 20   minimum node size before attempting split
  #   minbucket = 7   minimum terminal node size
  #   maxdepth = 5    maximum tree depth
  # -----------------------------------------------------------
  df_train_ct    <- df_train
  df_train_ct$y  <- factor(df_train_ct$y, levels = c(0, 1))
  df_test_ct     <- df_test
  df_test_ct$y   <- factor(df_test_ct$y,  levels = c(0, 1))
  
  fit_cart_k <- rpart::rpart(
    cart_formula,
    data    = df_train_ct,
    method  = "class",
    control = rpart::rpart.control(
      cp        = 0.0,
      minsplit  = 20,
      minbucket = 7,
      maxdepth  = 5
    )
  )
  pred_cart_class <- predict(fit_cart_k, newdata = df_test_ct, type = "class")
  yhat_cart       <- as.integer(as.character(pred_cart_class))
  
  # -----------------------------------------------------------
  # Model 6: RF (random forest single-level ensemble baseline)
  # ranger::ranger() with no random effects (Section 3.1.3).
  # classification = TRUE and probability = FALSE return hard
  # integer class predictions directly.
  #
  # Hyperparameters:
  #   num.trees = 500   standard ensemble size
  #   mtry = 2          ~sqrt(p) for p = 3 or 4 predictors
  #   min.node.size = 5 minimum terminal node size
  # -----------------------------------------------------------
  fit_rf_k <- ranger::ranger(
    formula       = rf_formula,
    data          = df_train,
    num.trees     = 500,
    mtry          = 2,
    min.node.size = 5,
    classification = TRUE,
    probability   = FALSE,
    seed          = rf_seed
  )
  yhat_rf <- as.integer(predict(fit_rf_k, data = df_test)$predictions)
  
  # -----------------------------------------------------------
  # Confusion matrices (one per model)
  # Stored in the return value for optional downstream use
  # (e.g. per-class breakdowns, convergence tables).
  # -----------------------------------------------------------
  cm_gmert  <- table(predicted = yhat_gmert, actual = df_test$y)
  cm_gmerf  <- table(predicted = yhat_gmerf, actual = df_test$y)
  cm_glmm   <- table(predicted = yhat_glmm,  actual = df_test$y)
  cm_logit  <- table(predicted = yhat_logit, actual = df_test$y)
  cm_cart   <- table(predicted = yhat_cart,  actual = df_test$y)
  cm_rf     <- table(predicted = yhat_rf,    actual = df_test$y)
  
  # -----------------------------------------------------------
  # Performance metrics for this fold
  # Accuracy, majority F1, minority F1, and relative bias
  # are computed for all six models (Section 3.4).
  # -----------------------------------------------------------
  metrics <- list(
    acc_gmert  = mean(yhat_gmert == df_test$y),
    acc_gmerf  = mean(yhat_gmerf == df_test$y),
    acc_glmm   = mean(yhat_glmm  == df_test$y),
    acc_logit  = mean(yhat_logit == df_test$y),
    acc_cart   = mean(yhat_cart  == df_test$y),
    acc_rf     = mean(yhat_rf    == df_test$y),
    
    f1_mag_gmert  = f1_fun(cm_gmert, majority = TRUE),
    f1_mag_gmerf  = f1_fun(cm_gmerf, majority = TRUE),
    f1_mag_glmm   = f1_fun(cm_glmm,  majority = TRUE),
    f1_mag_logit  = f1_fun(cm_logit, majority = TRUE),
    f1_mag_cart   = f1_fun(cm_cart,  majority = TRUE),
    f1_mag_rf     = f1_fun(cm_rf,    majority = TRUE),
    
    f1_min_gmert  = f1_fun(cm_gmert, majority = FALSE),
    f1_min_gmerf  = f1_fun(cm_gmerf, majority = FALSE),
    f1_min_glmm   = f1_fun(cm_glmm,  majority = FALSE),
    f1_min_logit  = f1_fun(cm_logit, majority = FALSE),
    f1_min_cart   = f1_fun(cm_cart,  majority = FALSE),
    f1_min_rf     = f1_fun(cm_rf,    majority = FALSE),
    
    bias_gmert  = bias_fun(yhat_gmert, df_test$y),
    bias_gmerf  = bias_fun(yhat_gmerf, df_test$y),
    bias_glmm   = bias_fun(yhat_glmm,  df_test$y),
    bias_logit  = bias_fun(yhat_logit, df_test$y),
    bias_cart   = bias_fun(yhat_cart,  df_test$y),
    bias_rf     = bias_fun(yhat_rf,    df_test$y)
  )
  
  list(
    fold_id   = k,
    train_idx = train_idx,
    test_idx  = test_idx,
    gmert_fit = fit_gmert_k,
    gmerf_fit = fit_gmerf_k,
    glmm_fit  = fit_glmm_k,
    logit_fit = fit_logit_k,
    ctree_fit = fit_cart_k,
    rf_fit    = fit_rf_k,
    metrics   = metrics,
    confusion = list(
      gmert = cm_gmert,
      gmerf = cm_gmerf,
      glmm  = cm_glmm,
      logit = cm_logit,
      cart  = cm_cart,
      rf    = cm_rf
    ),
    df_train = df_train,
    df_test  = df_test
  )
}


# ===============================================================
# SECTION 3: Cross-validation summary
# ===============================================================

# ---------------------------------------------------------------
# summarise_cv_bin()
#
# Aggregates per-fold metrics returned by cv_one_fold_bin()
# into a tidy summary (mean ± SD across folds), matching the
# reporting format used in the thesis (Section 4.1.1 and
# Appendix Table 9).
#
# Also builds a per-fold convergence diagnostic table for
# GMERT, reporting whether the outer PQL loop converged and
# the percentage of inner EM iterations that converged.
#
# Arguments:
#   cv_results  List of K outputs from cv_one_fold_bin().
#
# Returns:
#   A list with three elements:
#     $metrics_df   K-row data frame, one column per
#                   model-metric combination.
#     $cv_summary   6-row data frame with mean and SD of
#                   each metric across folds.
#     $fold_table   K-row convergence diagnostic table.
# ---------------------------------------------------------------
summarise_cv_bin <- function(cv_results) {
  
  # Stack per-fold metrics into a single data frame
  metrics_df <- do.call(rbind, lapply(cv_results, function(x) {
    data.frame(
      fold         = x$fold_id,
      acc_gmert    = x$metrics$acc_gmert,
      acc_gmerf    = x$metrics$acc_gmerf,
      acc_glmm     = x$metrics$acc_glmm,
      acc_logit    = x$metrics$acc_logit,
      acc_cart     = x$metrics$acc_cart,
      acc_rf       = x$metrics$acc_rf,
      f1_mag_gmert = x$metrics$f1_mag_gmert,
      f1_mag_gmerf = x$metrics$f1_mag_gmerf,
      f1_mag_glmm  = x$metrics$f1_mag_glmm,
      f1_mag_logit = x$metrics$f1_mag_logit,
      f1_mag_cart  = x$metrics$f1_mag_cart,
      f1_mag_rf    = x$metrics$f1_mag_rf,
      f1_min_gmert = x$metrics$f1_min_gmert,
      f1_min_gmerf = x$metrics$f1_min_gmerf,
      f1_min_glmm  = x$metrics$f1_min_glmm,
      f1_min_logit = x$metrics$f1_min_logit,
      f1_min_cart  = x$metrics$f1_min_cart,
      f1_min_rf    = x$metrics$f1_min_rf,
      bias_gmert   = x$metrics$bias_gmert,
      bias_gmerf   = x$metrics$bias_gmerf,
      bias_glmm    = x$metrics$bias_glmm,
      bias_logit   = x$metrics$bias_logit,
      bias_cart    = x$metrics$bias_cart,
      bias_rf      = x$metrics$bias_rf
    )
  }))
  
  # Compute mean and SD across folds for each model and metric
  cv_summary <- data.frame(
    model = c("GMERT", "GMERF", "GLMM", "Logit", "CART", "RF"),
    
    acc_mean = c(mean(metrics_df$acc_gmert),  mean(metrics_df$acc_gmerf),
                 mean(metrics_df$acc_glmm),   mean(metrics_df$acc_logit),
                 mean(metrics_df$acc_cart),   mean(metrics_df$acc_rf)),
    acc_sd   = c(sd(metrics_df$acc_gmert),    sd(metrics_df$acc_gmerf),
                 sd(metrics_df$acc_glmm),     sd(metrics_df$acc_logit),
                 sd(metrics_df$acc_cart),     sd(metrics_df$acc_rf)),
    
    f1_mag_mean = c(mean(metrics_df$f1_mag_gmert), mean(metrics_df$f1_mag_gmerf),
                    mean(metrics_df$f1_mag_glmm),  mean(metrics_df$f1_mag_logit),
                    mean(metrics_df$f1_mag_cart),  mean(metrics_df$f1_mag_rf)),
    f1_mag_sd   = c(sd(metrics_df$f1_mag_gmert),   sd(metrics_df$f1_mag_gmerf),
                    sd(metrics_df$f1_mag_glmm),    sd(metrics_df$f1_mag_logit),
                    sd(metrics_df$f1_mag_cart),    sd(metrics_df$f1_mag_rf)),
    
    f1_min_mean = c(mean(metrics_df$f1_min_gmert), mean(metrics_df$f1_min_gmerf),
                    mean(metrics_df$f1_min_glmm),  mean(metrics_df$f1_min_logit),
                    mean(metrics_df$f1_min_cart),  mean(metrics_df$f1_min_rf)),
    f1_min_sd   = c(sd(metrics_df$f1_min_gmert),   sd(metrics_df$f1_min_gmerf),
                    sd(metrics_df$f1_min_glmm),    sd(metrics_df$f1_min_logit),
                    sd(metrics_df$f1_min_cart),    sd(metrics_df$f1_min_rf)),
    
    bias_mean = c(mean(metrics_df$bias_gmert), mean(metrics_df$bias_gmerf),
                  mean(metrics_df$bias_glmm),  mean(metrics_df$bias_logit),
                  mean(metrics_df$bias_cart),  mean(metrics_df$bias_rf)),
    bias_sd   = c(sd(metrics_df$bias_gmert),   sd(metrics_df$bias_gmerf),
                  sd(metrics_df$bias_glmm),    sd(metrics_df$bias_logit),
                  sd(metrics_df$bias_cart),    sd(metrics_df$bias_rf))
  )
  
  # Per-fold GMERT convergence diagnostics.
  # converged_in stores a logical vector of inner-loop
  # convergence flags; its mean gives the proportion of inner
  # iterations that converged. converged_out is a single logical
  # indicating whether the outer PQL loop converged.
  fold_table <- do.call(rbind, lapply(cv_results, function(x) {
    inner_vec <- x$gmert_fit$converged_in
    inner_pct <- if (length(inner_vec)) mean(inner_vec) * 100 else NA_real_
    
    data.frame(
      Fold               = x$fold_id,
      Outer_Converged    = isTRUE(x$gmert_fit$converged_out),
      Inner_Conv_Percent = round(inner_pct, 1),
      F1_GMERT           = round(x$metrics$f1_mag_gmert, 4),
      F1_GMERF           = round(x$metrics$f1_mag_gmerf, 4),
      F1_GLMM            = round(x$metrics$f1_mag_glmm,  4),
      F1_Logit           = round(x$metrics$f1_mag_logit, 4),
      F1_CART            = round(x$metrics$f1_mag_cart,  4),
      F1_RF              = round(x$metrics$f1_mag_rf,    4)
    )
  }))
  
  fold_table <- fold_table[order(fold_table$Fold), ]
  
  list(
    metrics_df = metrics_df,
    cv_summary = cv_summary,
    fold_table = fold_table
  )
}


# ===============================================================
# SECTION 4: Parallel cluster setup
# ===============================================================

# One worker per fold, leaving one core free for the main process
n_workers <- max(1L, min(K_folds, parallel::detectCores() - 1L))
cl <- parallel::makeCluster(n_workers, type = "PSOCK")

# Load required packages on every worker.
# mlml provides fit_gmert_small(), fit_gmerf_small(), and their
# predict() methods; lme4 provides glmer().
parallel::clusterEvalQ(cl, {
  library(mlml)
  library(lme4)
  library(rpart)
  library(ranger)
  library(MASS)
  
  # Disable thread-level parallelism inside each worker to
  # prevent over-subscription across the PSOCK processes
  Sys.setenv(
    OMP_NUM_THREADS      = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS      = "1",
    OMP_THREAD_LIMIT     = "1"
  )
})

# Export all functions needed inside cv_one_fold_bin() to
# each worker; PSOCK workers start with a clean environment
parallel::clusterExport(
  cl,
  varlist = c(
    "cv_one_fold_bin",
    "make_cluster_folds",
    "f1_fun",
    "bias_fun",
    "summarise_cv_bin",
    # mlml internals required by fit_gmert_small / fit_gmerf_small
    "fit_gmert_small",
    "predict_gmert",
    "fit_gmerf_small",
    "predict_gmerf",
    "Ajnv_fun",
    "b_fun_small",
    "D_fun_small",
    "sigma_fun_small",
    "gll_fun"
  ),
  envir = environment()
)


# ===============================================================
# SECTION 5: Main CV loop — one scenario at a time
# ===============================================================

results <- list()

for (nm in names(df_list)) {
  
  df  <- df_list[[nm]]
  cfg <- scenario_cfg[[nm]]
  
  # Build folds once and share with all workers so every model
  # sees identical splits
  folds <- make_cluster_folds(df, K = K_folds, seed = seed_folds)
  
  # Re-export scenario-specific data and folds at each iteration
  parallel::clusterExport(cl, varlist = c("df", "folds"), envir = environment())
  
  # Reproducible parallel RNG: L'Ecuyer-CMRG streams ensure
  # each worker gets a statistically independent RNG sequence
  parallel::clusterSetRNGStream(cl, iseed = seed_cluster)
  
  # Dispatch one fold per worker; parLapply() blocks until all
  # K folds complete, then returns a list of K result bundles
  cv_res <- parallel::parLapply(
    cl,
    X   = seq_len(K_folds),
    fun = function(k, glmm_formula, logit_formula, cart_formula,
                   rf_formula, gmert_args, gmerf_args) {
      cv_one_fold_bin(
        k             = k,
        df            = df,
        folds         = folds,
        gmert_args    = gmert_args,
        gmerf_args    = gmerf_args,
        glmm_formula  = glmm_formula,
        logit_formula = logit_formula,
        cart_formula  = cart_formula,
        rf_formula    = rf_formula,
        rf_seed       = 42
      )
    },
    glmm_formula  = cfg$glmm_formula,
    logit_formula = cfg$logit_formula,
    cart_formula  = cfg$cart_formula,
    rf_formula    = cfg$rf_formula,
    gmert_args    = cfg$gmert_args,
    gmerf_args    = cfg$gmerf_args
  )
  
  # Aggregate per-fold results into mean ± SD summary
  summ <- summarise_cv_bin(cv_res)
  
  results[[nm]] <- list(
    metrics_df = summ$metrics_df,
    cv_summary = summ$cv_summary,
    fold_table = summ$fold_table
  )
  
  cat("Done:", nm, "\n")
}

parallel::stopCluster(cl)