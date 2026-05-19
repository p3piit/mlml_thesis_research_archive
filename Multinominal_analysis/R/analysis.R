# ===============================================================
# analysis_categorical.R
#
# 10-fold cross-validation analysis for the multinomial
# simulation study (Section 3.2 / Section 4.1.2 of the thesis).
#
# Six models are compared across four data-generating scenarios
# (linear, nonlinear, interactions, multicollinearity):
#
#   GMERT-cat   Generalised Mixed-Effects Regression Tree,
#               multinomial extension (Appendix B); fitted via
#               mlml::fit_gmert_cat()
#   GMERF-cat   Generalised Mixed-Effects Random Forest,
#               multinomial extension (Appendix B); fitted via
#               mlml::fit_gmerf_cat()
#   Multinom    Multinomial logit (GLM baseline); nnet::multinom()
#   MBLOGIT     Mixed-effects multinomial logit (GLMM baseline);
#               mclogit::mblogit() with a cluster random intercept
#   CART        Classification tree; rpart::rpart()
#   RF          Random forest; ranger::ranger()
#
# Evaluation metrics (Section 3.4):
#   Accuracy, macro F1, weighted F1, Brier score, cross-entropy
#
# Cross-validation design (Section 3.4):
#   K = 10 folds; observations within each cluster are
#   distributed across folds so that every fold contains
#   data from every cluster (split_by_cluster = FALSE).
#   Folds are constructed before the parallel loop and shared
#   with all workers so every model sees the exact same splits.
#
# Parallelisation:
#   One fold per worker (PSOCK cluster). Thread-level
#   parallelism inside each worker is disabled via environment
#   variables to prevent over-subscription.
#
# Inputs (expected in the calling environment):
#   df_list       Named list of data frames, one per scenario.
#                 Names must match keys in scenario_cfg.
#   scenario_cfg  Named list; each element holds model formulas
#                 and hyperparameter lists for one scenario:
#                   $multinom_formula, $mblogit_formula,
#                   $cart_formula, $rf_formula,
#                   $gmert_args, $gmerf_args
#   K_folds       Integer; number of CV folds (typically 10).
#   seed_folds    Integer; seed for fold construction.
#   seed_cluster  Integer; seed for the parallel RNG stream.
#
# Output:
#   results       Named list with one entry per scenario, each
#                 containing:
#                   $metrics_df   per-fold metric data frame
#                   $cv_summary   mean ± SD summary across folds
# ===============================================================


# ===============================================================
# SECTION 1: Helper functions
# ===============================================================

# ---------------------------------------------------------------
# make_cluster_folds()
#
# Constructs K stratified cross-validation folds that respect
# the clustered structure of the data.
#
# Two allocation strategies are available via split_by_cluster:
#
#   FALSE (default) — observation-level split:
#     Within each cluster, observations are shuffled and then
#     assigned to folds in a round-robin fashion. This ensures
#     that every fold contains observations from every cluster,
#     which is important for fitting mixed-effects models: a
#     cluster absent from the training set would produce
#     unestimable random effects.
#
#   TRUE — cluster-level split:
#     Entire clusters are assigned to folds. Clusters seen
#     during testing are not seen during training, which tests
#     out-of-sample generalisation to new clusters. Not used
#     in the simulation study but provided for completeness.
#
# Arguments:
#   df              Data frame with an 'id' column identifying
#                   cluster membership.
#   K               Number of folds. Default 10.
#   seed            RNG seed for reproducibility. Default 42.
#   split_by_cluster Logical; see above. Default FALSE.
#
# Returns:
#   A list of K integer vectors, each containing the row
#   indices assigned to that fold (used as the test set).
# ---------------------------------------------------------------
make_cluster_folds <- function(df, K = 10, seed = 42, split_by_cluster = FALSE) {
  set.seed(seed)
  folds <- vector("list", K)
  for (k in seq_len(K)) folds[[k]] <- integer(0)
  ids <- unique(as.character(df$id))
  
  if (split_by_cluster) {
    # Cluster-level split: shuffle cluster ids and assign each
    # cluster to one fold in round-robin order
    shuffled_ids <- sample(ids)
    assignments <- rep(seq_len(K), length.out = length(shuffled_ids))
    for (i in seq_along(shuffled_ids)) {
      k <- assignments[i]
      folds[[k]] <- c(folds[[k]], which(as.character(df$id) == shuffled_ids[i]))
    }
  } else {
    # Observation-level split: within each cluster, shuffle
    # observations and distribute them across folds round-robin.
    # This keeps every cluster represented in every fold.
    for (g in ids) {
      idx_g <- which(as.character(df$id) == g)
      if (length(idx_g) == 0) next
      shuffled_idx_g <- sample(idx_g)
      assignments <- rep(seq_len(K), length.out = length(idx_g))
      for (i in seq_along(shuffled_idx_g)) {
        k <- assignments[i]
        folds[[k]] <- c(folds[[k]], shuffled_idx_g[i])
      }
    }
  }
  
  folds
}

# ---------------------------------------------------------------
# acc_fun_mc()
#
# Overall classification accuracy for a multiclass outcome.
# Defined as the proportion of observations for which the
# predicted class matches the true class (Section 3.4).
#
# Arguments:
#   y_true  True class labels (factor or character).
#   y_pred  Predicted class labels (same type as y_true).
#
# Returns: scalar in [0, 1].
# ---------------------------------------------------------------
acc_fun_mc <- function(y_true, y_pred) {
  mean(as.character(y_true) == as.character(y_pred))
}

# ---------------------------------------------------------------
# f1_per_class()
#
# Computes the F1 score for each class separately, treating
# each class in turn as the "positive" class in a one-vs-rest
# binary comparison. Returns NA for any class with no true
# positives and no predicted positives (precision undefined)
# or no true positives and no false negatives (recall
# undefined), as defined in Section 3.4.
#
# Arguments:
#   y_true   True class labels (factor).
#   y_pred   Predicted class labels (factor or character).
#   classes  Character vector of class names. Defaults to
#            levels(y_true).
#
# Returns: named numeric vector of length length(classes).
# ---------------------------------------------------------------
f1_per_class <- function(y_true, y_pred, classes = levels(y_true)) {
  sapply(classes, function(cl) {
    tp <- sum(y_true == cl & y_pred == cl)
    fp <- sum(y_true != cl & y_pred == cl)
    fn <- sum(y_true == cl & y_pred != cl)
    
    # Return NA if precision or recall is undefined
    if ((tp + fp) == 0 || (tp + fn) == 0) return(NA_real_)
    
    prec <- tp / (tp + fp)
    rec  <- tp / (tp + fn)
    
    if ((prec + rec) == 0) return(NA_real_)
    2 * prec * rec / (prec + rec)
  })
}

# ---------------------------------------------------------------
# macro_f1_fun()
#
# Macro-averaged F1: unweighted mean of per-class F1 scores.
# Gives equal weight to every class regardless of its frequency,
# making it sensitive to performance on rare classes.
# NA values (classes with undefined F1) are excluded from the
# mean (Section 3.4).
#
# Arguments:
#   y_true, y_pred, classes  — as in f1_per_class().
#
# Returns: scalar in [0, 1].
# ---------------------------------------------------------------
macro_f1_fun <- function(y_true, y_pred, classes = levels(y_true)) {
  mean(f1_per_class(y_true, y_pred, classes = classes), na.rm = TRUE)
}

# ---------------------------------------------------------------
# weighted_f1_fun()
#
# Weighted-averaged F1: mean of per-class F1 scores weighted by
# the observed class frequencies in y_true. Gives more weight
# to frequent classes, making it closer to overall accuracy
# for imbalanced outcomes (Section 3.4).
# NA F1 values contribute 0 to the numerator (na.rm behaviour
# of sum) but their weight is still excluded from the denominator
# via the sum(wts) term — only non-NA classes contribute.
#
# Arguments:
#   y_true, y_pred, classes  — as in f1_per_class().
#
# Returns: scalar in [0, 1].
# ---------------------------------------------------------------
weighted_f1_fun <- function(y_true, y_pred, classes = levels(y_true)) {
  f1s <- f1_per_class(y_true, y_pred, classes = classes)
  wts <- as.numeric(table(factor(y_true, levels = classes)))
  sum(f1s * wts, na.rm = TRUE) / sum(wts)
}

# ---------------------------------------------------------------
# brier_score()
#
# Multiclass Brier score: mean squared error between the
# predicted probability matrix and the one-hot encoded true
# labels. Lower is better; the score is 0 for perfect
# probability predictions and at most 2 for K classes.
#
# Note: df_test is referenced directly from the calling scope
# rather than through y_true — this is a known quirk of the
# current implementation and should be unified in a future
# revision.
#
# Arguments:
#   prob_mat  N x K matrix of predicted class probabilities.
#   y_true    Length-N factor of true class labels.
#   classes   Character vector of K class names.
#
# Returns: scalar >= 0.
# ---------------------------------------------------------------
brier_score <- function(prob_mat, y_true, classes) {
  Y_oh <- model.matrix(~ 0 + factor(df_test$y, levels = classes))
  mean(rowSums((prob_mat - Y_oh)^2))
}

# ---------------------------------------------------------------
# cross_entropy()
#
# Mean cross-entropy (log loss) between predicted probabilities
# and the true class labels. Equivalent to the negative mean
# log-likelihood of the predicted distribution. Lower is better.
# Predicted probabilities are clipped at 1e-15 to avoid log(0).
#
# Arguments:
#   prob_mat  N x K matrix of predicted class probabilities.
#   y_true    Length-N factor or character of true class labels.
#   classes   Character vector of K class names used to index
#             the columns of prob_mat.
#
# Returns: scalar >= 0.
# ---------------------------------------------------------------
cross_entropy <- function(prob_mat, y_true, classes) {
  idx <- match(as.character(y_true), classes)
  probs <- prob_mat[cbind(seq_along(idx), idx)]
  -mean(log(pmax(probs, 1e-15)))
}


# ===============================================================
# SECTION 2: Single-fold cross-validation worker
# ===============================================================

# ---------------------------------------------------------------
# cv_one_fold_cat()
#
# Fits all six models on the training portion of fold k and
# evaluates them on the held-out test portion. This function
# is the unit of work dispatched to each parallel worker.
#
# All six models receive exactly the same training and test
# split (determined by folds[[k]]) to ensure a fair comparison.
#
# For GMERT-cat and GMERF-cat, predictions are made at the
# population level (random_effect = NULL), meaning the random
# effects are set to zero at test time. This corresponds to
# predicting for observations from clusters not seen during
# training, consistent with the cross-validation design where
# some test observations come from partially held-out clusters.
# The same marginal (population-level) prediction strategy is
# used for MBLOGIT (conditional = FALSE).
#
# For models that may not predict all K classes in a given fold
# (e.g. if a rare class is absent from the training set), the
# probability matrix is padded with zeros for the missing
# classes so that all downstream metric functions receive a
# consistently shaped N x K matrix.
#
# Arguments:
#   k                Fold index (integer in 1:K).
#   df               Full data frame (all folds).
#   folds            List of K integer index vectors from
#                    make_cluster_folds().
#   gmert_args       Named list of additional arguments passed
#                    to mlml::fit_gmert_cat().
#   gmerf_args       Named list of additional arguments passed
#                    to mlml::fit_gmerf_cat().
#   multinom_formula Formula for nnet::multinom().
#   mblogit_formula  Formula for mclogit::mblogit().
#   cart_formula     Formula for rpart::rpart().
#   rf_formula       Formula for ranger::ranger().
#   rf_seed          Integer seed for ranger's internal RNG.
#                    Default 42.
#
# Returns:
#   A list with elements:
#     fold_id    Integer fold index k.
#     train_idx  Integer vector of training row indices.
#     test_idx   Integer vector of test row indices.
#     metrics    Named list of scalar performance metrics.
#     df_train   Training data frame for this fold.
#     df_test    Test data frame for this fold.
# ---------------------------------------------------------------
cv_one_fold_cat <- function(k, df, folds,
                            gmert_args = list(),
                            gmerf_args = list(),
                            multinom_formula,
                            mblogit_formula,
                            cart_formula,
                            rf_formula,
                            rf_seed = 42) {
  
  # Split data into training and test sets for this fold
  test_idx  <- folds[[k]]
  train_idx <- setdiff(seq_len(nrow(df)), test_idx)
  
  df_train <- df[train_idx, ]
  df_test  <- df[test_idx, ]
  
  # Class labels used for padding and metric computation
  classes <- levels(df$y)
  
  # -----------------------------------------------------------
  # Model 1: GMERT-cat
  # Generalised Mixed-Effects Regression Tree for multinomial
  # outcomes, implemented in the mlml package (Appendix B).
  # random_effects = NULL uses a random intercept only,
  # matching the data-generating structure (Table 3).
  # predict() is called twice: once for class labels and once
  # for the full probability matrix (prob_saved = TRUE).
  # -----------------------------------------------------------
  fit_gmert_k <- do.call(
    fit_gmert_cat,
    c(list(
      df = df_train,
      id = "id",
      target = "y",
      random_effects = NULL
    ), gmert_args)
  )
  
  pred_gmert <- predict_gmert_cat(
    fit = fit_gmert_k,
    new_df = df_test,
    random_effect = NULL,   # population-level prediction
    id = "id"
  )
  
  prob_gmert <- predict_gmert_cat(
    fit = fit_gmert_k,
    new_df = df_test,
    random_effect = NULL,
    id = "id",
    prob_saved = TRUE       # return full N x K probability matrix
  )
  
  # -----------------------------------------------------------
  # Model 2: GMERF-cat
  # Generalised Mixed-Effects Random Forest for multinomial
  # outcomes (Appendix B). Same interface as GMERT-cat but
  # uses a random forest for the fixed-effects component,
  # which reduces variance at the cost of interpretability.
  # -----------------------------------------------------------
  fit_gmerf_k <- do.call(
    fit_gmerf_cat,
    c(list(
      df = df_train,
      id = "id",
      target = "y",
      random_effects = NULL
    ), gmerf_args)
  )
  
  pred_gmerf <- predict_gmerf_cat(
    fit = fit_gmerf_k,
    new_df = df_test,
    random_effect = NULL,
    id = "id"
  )
  
  prob_gmerf <- predict_gmerf_cat(
    fit = fit_gmerf_k,
    new_df = df_test,
    random_effect = NULL,
    id = "id",
    prob_saved = TRUE
  )
  
  # -----------------------------------------------------------
  # Model 3: Multinomial logit (GLM baseline)
  # Standard multinomial logistic regression via nnet::multinom.
  # No random effects; treats all observations as independent.
  # This is the parametric single-level baseline (Section 3.1.1).
  # trace = FALSE suppresses iteration output.
  # -----------------------------------------------------------
  fit_multinom_k <- nnet::multinom(
    formula = multinom_formula,
    data = df_train,
    trace = FALSE
  )
  
  prob_multinom <- predict(fit_multinom_k, newdata = df_test, type = "probs")
  prob_multinom <- as.matrix(prob_multinom)
  
  # Pad with zeros for any class absent from the training fold
  if (!all(classes %in% colnames(prob_multinom))) {
    tmp <- matrix(0, nrow = nrow(df_test), ncol = length(classes))
    colnames(tmp) <- classes
    tmp[, colnames(prob_multinom)] <- prob_multinom
    prob_multinom <- tmp
  }
  
  # Hard class predictions from the highest predicted probability
  pred_multinom <- classes[max.col(prob_multinom, ties.method = "first")]
  
  # -----------------------------------------------------------
  # Model 4: Mixed-effects multinomial logit (GLMM baseline)
  # mclogit::mblogit() fits a mixed-effects baseline-category
  # logit model with a cluster-level random intercept.
  # This is the parametric multilevel baseline (Section 3.1.2).
  # conditional = FALSE returns population-level probabilities
  # (random effects set to zero), matching the prediction
  # strategy used for GMERT-cat and GMERF-cat.
  # -----------------------------------------------------------
  fit_mblogit_k <- mclogit::mblogit(
    formula = mblogit_formula,
    random = ~ 1 | id,     # random intercept per cluster
    data = df_train
  )
  
  prob_mblogit <- predict(
    fit_mblogit_k,
    newdata = df_test,
    type = "response",
    conditional = FALSE     # population-level (marginal) prediction
  )
  prob_mblogit <- as.matrix(prob_mblogit)
  
  # Pad with zeros for any class absent from the training fold
  if (!all(classes %in% colnames(prob_mblogit))) {
    tmp <- matrix(0, nrow = nrow(df_test), ncol = length(classes))
    colnames(tmp) <- classes
    tmp[, colnames(prob_mblogit)] <- prob_mblogit
    prob_mblogit <- tmp
  }
  
  pred_mblogit <- classes[max.col(prob_mblogit, ties.method = "first")]
  
  # -----------------------------------------------------------
  # Model 5: CART (single-level tree baseline)
  # Classification tree via rpart::rpart(). No random effects;
  # treats all observations as independent. Serves as the
  # non-ensemble single-level machine-learning baseline
  # (Section 3.1.3).
  #
  # Hyperparameter choices:
  #   cp = 0.0       No cost-complexity pruning; the tree
  #                  is grown to the limits set by minsplit,
  #                  minbucket, and maxdepth.
  #   minsplit = 20  Minimum observations in a node before a
  #                  split is attempted.
  #   minbucket = 7  Minimum observations in any terminal node.
  #   maxdepth = 3   Maximum tree depth; limits model complexity
  #                  and prevents severe overfitting.
  # -----------------------------------------------------------
  fit_cart_k <- rpart::rpart(
    formula = cart_formula,
    data = df_train,
    method = "class",
    control = rpart::rpart.control(
      cp = 0.0,
      minsplit = 20,
      minbucket = 7,
      maxdepth = 3
    )
  )
  
  prob_cart <- predict(fit_cart_k, newdata = df_test, type = "prob")
  prob_cart <- as.matrix(prob_cart)
  
  # Pad with zeros for any class absent from the training fold
  if (!all(classes %in% colnames(prob_cart))) {
    tmp <- matrix(0, nrow = nrow(df_test), ncol = length(classes))
    colnames(tmp) <- classes
    tmp[, colnames(prob_cart)] <- prob_cart
    prob_cart <- tmp
  }
  
  pred_cart <- colnames(prob_cart)[max.col(prob_cart, ties.method = "first")]
  
  # -----------------------------------------------------------
  # Model 6: Random forest (single-level ensemble baseline)
  # ranger::ranger() fits a fast random forest. No random
  # effects; treats all observations as independent. Serves as
  # the ensemble single-level machine-learning baseline
  # (Section 3.1.3). probability = TRUE returns a full N x K
  # probability matrix rather than hard class predictions.
  #
  # Hyperparameter choices:
  #   num.trees = 500   Standard ensemble size balancing
  #                     variance reduction and compute cost.
  #   mtry = 3          Candidate predictors per split;
  #                     approximately sqrt(p) for p = 10.
  #   min.node.size = 10 Minimum terminal node size for
  #                     probability forests (ranger default
  #                     for probability = TRUE is 10).
  #   seed = rf_seed    Fixed for reproducibility within fold.
  # -----------------------------------------------------------
  fit_rf_k <- ranger::ranger(
    formula = rf_formula,
    data = df_train,
    num.trees = 500,
    mtry = 3,
    min.node.size = 10,
    classification = TRUE,
    probability = TRUE,
    seed = rf_seed
  )
  
  prob_rf <- predict(fit_rf_k, data = df_test)$predictions
  prob_rf <- as.matrix(prob_rf)
  
  # Pad with zeros for any class absent from the training fold
  if (!all(classes %in% colnames(prob_rf))) {
    tmp <- matrix(0, nrow = nrow(df_test), ncol = length(classes))
    colnames(tmp) <- classes
    tmp[, colnames(prob_rf)] <- prob_rf
    prob_rf <- tmp
  }
  
  pred_rf <- colnames(prob_rf)[max.col(prob_rf, ties.method = "first")]
  
  # -----------------------------------------------------------
  # Confusion matrices
  # One per model; stored in the return value for optional
  # downstream inspection (e.g. per-class breakdown).
  # -----------------------------------------------------------
  cm_gmert    <- table(predicted = pred_gmert,    actual = df_test$y)
  cm_gmerf    <- table(predicted = pred_gmerf,    actual = df_test$y)
  cm_multinom <- table(predicted = pred_multinom, actual = df_test$y)
  cm_mblogit  <- table(predicted = pred_mblogit,  actual = df_test$y)
  cm_cart     <- table(predicted = pred_cart,     actual = df_test$y)
  cm_rf       <- table(predicted = pred_rf,       actual = df_test$y)
  
  # -----------------------------------------------------------
  # Performance metrics for this fold
  # Accuracy, macro F1, weighted F1, Brier score, cross-entropy
  # are computed for all six models and collected in a flat
  # named list for easy extraction in summarise_cv_cat().
  # -----------------------------------------------------------
  metrics <- list(
    acc_gmert    = acc_fun_mc(df_test$y, pred_gmert),
    acc_gmerf    = acc_fun_mc(df_test$y, pred_gmerf),
    acc_multinom = acc_fun_mc(df_test$y, pred_multinom),
    acc_mblogit  = acc_fun_mc(df_test$y, pred_mblogit),
    acc_cart     = acc_fun_mc(df_test$y, pred_cart),
    acc_rf       = acc_fun_mc(df_test$y, pred_rf),
    
    macrof1_gmert    = macro_f1_fun(df_test$y, pred_gmert),
    macrof1_gmerf    = macro_f1_fun(df_test$y, pred_gmerf),
    macrof1_multinom = macro_f1_fun(df_test$y, pred_multinom),
    macrof1_mblogit  = macro_f1_fun(df_test$y, pred_mblogit),
    macrof1_cart     = macro_f1_fun(df_test$y, pred_cart),
    macrof1_rf       = macro_f1_fun(df_test$y, pred_rf),
    
    wf1_gmert    = weighted_f1_fun(df_test$y, pred_gmert),
    wf1_gmerf    = weighted_f1_fun(df_test$y, pred_gmerf),
    wf1_multinom = weighted_f1_fun(df_test$y, pred_multinom),
    wf1_mblogit  = weighted_f1_fun(df_test$y, pred_mblogit),
    wf1_cart     = weighted_f1_fun(df_test$y, pred_cart),
    wf1_rf       = weighted_f1_fun(df_test$y, pred_rf),
    
    brier_gmert    = brier_score(prob_gmert, df_test$y, classes),
    brier_gmerf    = brier_score(prob_gmerf, df_test$y, classes),
    brier_multinom = brier_score(prob_multinom, df_test$y, classes),
    brier_mblogit  = brier_score(prob_mblogit, df_test$y, classes),
    brier_cart     = brier_score(prob_cart, df_test$y, classes),
    brier_rf       = brier_score(prob_rf, df_test$y, classes),
    
    ce_gmert    = cross_entropy(prob_gmert, df_test$y, classes),
    ce_gmerf    = cross_entropy(prob_gmerf, df_test$y, classes),
    ce_multinom = cross_entropy(prob_multinom, df_test$y, classes),
    ce_mblogit  = cross_entropy(prob_mblogit, df_test$y, classes),
    ce_cart     = cross_entropy(prob_cart, df_test$y, classes),
    ce_rf       = cross_entropy(prob_rf, df_test$y, classes)
  )
  
  list(
    fold_id   = k,
    train_idx = train_idx,
    test_idx  = test_idx,
    metrics   = metrics,
    df_train  = df_train,
    df_test   = df_test
  )
}


# ===============================================================
# SECTION 3: Cross-validation summary
# ===============================================================

# ---------------------------------------------------------------
# summarise_cv_cat()
#
# Aggregates per-fold metrics returned by cv_one_fold_cat()
# into a tidy summary table (mean ± SD across folds), matching
# the reporting format used in the thesis (Section 4.1.2 and
# Appendix Table 10).
#
# Arguments:
#   cv_results  List of K outputs from cv_one_fold_cat().
#   classes     Character vector of class names (currently
#               unused in the summary but kept for potential
#               per-class breakdowns in future extensions).
#
# Returns:
#   A list with two elements:
#     $metrics_df   K x (n_metrics) data frame with one row
#                   per fold and one column per model-metric
#                   combination.
#     $cv_summary   6 x (n_metrics * 2 + 1) data frame with
#                   one row per model and columns for mean and
#                   SD of each metric across folds.
# ---------------------------------------------------------------
summarise_cv_cat <- function(cv_results, classes) {
  
  # Stack per-fold metrics into a single data frame (K rows)
  metrics_df <- do.call(rbind, lapply(cv_results, function(x) {
    data.frame(
      fold = as.double(x$fold_id),
      
      acc_gmert    = x$metrics$acc_gmert,
      acc_gmerf    = x$metrics$acc_gmerf,
      acc_multinom = x$metrics$acc_multinom,
      acc_mblogit  = x$metrics$acc_mblogit,
      acc_cart     = x$metrics$acc_cart,
      acc_rf       = x$metrics$acc_rf,
      
      macrof1_gmert    = x$metrics$macrof1_gmert,
      macrof1_gmerf    = x$metrics$macrof1_gmerf,
      macrof1_multinom = x$metrics$macrof1_multinom,
      macrof1_mblogit  = x$metrics$macrof1_mblogit,
      macrof1_cart     = x$metrics$macrof1_cart,
      macrof1_rf       = x$metrics$macrof1_rf,
      
      wf1_gmert    = x$metrics$wf1_gmert,
      wf1_gmerf    = x$metrics$wf1_gmerf,
      wf1_multinom = x$metrics$wf1_multinom,
      wf1_mblogit  = x$metrics$wf1_mblogit,
      wf1_cart     = x$metrics$wf1_cart,
      wf1_rf       = x$metrics$wf1_rf,
      
      brier_gmert    = x$metrics$brier_gmert,
      brier_gmerf    = x$metrics$brier_gmerf,
      brier_multinom = x$metrics$brier_multinom,
      brier_mblogit  = x$metrics$brier_mblogit,
      brier_cart     = x$metrics$brier_cart,
      brier_rf       = x$metrics$brier_rf,
      
      ce_gmert    = x$metrics$ce_gmert,
      ce_gmerf    = x$metrics$ce_gmerf,
      ce_multinom = x$metrics$ce_multinom,
      ce_mblogit  = x$metrics$ce_mblogit,
      ce_cart     = x$metrics$ce_cart,
      ce_rf       = x$metrics$ce_rf
    )
  }))
  
  # Compute mean and SD across folds for each model and metric
  cv_summary <- data.frame(
    model = c("GMERT-cat", "GMERF-cat", "Multinom", "MBLOGIT", "CART", "RF"),
    
    acc_mean = c(mean(metrics_df$acc_gmert),
                 mean(metrics_df$acc_gmerf),
                 mean(metrics_df$acc_multinom),
                 mean(metrics_df$acc_mblogit),
                 mean(metrics_df$acc_cart),
                 mean(metrics_df$acc_rf)),
    acc_sd = c(sd(metrics_df$acc_gmert),
               sd(metrics_df$acc_gmerf),
               sd(metrics_df$acc_multinom),
               sd(metrics_df$acc_mblogit),
               sd(metrics_df$acc_cart),
               sd(metrics_df$acc_rf)),
    
    macrof1_mean = c(mean(metrics_df$macrof1_gmert),
                     mean(metrics_df$macrof1_gmerf),
                     mean(metrics_df$macrof1_multinom),
                     mean(metrics_df$macrof1_mblogit),
                     mean(metrics_df$macrof1_cart),
                     mean(metrics_df$macrof1_rf)),
    macrof1_sd = c(sd(metrics_df$macrof1_gmert),
                   sd(metrics_df$macrof1_gmerf),
                   sd(metrics_df$macrof1_multinom),
                   sd(metrics_df$macrof1_mblogit),
                   sd(metrics_df$macrof1_cart),
                   sd(metrics_df$macrof1_rf)),
    
    wf1_mean = c(mean(metrics_df$wf1_gmert),
                 mean(metrics_df$wf1_gmerf),
                 mean(metrics_df$wf1_multinom),
                 mean(metrics_df$wf1_mblogit),
                 mean(metrics_df$wf1_cart),
                 mean(metrics_df$wf1_rf)),
    wf1_sd = c(sd(metrics_df$wf1_gmert),
               sd(metrics_df$wf1_gmerf),
               sd(metrics_df$wf1_multinom),
               sd(metrics_df$wf1_mblogit),
               sd(metrics_df$wf1_cart),
               sd(metrics_df$wf1_rf)),
    
    brier_mean = c(mean(metrics_df$brier_gmert),
                   mean(metrics_df$brier_gmerf),
                   mean(metrics_df$brier_multinom),
                   mean(metrics_df$brier_mblogit),
                   mean(metrics_df$brier_cart),
                   mean(metrics_df$brier_rf)),
    brier_sd = c(sd(metrics_df$brier_gmert),
                 sd(metrics_df$brier_gmerf),
                 sd(metrics_df$brier_multinom),
                 sd(metrics_df$brier_mblogit),
                 sd(metrics_df$brier_cart),
                 sd(metrics_df$brier_rf)),
    
    ce_mean = c(mean(metrics_df$ce_gmert),
                mean(metrics_df$ce_gmerf),
                mean(metrics_df$ce_multinom),
                mean(metrics_df$ce_mblogit),
                mean(metrics_df$ce_cart),
                mean(metrics_df$ce_rf)),
    ce_sd = c(sd(metrics_df$ce_gmert),
              sd(metrics_df$ce_gmerf),
              sd(metrics_df$ce_multinom),
              sd(metrics_df$ce_mblogit),
              sd(metrics_df$ce_cart),
              sd(metrics_df$ce_rf))
  )
  
  list(
    metrics_df = metrics_df,
    cv_summary = cv_summary
  )
}


# ===============================================================
# SECTION 4: Parallel cluster setup
# ===============================================================

# Determine the number of workers: one per fold up to
# (available cores - 1), leaving one core free for the main
# process. Using all cores risks starving the OS scheduler.
n_workers <- max(1L, min(K_folds, parallel::detectCores() - 1L))

# PSOCK cluster: spawns independent R sessions connected via
# sockets. Works on all platforms including Windows.
cl <- parallel::makeCluster(n_workers, type = "PSOCK")

# Load required packages on every worker.
# mlml must be installed and available on all worker sessions;
# it provides fit_gmert_cat(), fit_gmerf_cat(), and their
# predict() methods.
parallel::clusterEvalQ(cl, {
  library(mlml)
  library(nnet)
  library(mclogit)
  library(rpart)
  library(ranger)
  library(tidyverse)
  
  # Disable thread-level parallelism inside each worker to
  # prevent over-subscription: each worker already runs in its
  # own process, so multi-threading within a worker would
  # create n_workers * n_threads active threads simultaneously,
  # exceeding the number of physical cores and degrading
  # performance through context switching.
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OMP_THREAD_LIMIT = "1"
  )
})

# Export all functions defined in this script to each worker.
# Worker sessions start with a clean environment, so every
# function used inside cv_one_fold_cat() must be exported
# explicitly.
parallel::clusterExport(
  cl,
  varlist = c(
    "cv_one_fold_cat",
    "make_cluster_folds",
    "acc_fun_mc",
    "f1_per_class",
    "macro_f1_fun",
    "weighted_f1_fun",
    "brier_score",
    "cross_entropy",
    "summarise_cv_cat"
  ),
  envir = environment()
)


# ===============================================================
# SECTION 5: Main CV loop — one scenario at a time
# ===============================================================

results <- list()

for (nm in names(df_list)) {
  
  df  <- df_list[[nm]]        # data frame for this scenario
  cfg <- scenario_cfg[[nm]]   # model formulas and hyperparameters
  
  # Construct folds once per scenario and share them with all
  # workers so every model is evaluated on identical splits.
  folds <- make_cluster_folds(df, K = K_folds, seed = seed_folds)
  
  # Export the scenario-specific data and folds to all workers.
  # These must be re-exported at every iteration because df and
  # folds change with each scenario.
  parallel::clusterExport(cl, varlist = c("df", "folds"), envir = environment())
  
  # Set a reproducible parallel RNG stream. clusterSetRNGStream
  # uses L'Ecuyer-CMRG streams so each worker gets a
  # statistically independent stream derived from iseed.
  parallel::clusterSetRNGStream(cl, iseed = seed_cluster)
  
  # Dispatch one fold per worker. parLapply() blocks until all
  # K folds are complete, then returns a list of K results.
  # Model formulas and hyperparameter lists are passed as
  # additional arguments to avoid relying on lexical scoping
  # across the socket boundary.
  cv_res <- parallel::parLapply(
    cl,
    X = seq_len(K_folds),
    fun = function(k, multinom_formula, mblogit_formula, cart_formula,
                   rf_formula, gmert_args, gmerf_args) {
      cv_one_fold_cat(
        k             = k,
        df            = df,
        folds         = folds,
        gmert_args    = gmert_args,
        gmerf_args    = gmerf_args,
        multinom_formula = multinom_formula,
        mblogit_formula  = mblogit_formula,
        cart_formula     = cart_formula,
        rf_formula       = rf_formula,
        rf_seed          = 42
      )
    },
    multinom_formula = cfg$multinom_formula,
    mblogit_formula  = cfg$mblogit_formula,
    cart_formula     = cfg$cart_formula,
    rf_formula       = cfg$rf_formula,
    gmert_args       = cfg$gmert_args,
    gmerf_args       = cfg$gmerf_args
  )
  
  # Aggregate per-fold results into mean ± SD summary
  summ <- summarise_cv_cat(cv_res, classes = levels(df$y))
  
  results[[nm]] <- list(
    metrics_df = summ$metrics_df,
    cv_summary = summ$cv_summary
  )
  
  cat("Done:", nm, "\n")
}

# Shut down the worker processes cleanly
parallel::stopCluster(cl)