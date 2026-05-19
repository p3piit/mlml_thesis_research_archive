# ===============================================================
# Helpers

# cluster-aware folds
make_cluster_folds <- function(df, K = 10, seed = 42, split_by_cluster = FALSE) {
  set.seed(seed)
  folds <- vector("list", K)
  for (k in seq_len(K)) folds[[k]] <- integer(0)
  ids <- unique(as.character(df$id))

  if (split_by_cluster) {
    shuffled_ids <- sample(ids)
    assignments <- rep(seq_len(K), length.out = length(shuffled_ids))
    for (i in seq_along(shuffled_ids)) {
      k <- assignments[i]
      folds[[k]] <- c(folds[[k]], which(as.character(df$id) == shuffled_ids[i]))
    }
  } else {
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

# multiclass accuracy
acc_fun_mc <- function(y_true, y_pred) {
  mean(as.character(y_true) == as.character(y_pred))
}

# per-class F1
f1_per_class <- function(y_true, y_pred, classes = levels(y_true)) {
  sapply(classes, function(cl) {
    tp <- sum(y_true == cl & y_pred == cl)
    fp <- sum(y_true != cl & y_pred == cl)
    fn <- sum(y_true == cl & y_pred != cl)

    if ((tp + fp) == 0 || (tp + fn) == 0) return(NA_real_)

    prec <- tp / (tp + fp)
    rec  <- tp / (tp + fn)

    if ((prec + rec) == 0) return(NA_real_)
    2 * prec * rec / (prec + rec)
  })
}

# macro F1
macro_f1_fun <- function(y_true, y_pred, classes = levels(y_true)) {
  mean(f1_per_class(y_true, y_pred, classes = classes), na.rm = TRUE)
}

# weighted F1
weighted_f1_fun <- function(y_true, y_pred, classes = levels(y_true)) {
  f1s <- f1_per_class(y_true, y_pred, classes = classes)
  wts <- as.numeric(table(factor(y_true, levels = classes)))
  sum(f1s * wts, na.rm = TRUE) / sum(wts)
}

# Brier score for K classes
brier_score<- function(prob_mat, y_true, classes) {
  Y_oh <- model.matrix(~ 0 + factor(df_test$y, levels = classes))
  mean(rowSums((prob_mat - Y_oh)^2))
}

# Cross-entropy
cross_entropy <- function(prob_mat, y_true, classes) {
  idx <- match(as.character(y_true), classes)
  probs <- prob_mat[cbind(seq_along(idx), idx)]
  -mean(log(pmax(probs, 1e-15)))
}

# ===============================================================
# One fold
cv_one_fold_cat <- function(k, df, folds,
                            gmert_args = list(),
                            gmerf_args = list(),
                            multinom_formula,
                            mblogit_formula,
                            cart_formula,
                            rf_formula,
                            rf_seed = 42) {

  test_idx  <- folds[[k]]
  train_idx <- setdiff(seq_len(nrow(df)), test_idx)

  df_train <- df[train_idx, ]
  df_test  <- df[test_idx, ]

  classes <- levels(df$y)

  # -------------------------------------------------------------
  # GMERT-cat
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
    random_effect = NULL,
    id = "id"
  )

  prob_gmert <- predict_gmert_cat(
    fit = fit_gmert_k,
    new_df = df_test,
    random_effect = NULL,
    id = "id",
    prob_saved = TRUE
  )

  # -------------------------------------------------------------
  # GMERF-cat
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

  # -------------------------------------------------------------
  # Multinomial logit
  fit_multinom_k <- nnet::multinom(
    formula = multinom_formula,
    data = df_train,
    trace = FALSE
  )

  prob_multinom <- predict(fit_multinom_k, newdata = df_test, type = "probs")
  prob_multinom <- as.matrix(prob_multinom)

  if (!all(classes %in% colnames(prob_multinom))) {
    tmp <- matrix(0, nrow = nrow(df_test), ncol = length(classes))
    colnames(tmp) <- classes
    tmp[, colnames(prob_multinom)] <- prob_multinom
    prob_multinom <- tmp
  }

  pred_multinom <- classes[max.col(prob_multinom, ties.method = "first")]

  # -------------------------------------------------------------
  # Mixed-effects multinomial logit
  fit_mblogit_k <- mclogit::mblogit(
    formula = mblogit_formula,
    random = ~ 1 | id,
    data = df_train
  )

  prob_mblogit <- predict(
    fit_mblogit_k,
    newdata = df_test,
    type = "response",
    conditional = FALSE
  )
  prob_mblogit <- as.matrix(prob_mblogit)

  if (!all(classes %in% colnames(prob_mblogit))) {
    tmp <- matrix(0, nrow = nrow(df_test), ncol = length(classes))
    colnames(tmp) <- classes
    tmp[, colnames(prob_mblogit)] <- prob_mblogit
    prob_mblogit <- tmp
  }

  pred_mblogit <- classes[max.col(prob_mblogit, ties.method = "first")]

  # -------------------------------------------------------------
  # CART
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

  if (!all(classes %in% colnames(prob_cart))) {
    tmp <- matrix(0, nrow = nrow(df_test), ncol = length(classes))
    colnames(tmp) <- classes
    tmp[, colnames(prob_cart)] <- prob_cart
    prob_cart <- tmp
  }

  pred_cart <- colnames(prob_cart)[max.col(prob_cart, ties.method = "first")]

  # -------------------------------------------------------------
  # Random forest
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

  if (!all(classes %in% colnames(prob_rf))) {
    tmp <- matrix(0, nrow = nrow(df_test), ncol = length(classes))
    colnames(tmp) <- classes
    tmp[, colnames(prob_rf)] <- prob_rf
    prob_rf <- tmp
  }

  pred_rf <- colnames(prob_rf)[max.col(prob_rf, ties.method = "first")]

  # -------------------------------------------------------------
  # Confusion matrices
  cm_gmert    <- table(predicted = pred_gmert,    actual = df_test$y)
  cm_gmerf    <- table(predicted = pred_gmerf,    actual = df_test$y)
  cm_multinom <- table(predicted = pred_multinom, actual = df_test$y)
  cm_mblogit  <- table(predicted = pred_mblogit,  actual = df_test$y)
  cm_cart     <- table(predicted = pred_cart,     actual = df_test$y)
  cm_rf       <- table(predicted = pred_rf,       actual = df_test$y)

  # -------------------------------------------------------------
  # Metrics
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
    ce_rf       = cross_entropy(prob_rf, df_test$y, classes
  ))

  list(
    fold_id = k,
    train_idx = train_idx,
    test_idx = test_idx,
    metrics = metrics,
    df_train = df_train,
    df_test = df_test
  )
}

# ===============================================================
# Summaries
summarise_cv_cat <- function(cv_results, classes) {
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
# Parallel setup
n_workers <- max(1L, min(K_folds, parallel::detectCores() - 1L))
cl <- parallel::makeCluster(n_workers, type = "PSOCK")

parallel::clusterEvalQ(cl, {
  library(mlml)
  library(nnet)
  library(mclogit)
  library(rpart)
  library(ranger)
  library(tidyverse)
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OMP_THREAD_LIMIT = "1"
  )
})

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
# Run CV for each scenario
results <- list()

for (nm in names(df_list)) {

  df <- df_list[[nm]]
  cfg <- scenario_cfg[[nm]]

  folds <- make_cluster_folds(df, K = K_folds, seed = seed_folds)

  parallel::clusterExport(cl, varlist = c("df", "folds"), envir = environment())
  parallel::clusterSetRNGStream(cl, iseed = seed_cluster)

  cv_res <- parallel::parLapply(
    cl,
    X = seq_len(K_folds),
    fun = function(k, multinom_formula, mblogit_formula, cart_formula, rf_formula, gmert_args, gmerf_args) {
      cv_one_fold_cat(
        k = k,
        df = df,
        folds = folds,
        gmert_args = gmert_args,
        gmerf_args = gmerf_args,
        multinom_formula = multinom_formula,
        mblogit_formula = mblogit_formula,
        cart_formula = cart_formula,
        rf_formula = rf_formula,
        rf_seed = 42
      )
    },
    multinom_formula = cfg$multinom_formula,
    mblogit_formula = cfg$mblogit_formula,
    cart_formula = cfg$cart_formula,
    rf_formula = cfg$rf_formula,
    gmert_args = cfg$gmert_args,
    gmerf_args = cfg$gmerf_args
  )

  summ <- summarise_cv_cat(cv_res, classes = levels(df$y))

  results[[nm]] <- list(
    metrics_df = summ$metrics_df,
    cv_summary = summ$cv_summary
  )

  cat("Done:", nm, "\n")
}

parallel::stopCluster(cl)

