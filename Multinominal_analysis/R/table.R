# ===============================================================
# In this script we generate the table needed for the research report
# ===============================================================


# build helper to format mean ± sd across scenarios
make_metric_table <- function(results, metric_mean, metric_sd,
                              scen_names = c("s1","s2","s3","s4"),
                              scen_labels = c("Baseline","Nonlinear","Interaction","Multicollinearity")) {
  models_map <- c("GMERT-cat" = "GMERT",
                  "GMERF-cat" = "GMERF",
                  "Multinom"   = "Logit",
                  "MBLOGIT"    = "GLMM",
                  "CART"       = "CART",
                  "RF"         = "RF")
  out <- matrix(nrow = length(models_map), ncol = length(scen_names))
  rownames(out) <- unname(models_map)
  colnames(out) <- scen_labels

  for (j in seq_along(scen_names)) {
    cv <- results[[scen_names[j]]]$cv_summary
    for (i in seq_along(models_map)) {
      model_key <- names(models_map)[i]
      mm <- cv[cv$model == model_key, metric_mean]
      sdv <- cv[cv$model == model_key, metric_sd]
      out[i, j] <- sprintf("%.3f ± %.3f", mm, sdv)
    }
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

# ensure output dir
if (!dir.exists("tables")) dir.create("tables", recursive = TRUE)

# tables for accuracy, macro F1, weighted F1
acc_tab <- make_metric_table(results, "acc_mean", "acc_sd")
f1_tab  <- make_metric_table(results, "macrof1_mean", "macrof1_sd")
wf1_tab <- make_metric_table(results, "wf1_mean", "wf1_sd")
brier_tab <- make_metric_table(results, "brier_mean", "brier_sd")
ce_tab <- make_metric_table(results, "ce_mean", "ce_sd")


xtable(acc_tab, caption = "Mean accuracy (mean ± sd) across scenarios", label = "tab:acc_scores")
xtable(f1_tab,  caption = "Mean macro F1-score (mean ± sd) across scenarios", label = "tab:f1_scores")
xtable(wf1_tab, caption = "Mean weighted F1-score (mean ± sd) across scenarios", label = "tab:wf1_scores")
xtable(brier_tab, caption = "Mean Brier score (mean ± sd) across scenarios", label = "tab:brier_scores")
xtable(ce_tab, caption = "Mean cross-entropy (mean ± sd) across scenarios", label = "tab:brier_scores")

png("tables/acc_table.png", width = 600, height = 400)
grid.table(acc_tab)
dev.off()

png("tables/f1_table.png", width = 600, height = 400)
grid.table(f1_tab)
dev.off()

png("tables/wf1_table.png", width = 600, height = 400)
grid.table(wf1_tab)
dev.off()

png("tables/brier_table.png", width = 600, height = 400)
grid.table(brier_tab)
dev.off()

png("tables/ce_table.png", width = 600, height = 400)
grid.table(ce_tab)
dev.off()