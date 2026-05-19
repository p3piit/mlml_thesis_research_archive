# ===============================================================
# table.R
#
# Generates formatted summary tables for the multinomial
# simulation study (Appendix Table 10 and related tables).
#
# For each of three evaluation metrics, a 6 × 4 table is
# produced with models as rows and simulation scenarios as
# columns. Each cell contains "mean ± SD" across the 10 CV
# folds, matching the reporting format used throughout the
# thesis (Section 3.4).
#
# Three metrics are tabulated:
#   Accuracy       — overall classification accuracy
#   Macro F1       — unweighted mean F1 across classes
#   Weighted F1    — frequency-weighted mean F1
#
# Each table is output in two formats:
#   LaTeX   via xtable::xtable(), for direct inclusion in the
#           thesis manuscript.
#   PNG     via gridExtra::grid.table(), for the research
#           archive and visual inspection.
#
# Input:
#   results   Named list produced by the analysis script;
#             must contain elements s1, s2, s3, s4, each
#             with a $cv_summary data frame.
#
# Output (written to tables/):
#   acc_table.png, f1_table.png, wf1_table.png
#
# LaTeX output is printed to the console and can be redirected
# or captured with print(xtable(...), file = "tables/....tex").
#
# Dependencies: xtable, gridExtra
# ===============================================================


# ---------------------------------------------------------------
# make_metric_table()
#
# Builds a formatted 6 × 4 character matrix (returned as a
# data frame) containing "mean ± SD" strings for one metric
# across all models and scenarios.
#
# Arguments:
#   results      Named list of scenario results; each element
#                must contain $cv_summary with columns `model`,
#                and the two metric columns named below.
#   metric_mean  String; column name of the mean metric in
#                cv_summary (e.g. "acc_mean").
#   metric_sd    String; column name of the SD metric in
#                cv_summary (e.g. "acc_sd").
#   scen_names   Character vector of keys into `results`.
#                Default: c("s1","s2","s3","s4").
#   scen_labels  Character vector of column labels for the
#                output table. Default: the four scenario names
#                used throughout the thesis.
#
# Returns:
#   A data frame with nrow = 6 (one per model) and
#   ncol = 4 (one per scenario). Each cell is a character
#   string of the form "x.xxx ± x.xxx".
# ---------------------------------------------------------------
make_metric_table <- function(results, metric_mean, metric_sd,
                              scen_names  = c("s1","s2","s3","s4"),
                              scen_labels = c("Baseline","Nonlinear",
                                              "Interaction","Multicollinearity")) {
  
  # Mapping from internal cv_summary model identifiers to the
  # display names used in the thesis tables (Appendix Table 10)
  models_map <- c("GMERT-cat" = "GMERT",
                  "GMERF-cat" = "GMERF",
                  "Multinom"   = "Logit",
                  "MBLOGIT"    = "GLMM",
                  "CART"       = "CART",
                  "RF"         = "RF")
  
  # Pre-allocate output matrix: rows = models, columns = scenarios
  out <- matrix(nrow = length(models_map), ncol = length(scen_names))
  rownames(out) <- unname(models_map)
  colnames(out) <- scen_labels
  
  # Fill one cell at a time by looking up each model row within
  # the cv_summary of each scenario
  for (j in seq_along(scen_names)) {
    cv <- results[[scen_names[j]]]$cv_summary
    
    for (i in seq_along(models_map)) {
      model_key <- names(models_map)[i]
      
      mm  <- cv[cv$model == model_key, metric_mean]
      sdv <- cv[cv$model == model_key, metric_sd]
      
      # Format as "x.xxx ± x.xxx" — three decimal places matches
      # the precision used in the thesis (Appendix Table 10)
      out[i, j] <- sprintf("%.3f \u00b1 %.3f", mm, sdv)
    }
  }
  
  as.data.frame(out, stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------
if (!dir.exists("tables")) dir.create("tables", recursive = TRUE)


# ---------------------------------------------------------------
# Build one table per metric
# ---------------------------------------------------------------
acc_tab <- make_metric_table(results, "acc_mean",     "acc_sd")
f1_tab  <- make_metric_table(results, "macrof1_mean", "macrof1_sd")
wf1_tab <- make_metric_table(results, "wf1_mean",     "wf1_sd")


# ---------------------------------------------------------------
# LaTeX output via xtable
# ---------------------------------------------------------------
xtable(acc_tab,
       caption = "Mean accuracy (mean \u00b1 SD) across scenarios",
       label   = "tab:acc_scores")

xtable(f1_tab,
       caption = "Mean macro F1-score (mean \u00b1 SD) across scenarios",
       label   = "tab:f1_scores")

xtable(wf1_tab,
       caption = "Mean weighted F1-score (mean \u00b1 SD) across scenarios",
       label   = "tab:wf1_scores")


# ---------------------------------------------------------------
# PNG output via gridExtra::grid.table()
# ---------------------------------------------------------------
png("tables/acc_table.png",  width = 600, height = 400)
grid.table(acc_tab)
dev.off()

png("tables/f1_table.png",   width = 600, height = 400)
grid.table(f1_tab)
dev.off()

png("tables/wf1_table.png",  width = 600, height = 400)
grid.table(wf1_tab)
dev.off()