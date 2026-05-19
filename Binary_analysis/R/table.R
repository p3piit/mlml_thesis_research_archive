# ===============================================================
# table.R
#
# Generates formatted summary tables for the binary simulation
# study (Appendix Table 9 and related tables).
#
# For each of four evaluation metrics, a 6 × 4 table is
# produced with models as rows and simulation scenarios as
# columns. Each cell contains "mean ± SD" across the 10 CV
# folds, matching the reporting format used throughout the
# thesis (Section 3.4).
#
# Four metrics are tabulated:
#   Accuracy        — overall classification accuracy
#   F1 Majority     — F1 score for the majority class
#   F1 Minority     — F1 score for the minority class
#   Bias            — relative prediction bias for class 1
#
# Each table is output in two formats:
#   LaTeX   via xtable::xtable(), for direct inclusion in the
#           thesis manuscript.
#   PNG     via gridExtra::grid.table(), for the research
#           archive and visual inspection.
#
# Input:
#   results   Named list produced by analysis_binary.R;
#             must contain elements s1, s2, s3, s4, each
#             with a $cv_summary data frame.
#
# Output (written to tables/):
#   acc_table.png, f1_mag_table.png,
#   f1_min_table.png, bias_table.png
#
# LaTeX output is printed to the console and can be redirected
# with print(xtable(...), file = "tables/....tex").
#
# Dependencies: xtable, gridExtra
# ===============================================================


# ---------------------------------------------------------------
# make_metric_table_bin()
#
# Builds a formatted 6 × 4 character matrix (returned as a
# data frame) containing "mean ± SD" strings for one metric
# across all models and scenarios.
#
# Row order follows models_map below, which maps internal model
# identifiers (as stored in cv_summary$model) to the display
# names used in the thesis tables (Appendix Table 9).
#
# Arguments:
#   results      Named list of scenario results; each element
#                must contain $cv_summary with columns `model`
#                and the two metric columns named below.
#   metric_mean  String; column name of the mean metric in
#                cv_summary (e.g. "acc_mean").
#   metric_sd    String; column name of the SD metric in
#                cv_summary (e.g. "acc_sd").
#   scen_names   Character vector of keys into `results`.
#                Default: c("s1","s2","s3","s4").
#   scen_labels  Character vector of column labels for the
#                output table.
#
# Returns:
#   A data frame with nrow = 6 (one per model) and
#   ncol = 4 (one per scenario). Each cell is a character
#   string of the form "x.xxx ± x.xxx".
# ---------------------------------------------------------------
make_metric_table_bin <- function(results, metric_mean, metric_sd,
                                  scen_names  = c("s1","s2","s3","s4"),
                                  scen_labels = c("Baseline","Nonlinear",
                                                  "Interaction","Multicollinearity")) {
  
  # Mapping from internal cv_summary model identifiers to the
  # display names used in the thesis tables (Appendix Table 9)
  models_map <- c("GMERT" = "GMERT",
                  "GMERF" = "GMERF",
                  "GLMM"  = "GLMM",
                  "Logit" = "Logit",
                  "CART"  = "CART",
                  "RF"    = "RF")
  
  # Pre-allocate output matrix: rows = models, columns = scenarios
  out <- matrix(nrow = length(models_map), ncol = length(scen_names))
  rownames(out) <- unname(models_map)
  colnames(out) <- scen_labels
  
  for (j in seq_along(scen_names)) {
    cv <- results[[scen_names[j]]]$cv_summary
    
    for (i in seq_along(models_map)) {
      model_key <- names(models_map)[i]
      
      mm  <- cv[cv$model == model_key, metric_mean]
      sdv <- cv[cv$model == model_key, metric_sd]
      
      # Format as "x.xxx ± x.xxx" — three decimal places matches
      # the precision used in the thesis (Appendix Table 9)
      out[i, j] <- sprintf("%.3f \u00b1 %.3f", mm, sdv)
    }
  }
  
  as.data.frame(out, stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------
if (!dir.exists(here("Binary_analysis", "tables"))) {
  dir.create(here("Binary_analysis", "tables"), recursive = TRUE)
}

# ---------------------------------------------------------------
# Build one table per metric
# ---------------------------------------------------------------
acc_tab_bin     <- make_metric_table_bin(results, "acc_mean",     "acc_sd")
f1_mag_tab_bin  <- make_metric_table_bin(results, "f1_mag_mean",  "f1_mag_sd")
f1_min_tab_bin  <- make_metric_table_bin(results, "f1_min_mean",  "f1_min_sd")
bias_tab_bin    <- make_metric_table_bin(results, "bias_mean",    "bias_sd")


# ---------------------------------------------------------------
# LaTeX output via xtable
# ---------------------------------------------------------------
xtable(acc_tab_bin,
       caption = "Mean accuracy (mean \u00b1 SD) across scenarios — binary simulation",
       label   = "tab:acc_scores_bin")

xtable(f1_mag_tab_bin,
       caption = "Mean majority-class F1 (mean \u00b1 SD) across scenarios — binary simulation",
       label   = "tab:f1_mag_scores_bin")

xtable(f1_min_tab_bin,
       caption = "Mean minority-class F1 (mean \u00b1 SD) across scenarios — binary simulation",
       label   = "tab:f1_min_scores_bin")

xtable(bias_tab_bin,
       caption = "Mean relative bias (mean \u00b1 SD) across scenarios — binary simulation",
       label   = "tab:bias_scores_bin")


# ---------------------------------------------------------------
# PNG output via gridExtra::grid.table()
# ---------------------------------------------------------------
png(here("Binary_analysis", "tables","acc_table.png"),    
    width = 600, 
    height = 400)
grid.table(acc_tab_bin)
dev.off()

png(here("Binary_analysis", "tables","f1_mag_table.png"), 
    width = 600, 
    height = 400)
grid.table(f1_mag_tab_bin)
dev.off()

png(here("Binary_analysis", "tables","f1_min_table.png"), 
    width = 600, 
    height = 400)
grid.table(f1_min_tab_bin)
dev.off()

png(here("Binary_analysis", "tables","bias_table.png"),   
    width = 600, 
    height = 400)
grid.table(bias_tab_bin)
dev.off()


