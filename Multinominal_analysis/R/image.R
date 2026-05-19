# ===============================================================
# image.R
#
# Builds the summary visualization for the multinomial
# simulation study (Figure 2 of the thesis).
#
# Takes the per-scenario CV summaries produced by the analysis
# script (stored in `results`), combines them into a single
# long-format data frame, and produces a faceted dot-and-
# errorbar plot showing mean ± SD across folds for three
# metrics (Accuracy, Macro F1, Weighted F1) broken down by
# model family (GLM, CART, RF) and nesting structure
# (With / Without Nesting) across the four simulation
# scenarios.
#
# Input:
#   results   Named list produced by the analysis script;
#             must contain elements s1, s2, s3, s4, each
#             with a $cv_summary data frame (one row per
#             model, columns: model, acc_mean, acc_sd,
#             macrof1_mean, macrof1_sd, wf1_mean, wf1_sd).
#
# Output:
#   images/graph.png   300 dpi PNG, 8 × 6 inches.
#
# Dependencies: dplyr, tidyr, ggplot2 (all via tidyverse)
# ===============================================================


# ---------------------------------------------------------------
# SECTION 1: Build unified long-format dataset
#
# Stack the four scenario summaries, add scenario labels, recode
# model names and nesting indicators, then pivot to long format
# so that each row represents one (scenario, model, nesting,
# metric) combination. This is the tidy structure expected by
# ggplot2.
# ---------------------------------------------------------------
plot_data <- results$s1$cv_summary %>%
  rbind(results$s2$cv_summary) %>%
  rbind(results$s3$cv_summary) %>%
  rbind(results$s4$cv_summary) %>%
  
  mutate(
    # Scenario labels: 6 models per scenario, stacked in order
    # s1 = Baseline, s2 = Nonlinear, s3 = Interaction, s4 = Multicollinearity
    scenario = c(rep("Baseline",         6),
                 rep("Nonlinear",        6),
                 rep("Interaction",      6),
                 rep("Multicollinearity",6)),
    
    # Nesting indicator: multilevel models (GMERT-cat, GMERF-cat,
    # MBLOGIT) are "With Nesting"; single-level models (Multinom,
    # CART, RF) are "Without Nesting"
    nesting = ifelse(model %in% c("GMERT-cat", "MBLOGIT", "GMERF-cat"),
                     "With Nesting", "Without Nesting"),
    
    # Collapse the six model names into three display families
    # so the x-axis groups parametric (GLM), tree (CART), and
    # ensemble (RF) approaches regardless of nesting — nesting
    # is encoded separately via the shape aesthetic
    model = ifelse(model %in% c("MBLOGIT", "Multinom"), "GLM",
                   ifelse(model %in% c("GMERT-cat", "CART"),  "CART", "RF")),
    
    # Convert nesting to an ordered factor to control legend
    # and shape ordering in ggplot2
    nesting = factor(
      nesting,
      levels = c("Without Nesting", "With Nesting"),
      labels = c("Without Nesting", "With Nesting")
    )
  ) %>%
  
  # Pivot the three mean-metric columns to long format.
  # After this step each original row becomes three rows,
  # one per metric. The SD columns (acc_sd, macrof1_sd,
  # wf1_sd) remain wide and are matched to their metric
  # in the next mutate() step.
  pivot_longer(
    cols      = c(acc_mean, macrof1_mean, wf1_mean),
    names_to  = "metric",
    values_to = "mean"
  ) %>%
  
  mutate(
    # Attach the correct SD to each pivoted row by matching
    # the original metric name. case_when() is used rather
    # than a join because the SD columns are still present
    # as wide columns after pivot_longer().
    sd = case_when(
      metric == "acc_mean"     ~ acc_sd,
      metric == "macrof1_mean" ~ macrof1_sd,
      metric == "wf1_mean"     ~ wf1_sd
    ),
    
    # Replace internal column-name strings with display labels
    # used in the plot legend
    metric = recode(metric,
                    acc_mean     = "Accuracy",
                    macrof1_mean = "Macro F1",
                    wf1_mean     = "Weighted F1")
  )


# ---------------------------------------------------------------
# SECTION 2: Faceted dot-and-errorbar plot
#
# One panel per simulation scenario (2 × 2 grid).
# Within each panel:
#   x-axis   Model family (GLM, CART, RF)
#   y-axis   Mean metric value ± 1 SD across 10 CV folds
#   color    Metric (Accuracy, Macro F1, Weighted F1)
#   shape    Nesting (circle = without, triangle = with)
#
# Points for the same x-position are dodged horizontally to
# prevent overlap across the three metrics.
# ---------------------------------------------------------------
ggplot(
  plot_data,
  aes(x = model, y = mean, color = metric, shape = nesting)
) +
  # Dodged points: one point per (metric × nesting) combination
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  
  # Error bars show ± 1 SD across the 10 CV folds,
  # matching the reporting format in the thesis (Section 3.4)
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.6),
    width     = 0.15,
    linewidth = 0.6
  ) +
  
  # One panel per scenario; 2 rows keeps the aspect ratio
  # consistent with the binary simulation figure (Figure 1)
  facet_wrap(~ scenario, nrow = 2) +
  
  labs(
    x     = NULL,
    y     = "Mean ± SD",
    color = "Metric",
    shape = "Model specification"
  ) +
  
  theme_minimal(base_size = 13) +
  
  # Metric legend first, then model specification legend
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  
  theme(
    legend.position   = "top",
    legend.box        = "vertical",
    legend.box.just   = "left",
    legend.justification = "left",
    legend.title      = element_text(size = 11, face = "bold"),
    legend.text       = element_text(size = 10),
    legend.spacing.y  = unit(2, "pt"),
    strip.text        = element_text(face = "bold"),
    panel.grid.minor  = element_blank(),
    # Visible panel border makes facet boundaries explicit
    panel.border      = element_rect(colour = "grey50", fill = NA, linewidth = 0.8)
  )


# ---------------------------------------------------------------
# SECTION 3: Save figure
#
# Written to images/ relative to the working directory.
# The directory must exist before running this script.
# 300 dpi at 8 × 6 inches matches typical journal submission
# requirements and the resolution used for the thesis figures.
# ---------------------------------------------------------------
ggsave(
  filename = "images/graph.png",
  width    = 8,
  height   = 6,
  dpi      = 300
)