# =============================================================================
# 03_make_simulation_figure_from_canonical_results.R
#
# Reproduces Figure 3 (LaTeX label: fig:3) of the manuscript:
#   "Relative variance of the calibration estimator for (a) beta_x3 and
#    (b) beta_x1:x3 compared to the s1 estimator with increased
#    external sample size."
#
# FAST ROUTE: reads pre-saved canonical Monte Carlo summaries from three
# n2 folders and regenerates the figure. Does NOT re-run the simulation.
#
# Run from:  simulation/
#   setwd("path/to/simulation")
#   source("code/03_make_simulation_figure_from_canonical_results.R")
#
# Outputs (written to figures/):
#   Relative Variance.png
#   Relative Variance.pdf
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

OUT_DIR <- "figures"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Read empirical variance files from the three canonical n2 folders
#    (all 7 CSVs per folder are read to match the original script; only
#     beta_empirical_variances is used for the relative variance figure)
# ---------------------------------------------------------------------------
read_n2_data <- function(n2_val, folder) {
  i <- 1
  fnames <- c(
    paste0("beta_means_fit_",                i, "_n2_", n2_val, ".csv"),
    paste0("beta_empirical_variances_fit_",  i, "_n2_", n2_val, ".csv"),
    paste0("beta_variance_candidate1_fit_",  i, "_n2_", n2_val, ".csv"),
    paste0("ratio_candidate1_fit_",          i, "_n2_", n2_val, ".csv"),
    paste0("beta_bias_fit_",                 i, "_n2_", n2_val, ".csv"),
    paste0("beta_bias_true_fit_",            i, "_n2_", n2_val, ".csv"),
    paste0("beta_coverage_fit_",             i, "_n2_", n2_val, ".csv")
  )
  vnames <- c("beta.means", "beta.evars", "beta.var1", "ratio1",
              "beta.bias",  "beta.bias.true", "beta.coverage")
  ml <- lapply(fnames, function(f)
    read.csv(file.path(folder, f), row.names = 1, check.names = TRUE))
  names(ml) <- vnames
  ml
}

data_2000  <- read_n2_data(2000,  "results/canonical/figure/n2_2000")
data_6000  <- read_n2_data(6000,  "results/canonical/figure/n2_6000")
data_20000 <- read_n2_data(20000, "results/canonical/figure/n2_20000")

# ---------------------------------------------------------------------------
# 2. Compute relative variances for beta_x3 and beta_x1:x3
#    Baseline (denominator): row 2 = s1.wt.f (s1-only full-model estimator)
# ---------------------------------------------------------------------------
add_relvars <- function(evars, n2_val) {
  evars$sample_size      <- n2_val
  evars$method           <- rownames(evars)
  # R reads "beta.x1:x3" column as "beta.x1.x3" due to check.names = TRUE
  evars$relative_var.x3    <- evars$beta.x3    / evars$beta.x3[2]
  evars$relative_var.x1.x3 <- evars$beta.x1.x3 / evars$beta.x1.x3[2]
  evars
}

ev2    <- add_relvars(data_2000$beta.evars,  2000)
ev6    <- add_relvars(data_6000$beta.evars,  6000)
ev20   <- add_relvars(data_20000$beta.evars, 20000)

# ---------------------------------------------------------------------------
# 3. Subset rows 9-16 (calibration methods) and assign manuscript labels
#    Row order (1-indexed):
#      9  = s1.cwt.i.r  (Pooled, Reduced)
#      10 = s1.cwt.i.f  (Pooled, Full)
#      11 = s1.cwt.i.m1 (Pooled, Surrogate rho~0.80)   [code: m1, a=0.67]
#      12 = s1.cwt.i.m2 (Pooled, Surrogate rho~0.50)   [code: m2, a=1.55]
#      13 = s1.cwts.r   (External, Reduced)
#      14 = s1.cwts.f   (External, Full)
#      15 = s1.cwts.m1  (External, Surrogate rho~0.80)
#      16 = s1.cwts.m2  (External, Surrogate rho~0.50)
#
# estimator label: "Individual-Level" = Pooled-Sample
#                  "Model-Assisted"   = External-Sample
# model label:     Surrogate Model 1 = m1 (rho~0.80)
#                  Surrogate Model 2 = m2 (rho~0.50)
# ---------------------------------------------------------------------------
method_labels <- c(
  "Individual-Level,Reduced Model",
  "Individual-Level,Full Model",
  "Individual-Level,Surrogate Model 1",
  "Individual-Level,Surrogate Model 2",
  "Model-Assisted,Reduced Model",
  "Model-Assisted,Full Model",
  "Model-Assisted,Surrogate Model 1",
  "Model-Assisted,Surrogate Model 2"
)

prepare_revars <- function(ev) {
  rv <- ev[9:16, c("sample_size","method","relative_var.x3","relative_var.x1.x3")]
  rv$method <- method_labels
  rv
}

revars_combined <- rbind(
  prepare_revars(ev2),
  prepare_revars(ev6),
  prepare_revars(ev20)
)

# ---------------------------------------------------------------------------
# 4. Pivot to long format
# ---------------------------------------------------------------------------
revars_long <- revars_combined %>%
  pivot_longer(
    cols      = starts_with("relative_var"),
    names_to  = "coefficient",
    values_to = "relative_variance"
  ) %>%
  separate(method, into = c("estimator", "model"), sep = ",")

# ---------------------------------------------------------------------------
# 5. Produce manuscript Figure 3 (Relative Variance vs External Sample Size)
# ---------------------------------------------------------------------------
fig3 <- ggplot(revars_long,
               aes(sample_size, relative_variance,
                   color = model, linetype = estimator)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~ coefficient, scales = "fixed") +

  scale_x_continuous(
    breaks      = sort(unique(revars_long$sample_size)),
    minor_breaks = NULL
  ) +

  labs(
    x        = expression("External Sample Size (" * n[2] * ")"),
    y        = "Relative Variance",
    color    = "Surrogate Model",
    linetype = "Estimator Type"
  ) +

  scale_color_manual(
    values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),
    labels = c(
      "Full Model",
      "Reduced Model",
      expression("Surrogate model with " * X[3]^"*" * " (" * rho == 0.8 * ")"),
      expression("Surrogate model with " * X[3]^"*" * " (" * rho == 0.5 * ")")
    )
  ) +

  scale_linetype_manual(
    values = c("solid", "dashed"),
    labels = c("Pooled-Sample", "External-Sample")
  ) +

  theme_minimal() +
  theme(
    plot.title       = element_blank(),
    strip.text       = element_blank(),
    strip.background = element_blank(),
    panel.spacing    = unit(1, "lines"),
    panel.grid.minor.x = element_blank(),
    axis.title  = element_text(size = 18),
    axis.text   = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )

# ---------------------------------------------------------------------------
# 6. Save figure
# ---------------------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Relative Variance.png"),
       plot = fig3, width = 10, height = 5, dpi = 300)
ggsave(file.path(OUT_DIR, "Relative Variance.pdf"),
       plot = fig3, width = 10, height = 5)

cat("Figure 3 saved to:\n")
cat("  ", file.path(OUT_DIR, "Relative Variance.png"), "\n")
cat("  ", file.path(OUT_DIR, "Relative Variance.pdf"), "\n")
