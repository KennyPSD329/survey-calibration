# =============================================================================
# 01_make_table1_from_canonical_results.R
#
# Reproduces Table 1 (LaTeX label: tab:bias_vr_full) of the manuscript.
# Caption: "Simulation results for the calibration estimators of regression
# coefficients using different surrogate prediction models (n2 = 6000)."
#
# FAST ROUTE: reads pre-saved canonical Monte Carlo summaries; does NOT
# re-run the simulation.
#
# Run from:  simulation/
#   setwd("path/to/simulation")
#   source("code/01_make_table1_from_canonical_results.R")
#
# Outputs (written to results/):
#   table1_reproduced.csv
#   table1_reproduced.xlsx
#   table1_reproduced.tex
# =============================================================================

library(writexl)

DATA_DIR <- "results/canonical/table1"
OUT_DIR  <- "results"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Load canonical result files (n2 = 6000, gamma2[3] = 0.1, fit = 1)
# ---------------------------------------------------------------------------
beta_means <- read.csv(file.path(DATA_DIR, "beta_means_fit_1_n2_6000.csv"),
                       row.names = 1, check.names = TRUE)
beta_evars <- read.csv(file.path(DATA_DIR, "beta_empirical_variances_fit_1_n2_6000.csv"),
                       row.names = 1, check.names = TRUE)
ratio1     <- read.csv(file.path(DATA_DIR, "ratio_candidate1_fit_1_n2_6000.csv"),
                       row.names = 1, check.names = TRUE)
coverage   <- read.csv(file.path(DATA_DIR, "beta_coverage_fit_1_n2_6000.csv"),
                       row.names = 1, check.names = TRUE)
bias_raw   <- read.csv(file.path(DATA_DIR, "beta_bias_true_fit_1_n2_6000.csv"),
                       row.names = 1, check.names = TRUE)

# beta_bias_true was saved with numeric row indices (1-16); assign method names
all_methods <- c(
  "s1.wt.r",    "s1.wt.f",    "s1.wt.m1",    "s1.wt.m2",
  "s2.wt.r",    "s2.wt.f",    "s2.wt.m1",    "s2.wt.m2",
  "s1.cwt.i.r", "s1.cwt.i.f", "s1.cwt.i.m1", "s1.cwt.i.m2",
  "s1.cwts.r",  "s1.cwts.f",  "s1.cwts.m1",  "s1.cwts.m2"
)
rownames(bias_raw) <- all_methods
colnames(bias_raw) <- colnames(beta_means)   # V1-V5 -> beta.(Intercept), x1, x2, x3, x1.x3

# ---------------------------------------------------------------------------
# 2. True beta values used in the simulation setting (from the DGP).
#
# The true parameter vector is:
#   beta = (beta0, beta_x1, beta_x2, beta_x3, beta_x1:x3)
#        = (-3.5, 1.0, 0.5, 0.5, 0.2)
#
# Bias(%) for Table 1 is computed as:
#   Bias(%) = 100 * (mean_estimate - true_FP_glm) / beta_sim_true
# where:
#   - numerator  = bias_raw from canonical CSV = mean(beta_hat) - true_FP_glm
#                  (true_FP_glm = glm fit on the full N=1,000,000 population)
#   - denominator = simulation-setting beta below (not true_FP_glm), to match
#                   the manuscript Table 1 numbers exactly.
#
# Note: true_FP_glm is very close to but not identical to beta_sim_true
# (e.g., true_FP_x1 ~ 1.009, true_FP_x1x3 ~ 0.190) because cluster effects
# shift the finite-population regression slightly. Dividing by beta_sim_true
# reproduces the manuscript Bias(%) values exactly.
# ---------------------------------------------------------------------------
beta_sim_true <- c(
  "beta..Intercept." = -3.5,
  "beta.x1"          =  1.0,
  "beta.x2"          =  0.5,
  "beta.x3"          =  0.5,
  "beta.x1.x3"       =  0.2
)

# Verification: derive true_FP_glm from canonical files and show comparison
# (for cross-checking only; not used as denominator in Bias computation)
true_FP_glm <- as.numeric(beta_means["s1.wt.f", ]) - as.numeric(bias_raw["s1.wt.f", ])
names(true_FP_glm) <- colnames(beta_means)
cat("=== Verification: sim-setting betas vs. population-fitted betas ===\n")
cat("  (small differences due to cluster effects; Bias denominator uses sim-setting)\n")
comparison <- rbind(
  "beta_sim_true" = beta_sim_true[colnames(beta_means)],
  "true_FP_glm"   = round(true_FP_glm, 6)
)
print(comparison)
cat("\n")

# ---------------------------------------------------------------------------
# 3. Table 1 columns of interest (excluding intercept)
# ---------------------------------------------------------------------------
target_coefs <- c("beta.x1", "beta.x2", "beta.x3", "beta.x1.x3")
coef_labels  <- c("beta_x1", "beta_x2", "beta_x3", "beta_x1:x3")

# ---------------------------------------------------------------------------
# 4. Helper function: extract Table 1 quantities for a single method
# ---------------------------------------------------------------------------
get_quantities <- function(meth) {
  # Bias(%) = (mean_estimate - true_FP_glm) / beta_sim_true * 100
  # bias_raw stores the numerator (mean_estimate - true_FP_glm) in absolute form.
  # beta_sim_true (1.0, 0.5, 0.5, 0.2) is the denominator to match manuscript values.
  bias_pct <- as.numeric(bias_raw[meth, target_coefs]) /
                beta_sim_true[target_coefs] * 100

  vr     <- as.numeric(ratio1[meth, target_coefs])
  cr_pct <- as.numeric(coverage[meth, target_coefs]) * 100
  relvar <- as.numeric(beta_evars[meth, target_coefs]) /
              as.numeric(beta_evars["s1.wt.f", target_coefs])

  data.frame(
    Coef     = target_coefs,
    Bias_pct = round(bias_pct, 2),
    VR       = round(vr,       2),
    CR_pct   = round(cr_pct,   1),
    RelVar   = round(relvar,   2)
  )
}

# ---------------------------------------------------------------------------
# 5. Surrogate model block definitions
#
# NOTE on code vs. manuscript surrogate labeling:
#   Code m1 (a = 0.67) => rho(X3, X3*) ~ 0.80 => manuscript X3*(2)  (Block 3)
#   Code m2 (a = 1.55) => rho(X3, X3*) ~ 0.50 => manuscript X3*(1)  (Block 2)
# ---------------------------------------------------------------------------
blocks <- list(
  list(label    = "Reduced: Y~X1+X2",
       pooled   = "s1.cwt.i.r",
       external = "s1.cwts.r"),
  list(label    = "Surrogate X3*(1) rho~0.50 [code: m2]",
       pooled   = "s1.cwt.i.m2",
       external = "s1.cwts.m2"),
  list(label    = "Surrogate X3*(2) rho~0.80 [code: m1]",
       pooled   = "s1.cwt.i.m1",
       external = "s1.cwts.m1"),
  list(label    = "Full: Y~X1+X2+X3+X1:X3",
       pooled   = "s1.cwt.i.f",
       external = "s1.cwts.f")
)

# ---------------------------------------------------------------------------
# 6. Build Table 1
# ---------------------------------------------------------------------------
rows <- list()
for (blk in blocks) {
  p <- get_quantities(blk$pooled)
  e <- get_quantities(blk$external)
  for (ci in seq_along(target_coefs)) {
    cc <- target_coefs[ci]
    pr <- p[p$Coef == cc, ]
    er <- e[e$Coef == cc, ]
    rows[[length(rows) + 1]] <- data.frame(
      SurrogateModel    = blk$label,
      Coefficient       = coef_labels[ci],
      Pooled_Bias_pct   = if (nrow(pr) == 1) pr$Bias_pct else NA_real_,
      Pooled_VR         = if (nrow(pr) == 1) pr$VR       else NA_real_,
      Pooled_CR_pct     = if (nrow(pr) == 1) pr$CR_pct   else NA_real_,
      Pooled_RelVar     = if (nrow(pr) == 1) pr$RelVar   else NA_real_,
      External_Bias_pct = if (nrow(er) == 1) er$Bias_pct else NA_real_,
      External_VR       = if (nrow(er) == 1) er$VR       else NA_real_,
      External_CR_pct   = if (nrow(er) == 1) er$CR_pct   else NA_real_,
      External_RelVar   = if (nrow(er) == 1) er$RelVar   else NA_real_
    )
  }
}
table1 <- do.call(rbind, rows)

# ---------------------------------------------------------------------------
# 7. Save outputs
# ---------------------------------------------------------------------------
write.csv(table1,  file.path(OUT_DIR, "table1_reproduced.csv"),  row.names = FALSE)
write_xlsx(table1, file.path(OUT_DIR, "table1_reproduced.xlsx"))

# --- LaTeX output (simplified row format for copy-paste into manuscript) ---
tex_lines <- c(
  "% Table 1 reproduced -- generated by 01_make_table1_from_canonical_results.R",
  "% Columns: SurrogateModel | Coef | Pooled Bias(%) VR CR% RelVar | External Bias(%) VR CR% RelVar",
  "% Paste rows into the tabular environment in Main.tex."
)
prev_blk <- ""
for (k in seq_len(nrow(table1))) {
  r <- table1[k, ]
  blk_lbl <- as.character(r$SurrogateModel)
  if (blk_lbl != prev_blk) {
    tex_lines <- c(tex_lines,
                   paste0("\\midrule"),
                   paste0("\\multicolumn{10}{l}{\\textbf{Surrogate model}: ",
                          blk_lbl, "} \\\\"))
    prev_blk <- blk_lbl
  }
  coef_tex <- switch(as.character(r$Coefficient),
    "beta_x1"   = "$\\beta_{x_1}$",
    "beta_x2"   = "$\\beta_{x_2}$",
    "beta_x3"   = "$\\beta_{x_3}$",
    "beta_x1:x3"= "$\\beta_{x_1:x_3}$",
    as.character(r$Coefficient))
  tex_lines <- c(tex_lines, paste0(
    coef_tex,
    " & ", formatC(r$Pooled_Bias_pct,   format="f", digits=2),
    " & ", formatC(r$Pooled_VR,         format="f", digits=2),
    " & ", formatC(r$Pooled_CR_pct,     format="f", digits=1),
    " & ", formatC(r$Pooled_RelVar,     format="f", digits=2),
    " && ", formatC(r$External_Bias_pct, format="f", digits=2),
    " & ", formatC(r$External_VR,        format="f", digits=2),
    " & ", formatC(r$External_CR_pct,    format="f", digits=1),
    " & ", formatC(r$External_RelVar,    format="f", digits=2),
    " \\\\"))
}
writeLines(tex_lines, file.path(OUT_DIR, "table1_reproduced.tex"))

cat("\n=== Table 1 reproduced ===\n")
print(table1, row.names = FALSE)
cat("\nOutputs written to:\n")
cat("  ", file.path(OUT_DIR, "table1_reproduced.csv"),  "\n")
cat("  ", file.path(OUT_DIR, "table1_reproduced.xlsx"), "\n")
cat("  ", file.path(OUT_DIR, "table1_reproduced.tex"),  "\n")
