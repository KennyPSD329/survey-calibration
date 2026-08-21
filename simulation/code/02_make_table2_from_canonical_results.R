# =============================================================================
# 02_make_table2_from_canonical_results.R
#
# Reproduces Table 2 (LaTeX label: tab:simulation_theta) of the manuscript.
# Caption: "Comparison of calibration and plug-in estimators (n2 = 6000)."
#
# FAST ROUTE: reads pre-saved canonical Monte Carlo summaries; does NOT
# re-run the simulation.
#
# Run from:  simulation/
#   setwd("path/to/simulation")
#   source("code/02_make_table2_from_canonical_results.R")
#
# Outputs (written to results/):
#   table2_reproduced.csv
#   table2_reproduced.xlsx
#   table2_reproduced.tex
#
# ---------------------------------------------------------------------------
# Excel file structure (beta_means and beta_evars xlsx):
#
# Both xlsx files have 9 columns after auto-rename by readxl:
#   Col 1  (...1)           : method labels (rownames)
#   Cols 2-5 (*.x1...2 etc.): fit=1 results (x1, x2, x3, x1:x3)
#   Cols 6-9 (*.x1...6 etc.): pre-computed summary quantities
#
# For beta_means xlsx:
#   Cols 6-9 = (mean - true_FP_glm) / true_FP_glm  [relative bias as fraction]
#   -> Bias(%) = col_value * 100
#
# For beta_evars xlsx:
#   Cols 6-9 = RV = empirical_var(IF_estimator) / empirical_var(calib_estimator)
#              (NA for calibration rows like *.beta.i and *.beta.s)
#   -> RV = col_value directly (round to 2 decimal places)
#
# Confirmed by spot-checks against manuscript Table 2 values.
# =============================================================================

library(readxl)
library(writexl)

T2_DIR  <- "results/canonical/table2"
OUT_DIR <- "results"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Read xlsx files and clean up structure
# ---------------------------------------------------------------------------
read_t2_xlsx <- function(fname) {
  df <- suppressMessages(read_excel(file.path(T2_DIR, fname)))
  df <- as.data.frame(df)
  # Drop trailing NA rows (xlsx may have 2 extra blank rows)
  df <- df[!is.na(df[[1]]), ]
  # Column layout: col1 = method label; cols 2-5 = fit-1 quantities;
  #                cols 6-9 = pre-computed summary (bias or RV)
  rownames(df) <- df[[1]]
  df <- df[, -1]   # drop the method-name column; rownames now carry it
  # Give clean names to make indexing unambiguous
  colnames(df) <- c("x1_fit1", "x2_fit1", "x3_fit1", "xi_fit1",
                    "x1_sum",  "x2_sum",  "x3_sum",  "xi_sum")
  df
}

means_df <- read_t2_xlsx("beta_means_fit_1_n2_6000.xlsx")
evars_df <- read_t2_xlsx("beta_empirical_variances_fit_1_n2_6000.xlsx")

cat("=== Methods in Table 2 xlsx ===\n"); print(rownames(means_df))

# ---------------------------------------------------------------------------
# 2. Helper: Bias(%) and RV for a single IF estimator row and coefficient
#
#   Bias(%) = means_df[meth, "x*_sum"] * 100
#             because means_df cols 6-9 store (mean - true_FP_glm)/true_FP_glm
#
#   RV      = evars_df[meth, "x*_sum"]  (pre-computed ratio; NA for calib rows)
# ---------------------------------------------------------------------------
coef_sums <- c(
  "beta.x1"    = "x1_sum",
  "beta.x2"    = "x2_sum",
  "beta.x3"    = "x3_sum",
  "beta.x1.x3" = "xi_sum"
)

get_bias_pct <- function(meth, coef_key) {
  if (!meth %in% rownames(means_df)) return(NA_real_)
  col <- coef_sums[coef_key]
  v   <- as.numeric(means_df[meth, col])
  round(v * 100, 2)
}

get_rv <- function(meth, coef_key) {
  if (!meth %in% rownames(evars_df)) return(NA_real_)
  col <- coef_sums[coef_key]
  v   <- as.numeric(evars_df[meth, col])
  round(v, 2)
}

# ---------------------------------------------------------------------------
# 3. Surrogate block definitions (4 blocks)
#
# NOTE on m1/m2 to manuscript surrogate mapping:
#   Code m1 (a = 0.67) => rho(X3, X3*) ~ 0.80 => manuscript X3*(2)  (Block 3)
#   Code m2 (a = 1.55) => rho(X3, X3*) ~ 0.50 => manuscript X3*(1)  (Block 2)
# ---------------------------------------------------------------------------
blocks <- list(
  list(label = "Reduced: Y~X1+X2",                     sc = "r"),
  list(label = "Surrogate X3*(1) rho~0.50 [code: m2]", sc = "m2"),
  list(label = "Surrogate X3*(2) rho~0.80 [code: m1]", sc = "m1"),
  list(label = "Full: Y~X1+X2+X3+X1:X3",               sc = "f")
)

# Table 2 IF estimator methods per scenario:
#   Pooled beta_hat_p(Theta)     : *.beta.if.S.pop
#   Pooled beta_hat_p(Theta_hat) : *.beta.if.S
#   External beta_hat_e(Theta)   : *.beta.if.s2.pop
#   External beta_hat_e(Theta_hat): *.beta.if.s2

target_coefs  <- c("beta.x1", "beta.x2", "beta.x3", "beta.x1.x3")
coef_labels   <- c("beta_x1", "beta_x2", "beta_x3", "beta_x1:x3")

# ---------------------------------------------------------------------------
# 4. Build Table 2
# ---------------------------------------------------------------------------
rows2 <- list()
for (blk in blocks) {
  sc <- blk$sc
  m_p_Th     <- paste0(sc, ".beta.if.S.pop")    # Pooled, beta_hat_p(Theta)
  m_p_Thhat  <- paste0(sc, ".beta.if.S")         # Pooled, beta_hat_p(Theta_hat)
  m_e_Th     <- paste0(sc, ".beta.if.s2.pop")   # External, beta_hat_e(Theta)
  m_e_Thhat  <- paste0(sc, ".beta.if.s2")        # External, beta_hat_e(Theta_hat)

  for (ci in seq_along(target_coefs)) {
    ck <- target_coefs[ci]
    rows2[[length(rows2) + 1]] <- data.frame(
      SurrogateModel                = blk$label,
      Coefficient                   = coef_labels[ci],
      Pooled_pTheta_Bias_pct        = get_bias_pct(m_p_Th,    ck),
      Pooled_pThetahat_Bias_pct     = get_bias_pct(m_p_Thhat, ck),
      Pooled_pTheta_RV              = get_rv(m_p_Th,    ck),
      Pooled_pThetahat_RV           = get_rv(m_p_Thhat, ck),
      External_eTheta_Bias_pct      = get_bias_pct(m_e_Th,    ck),
      External_eThetahat_Bias_pct   = get_bias_pct(m_e_Thhat, ck),
      External_eTheta_RV            = get_rv(m_e_Th,    ck),
      External_eThetahat_RV         = get_rv(m_e_Thhat, ck)
    )
  }
}
table2 <- do.call(rbind, rows2)

# ---------------------------------------------------------------------------
# 5. Save CSV and xlsx
# ---------------------------------------------------------------------------
write.csv(table2,  file.path(OUT_DIR, "table2_reproduced.csv"),  row.names = FALSE)
write_xlsx(table2, file.path(OUT_DIR, "table2_reproduced.xlsx"))

# ---------------------------------------------------------------------------
# 6. LaTeX output (simplified row format for copy-paste)
# ---------------------------------------------------------------------------
tex_lines <- c(
  "% Table 2 reproduced -- generated by 02_make_table2_from_canonical_results.R",
  "% Columns: Coef | Pooled Bias%(Theta) Bias%(Theta_hat) RV(Theta) RV(Theta_hat)",
  "%               | External Bias%(Theta) Bias%(Theta_hat) RV(Theta) RV(Theta_hat)"
)
prev_blk <- ""
for (k in seq_len(nrow(table2))) {
  r <- table2[k, ]
  blk_lbl <- as.character(r$SurrogateModel)
  if (blk_lbl != prev_blk) {
    tex_lines <- c(tex_lines,
                   "\\midrule",
                   paste0("\\multicolumn{12}{l}{\\textbf{Surrogate model}: ",
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
    " & ", formatC(r$Pooled_pTheta_Bias_pct,      format="f", digits=2),
    " & ", formatC(r$Pooled_pThetahat_Bias_pct,   format="f", digits=2),
    " && ", formatC(r$Pooled_pTheta_RV,            format="f", digits=2),
    " & ", formatC(r$Pooled_pThetahat_RV,          format="f", digits=2),
    " && ", formatC(r$External_eTheta_Bias_pct,    format="f", digits=2),
    " & ", formatC(r$External_eThetahat_Bias_pct,  format="f", digits=2),
    " && ", formatC(r$External_eTheta_RV,          format="f", digits=2),
    " & ", formatC(r$External_eThetahat_RV,        format="f", digits=2),
    " \\\\"))
}
writeLines(tex_lines, file.path(OUT_DIR, "table2_reproduced.tex"))

cat("\n=== Table 2 reproduced ===\n")
print(table2, row.names = FALSE)
cat("\nOutputs written to:\n")
cat("  ", file.path(OUT_DIR, "table2_reproduced.csv"),  "\n")
cat("  ", file.path(OUT_DIR, "table2_reproduced.xlsx"), "\n")
cat("  ", file.path(OUT_DIR, "table2_reproduced.tex"),  "\n")
