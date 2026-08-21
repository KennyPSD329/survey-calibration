# ==============================================================================
# 02_run_data_application_models.R
#
# Sources helper functions, loads processed data from data/processed/, runs
# all calibration models, and writes results to results/.
#
# Run from the application/ root directory after Script 01.
# ==============================================================================

rm(list = ls())
library(survey)
library(dplyr)
library(broom)
library(writexl)

source("code/subfunctions20250417_V1_3.R")
source("code/beta_est_var_Biometrics_2.R")

if (!file.exists("data/processed/nhanes_analysis.rds")) {
  stop("Run 01_prepare_data_application_data.R first.")
}

# ==============================================================================
# 1. Load processed data
# ==============================================================================
cat("Loading processed data...\n")
nhanes_all <- readRDS("data/processed/nhanes_analysis.rds")
nhis       <- readRDS("data/processed/nhis_analysis.rds")

# ==============================================================================
# 2. Set up NHANES internal / external / combined datasets
# ==============================================================================

# Full 3-cycle dataset (2001-02 external, 2003-04 internal, 2005-06 external)
nhanes_nhanes          <- nhanes_all
nhanes_nhanes$internal <- as.numeric(nhanes_nhanes$source == "internal")

# External cycles only (2001-02 and 2005-06)
nhanes_ext <- nhanes_all[nhanes_all$source == "external", ]

# Internal cycle only (2003-04) — also used as internal sample in NHIS analysis
nhanes <- nhanes_all[nhanes_all$source == "internal", ]

# ==============================================================================
# 3. Survey designs for NHANES-only analysis
# ==============================================================================
# WTINT2YR  : 2-year weight (internal 03-04 only)
# WTINT2YR3 : WTINT2YR/2 for 2-cycle external pooled design
# WTINT2YR4 : WTINT2YR/3 for 3-cycle combined design
ds.nhanes         <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA,
                               weight = ~WTINT2YR,  nest = TRUE, data = nhanes)
ds.nhanes_ext     <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA,
                               weight = ~WTINT2YR3, nest = TRUE, data = nhanes_ext)
ds.combined_nhanes <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA,
                                weight = ~WTINT2YR4, nest = TRUE, data = nhanes_nhanes)

# ==============================================================================
# 4. Set up NHIS combined dataset and survey designs
# ==============================================================================

# For NHIS analysis, nhanes (internal 03-04) gets source = 'nhanes'
nhanes_for_nhis        <- nhanes
nhanes_for_nhis$source <- "nhanes"
nhanes_for_nhis        <- nhanes_for_nhis[, !(names(nhanes_for_nhis) %in%
                                               c("WTINT2YR4", "WTINT2YR3", "cycle"))]

nhis_nhanes <- rbind(nhanes_for_nhis, nhis)
n_combined  <- nrow(nhis_nhanes)

# Proportional weight adjustment for combined NHANES+NHIS design
nhis_nhanes$WTINT2YR4 <- c(
  nhanes_for_nhis$WTINT2YR * nrow(nhanes_for_nhis) / n_combined,
  nhis$WTINT2YR            * nrow(nhis)             / n_combined
)
nhis_nhanes$internal <- as.numeric(nhis_nhanes$source == "nhanes")

ds.nhis          <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA,
                              weight = ~WTINT2YR,  nest = TRUE, data = nhis)
ds.combined_nhis <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA,
                              weight = ~WTINT2YR4, nest = TRUE, data = nhis_nhanes)

# ==============================================================================
# 5. Model formulas
# ==============================================================================
mort.fit   <- as.formula(
  "mortality~RIDAGEYR + RIAGENDR + as.factor(race4) + phys.all +
   as.factor(smk3) + as.factor(hvy.drk) + DXDTOFAT"
)
mort.fit.i <- as.formula(
  "mortality~RIDAGEYR + RIAGENDR + as.factor(race4) + phys.all +
   as.factor(smk3) + as.factor(hvy.drk) + slfBMI"
)

# ==============================================================================
# 6. Preliminary model comparison across cycles (8 covariate models)
# ==============================================================================
cat("Fitting preliminary multi-model comparison...\n")

m_list <- list(
  m1 = mortality ~ slfBMI,
  m2 = mortality ~ slfBMI + RIDAGEYR,
  m3 = mortality ~ slfBMI + RIDAGEYR + height,
  m4 = mortality ~ slfBMI + RIDAGEYR + height + RIAGENDR,
  m5 = mortality ~ slfBMI + RIDAGEYR + height + RIAGENDR + as.factor(race4),
  m6 = mortality ~ slfBMI + RIDAGEYR + height + RIAGENDR + as.factor(race4) + phys.all,
  m7 = mortality ~ slfBMI + RIDAGEYR + height + RIAGENDR + as.factor(race4) + phys.all + as.factor(smk3),
  m8 = mortality ~ slfBMI + RIDAGEYR + height + RIAGENDR + as.factor(race4) + phys.all + as.factor(smk3) + as.factor(hvy.drk)
)

fit_all <- function(m) svyglm(m, ds.nhanes, family = "binomial")
fits_nhanes          <- lapply(m_list, fit_all)
fits_nhanes_ext      <- lapply(m_list, function(m) svyglm(m, ds.nhanes_ext,      family = "binomial"))
fits_nhis            <- lapply(m_list, function(m) svyglm(m, ds.nhis,            family = "binomial"))
fits_combined_nhanes <- lapply(m_list, function(m) svyglm(m, ds.combined_nhanes, family = "binomial"))
fits_combined_nhis   <- lapply(m_list, function(m) svyglm(m, ds.combined_nhis,   family = "binomial"))

get_coefs  <- function(fl) lapply(fl, coef)
get_pvals  <- function(fl) lapply(fl, function(f) {
  tb <- as.data.frame(summary(f)$coefficients)
  pv <- tb[, grep("^Pr", colnames(tb))[1], drop = TRUE]
  names(pv) <- rownames(tb); pv
})

coefs_nhanes          <- get_coefs(fits_nhanes)
coefs_nhanes_ext      <- get_coefs(fits_nhanes_ext)
coefs_nhis            <- get_coefs(fits_nhis)
coefs_combined_nhanes <- get_coefs(fits_combined_nhanes)
coefs_combined_nhis   <- get_coefs(fits_combined_nhis)
pvals_nhanes          <- get_pvals(fits_nhanes)
pvals_nhanes_ext      <- get_pvals(fits_nhanes_ext)
pvals_nhis            <- get_pvals(fits_nhis)
pvals_combined_nhanes <- get_pvals(fits_combined_nhanes)
pvals_combined_nhis   <- get_pvals(fits_combined_nhis)

alpha <- 0.05
deviations <- mapply(function(b1, b2, b3, b4, b5, p1, p2, p3, p4, p5) {
  all.names <- Reduce(union, list(names(b1), names(b2), names(b3), names(b4), names(b5)))
  b1 <- b1[all.names]; b2 <- b2[all.names]; b3 <- b3[all.names]
  b4 <- b4[all.names]; b5 <- b5[all.names]
  p1 <- p1[all.names]; p2 <- p2[all.names]; p3 <- p3[all.names]
  p4 <- p4[all.names]; p5 <- p5[all.names]
  data.frame(
    NHANES           = b1,
    NHANES_EXT       = b2,
    NHIS             = b3,
    NHANES_COMBINE   = b4,
    NHIS_COMBINE     = b5,
    Dev_NHANES_EXT   = b2 - b1,
    Dev_NHIS         = b3 - b1,
    Dev_NHANES_COMBINE = b4 - b1,
    Dev_NHIS_COMBINE = b5 - b1,
    p_NHANES         = p1,
    p_NHANES_EXT     = p2,
    p_NHIS           = p3,
    p_NHANES_COMBINED = p4,
    p_NHIS_COMBINED  = p5,
    Sig_NHANES       = ifelse(p1 < alpha, "Yes", "No"),
    Sig_NHANES_EXT   = ifelse(p2 < alpha, "Yes", "No"),
    Sig_NHIS         = ifelse(p3 < alpha, "Yes", "No"),
    Sig_NHANES_COMBINED = ifelse(p4 < alpha, "Yes", "No"),
    Sig_NHIS_COMBINED = ifelse(p5 < alpha, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}, coefs_nhanes, coefs_nhanes_ext, coefs_nhis, coefs_combined_nhanes, coefs_combined_nhis,
   pvals_nhanes, pvals_nhanes_ext, pvals_nhis, pvals_combined_nhanes, pvals_combined_nhis,
   SIMPLIFY = FALSE)

result_table <- bind_rows(
  lapply(seq_along(m_list), function(i) deviations[[i]] %>% mutate(Model = paste0("m", i)))
) %>% mutate(across(where(is.numeric), ~round(., 3)))

print(result_table)
write.csv(result_table, "results/Model_Bias_Final.csv", row.names = TRUE)

rm(coefs_nhanes, coefs_nhanes_ext, coefs_nhis, coefs_combined_nhanes, coefs_combined_nhis,
   pvals_nhanes, pvals_nhanes_ext, pvals_nhis, pvals_combined_nhanes, pvals_combined_nhis,
   fits_nhanes, fits_nhanes_ext, fits_nhis, fits_combined_nhanes, fits_combined_nhis, deviations)

# ==============================================================================
# 7. NHANES-only calibration analysis
# ==============================================================================
cat("Running NHANES-only calibration...\n")

# Row indices: external rows come first in nhanes_nhanes (rbind order: ext01, ext05, int03)
matched_idx     <- (nrow(nhanes_ext) + 1):nrow(nhanes_nhanes)
matched_idx_ext <- 1:nrow(nhanes_ext)

# Gold standard: full 3-cycle combined
mort.FAT.GOLD <- svyglm(mort.fit, family = "binomial", design = ds.combined_nhanes)

# NHANES external only
mort.FAT.NHANES_EXT <- svyglm(mort.fit, family = "binomial", design = ds.nhanes_ext)

# Internal-only (naive)
s1.fit      <- svyglm(mort.fit, ds.nhanes, family = "binomial")
beta.1      <- coef(s1.fit)
beta.var.1  <- diag(vcov(s1.fit))

# Pooled calibration via beta.est.var (from beta_est_var_Biometrics_2.R)
calib.out.nhanes <- beta.est.var(
  ds        = ds.combined_nhanes,
  fit       = mort.fit,
  fit.i     = mort.fit.i,
  sub       = "internal",
  sub.wt    = "WTINT2YR",
  structure = "combine"
)

# Combined calibration (influence-function GREG on internal sample)
S.fit        <- svyglm(mort.fit.i, ds.combined_nhanes, family = "binomial")
X.full       <- model.matrix(S.fit)
fitted.full  <- fitted(S.fit)
X.sub        <- X.full[matched_idx, ]
p.sub        <- fitted.full[matched_idx]
y.sub        <- ds.nhanes$variables[[all.vars(mort.fit.i)[1]]]
ui_nhanes    <- (y.sub - p.sub) * X.sub

ds.nhanes$variables$cwt.i <- weights(ds.nhanes) *
  greg_f(wt0 = weights(ds.nhanes), v.mtx0 = ui_nhanes,
         VS.hat = rep(0, ncol(X.full)))$f
ds.nhanes$variables$cwt.i[ds.nhanes$variables$cwt.i < 0] <- 0
ds.nhanes.i <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                         weights = ~cwt.i, nest = TRUE,
                         data = ds.nhanes$variables)
beta.i.nhanes <- coef(svyglm(mort.fit, ds.nhanes.i, family = "binomial"))

beta.var.i.fit.nhanes <- calibrate_beta_variance(
  ds_s         = ds.nhanes.i,
  ds_S         = ds.combined_nhanes,
  formula      = mort.fit.i,
  formula_i    = mort.fit,
  u_star       = ui_nhanes,
  wt0          = weights(ds.nhanes),
  structure    = "combined",
  ds2          = ds.nhanes_ext,
  n1           = nrow(ds.nhanes$variables),
  n2           = nrow(ds.nhanes_ext$variables),
  ds1          = ds.nhanes,
  matched_idx  = matched_idx_ext
)
beta.var.i.nhanes <- diag(beta.var.i.fit.nhanes$var_beta_calib)

# Model-based calibration
s2.fit.nhanes <- svyglm(mort.fit.i, ds.nhanes_ext, family = "binomial")
x.s1 <- model.matrix(lm(mort.fit.i, data = ds.nhanes$variables))
odds.hat <- exp(x.s1 %*% coef(s2.fit.nhanes))
p.hat    <- odds.hat / (1 + odds.hat)
y        <- ds.nhanes$variables[[all.vars(mort.fit)[1]]]
ui_mb_nhanes <- c(y - p.hat) * x.s1

ds.nhanes$variables$cwt.s <- weights(ds.nhanes) *
  greg_f(wt0 = weights(ds.nhanes), v.mtx0 = ui_mb_nhanes,
         VS.hat = rep(0, ncol(x.s1)))$f
ds.nhanes$variables$cwt.s[ds.nhanes$variables$cwt.s < 0] <- 0
ds.nhanes.s <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                         weights = ~cwt.s, nest = TRUE,
                         data = ds.nhanes$variables)
beta.s.nhanes <- coef(svyglm(mort.fit, ds.nhanes.s, family = "binomial"))

beta.var.s1.fit.nhanes <- calibrate_beta_variance(
  ds_s         = ds.nhanes.s,
  ds_S         = ds.nhanes,
  formula      = mort.fit.i,
  formula_i    = mort.fit,
  u_star       = ui_mb_nhanes,
  wt0          = weights(ds.nhanes),
  structure    = "model-based",
  var_beta_star = vcov(s2.fit.nhanes),
  data         = ds.nhanes$variables,
  beta         = coef(s2.fit.nhanes)
)
beta.var.s1.nhanes <- diag(beta.var.s1.fit.nhanes$var_beta_calib)

# ==============================================================================
# 8. NHIS calibration analysis
# ==============================================================================
cat("Running NHIS calibration...\n")

matched_idx_nhis     <- 1:nrow(nhanes_for_nhis)
matched_idx_ext_nhis <- (nrow(nhanes_for_nhis) + 1):nrow(nhis_nhanes)

# Refit s1.fit with the updated ds.nhanes (may have cwt.i/cwt.s in variables now)
s1.fit <- svyglm(mort.fit, ds.nhanes, family = "binomial")

# Pooled calibration (NHIS)
calib.out.nhis <- beta.est.var(
  ds        = ds.combined_nhis,
  fit       = mort.fit,
  fit.i     = mort.fit.i,
  sub       = "internal",
  sub.wt    = "WTINT2YR",
  structure = "combine"
)

# Combined calibration (NHIS)
S.fit.nhis    <- svyglm(mort.fit.i, ds.combined_nhis, family = "binomial")
X.full.nhis   <- model.matrix(S.fit.nhis)
fitted.nhis   <- fitted(S.fit.nhis)
X.sub.nhis    <- X.full.nhis[matched_idx_nhis, ]
p.sub.nhis    <- fitted.nhis[matched_idx_nhis]
y.sub.nhis    <- ds.nhanes$variables[[all.vars(mort.fit.i)[1]]]
ui_nhis       <- (y.sub.nhis - p.sub.nhis) * X.sub.nhis

ds.nhanes$variables$cwt.i <- weights(ds.nhanes) *
  greg_f(wt0 = weights(ds.nhanes), v.mtx0 = ui_nhis,
         VS.hat = rep(0, ncol(X.full.nhis)))$f
ds.nhanes$variables$cwt.i[ds.nhanes$variables$cwt.i < 0] <- 0
ds.nhanes.i.nhis <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                               weights = ~cwt.i, nest = TRUE,
                               data = ds.nhanes$variables)
beta.i.nhis <- coef(svyglm(mort.fit, ds.nhanes.i.nhis, family = "binomial"))

beta.var.i.fit.nhis <- calibrate_beta_variance(
  ds_s         = ds.nhanes.i.nhis,
  ds_S         = ds.combined_nhis,
  formula      = mort.fit.i,
  formula_i    = mort.fit,
  u_star       = ui_nhis,
  wt0          = weights(ds.nhanes),
  structure    = "combined",
  ds2          = ds.nhis,
  n1           = nrow(ds.nhanes$variables),
  n2           = nrow(ds.nhis$variables),
  ds1          = ds.nhanes,
  matched_idx  = matched_idx_ext_nhis
)
beta.var.i.nhis <- diag(beta.var.i.fit.nhis$var_beta_calib)

# Model-based calibration (NHIS)
s2.fit.nhis <- svyglm(mort.fit.i, ds.nhis, family = "binomial")
x.s1.nhis   <- model.matrix(lm(mort.fit.i, data = ds.nhanes$variables))
odds.hat.nhis <- exp(x.s1.nhis %*% coef(s2.fit.nhis))
p.hat.nhis    <- odds.hat.nhis / (1 + odds.hat.nhis)
y.nhis        <- ds.nhanes$variables[[all.vars(mort.fit)[1]]]
ui_mb_nhis    <- c(y.nhis - p.hat.nhis) * x.s1.nhis

ds.nhanes$variables$cwt.s <- weights(ds.nhanes) *
  greg_f(wt0 = weights(ds.nhanes), v.mtx0 = ui_mb_nhis,
         VS.hat = rep(0, ncol(x.s1.nhis)))$f
ds.nhanes$variables$cwt.s[ds.nhanes$variables$cwt.s < 0] <- 0
ds.nhanes.s.nhis <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                               weights = ~cwt.s, nest = TRUE,
                               data = ds.nhanes$variables)
beta.s.nhis <- coef(svyglm(mort.fit, ds.nhanes.s.nhis, family = "binomial"))

beta.var.s1.fit.nhis <- calibrate_beta_variance(
  ds_s          = ds.nhanes.s.nhis,
  ds_S          = ds.nhanes,
  formula       = mort.fit.i,
  formula_i     = mort.fit,
  u_star        = ui_mb_nhis,
  wt0           = weights(ds.nhanes),
  structure     = "model-based",
  var_beta_star = vcov(s2.fit.nhis),
  data          = ds.nhanes$variables,
  beta          = coef(s2.fit.nhis)
)
beta.var.s1.nhis <- diag(beta.var.s1.fit.nhis$var_beta_calib)

# ==============================================================================
# 9. Assemble results matrix (table4)
# ==============================================================================
cat("Assembling results...\n")

results <- cbind(
  coef(mort.FAT.GOLD),              diag(vcov(mort.FAT.GOLD)),
  coef(mort.FAT.NHANES_EXT),        diag(vcov(mort.FAT.NHANES_EXT)),
  beta.1,                            beta.var.1,
  calib.out.nhanes$beta.est,        calib.out.nhanes$var.all1[seq_along(calib.out.nhanes$beta.est)],
  beta.i.nhanes,                     beta.var.i.nhanes,
  beta.s.nhanes,                     beta.var.s1.nhanes,
  calib.out.nhis$beta.est,          calib.out.nhis$var.all1[seq_along(calib.out.nhis$beta.est)],
  beta.i.nhis,                       beta.var.i.nhis,
  beta.s.nhis,                       beta.var.s1.nhis
)

table4 <- results
colnames(table4) <- c(
  "Gold.beta",           "Gold.var",
  "NHANES.ext.beta",     "NHANES.ext.beta.var",
  "NHANES05.beta",       "NHANES05.beta.beta.var",
  "calib.NHANES.beta",   "calib.NHANES.var",
  "Combined.NHANES.beta","Combined.NHANES.beta.var",
  "M-Based.NHANES.beta", "M-Based.NHANES.beta.var",
  "calib.NHIS.beta",     "calib.NHIS.var",
  "Combined.NHIS.beta",  "Combined.NHIS.beta.var",
  "M-Based.NHIS.beta",   "M-Based.NHIS.beta.var"
)

# ==============================================================================
# 10. Build CI table (tidy format, display-scaled)
# ==============================================================================
z <- 1.96

# Mapping: estimator name -> (beta column, var column, display label)
ci_map <- list(
  list(beta = "Gold.beta",            var = "Gold.var",                label = "Benchmark: NHANES 01-06"),
  list(beta = "NHANES05.beta",        var = "NHANES05.beta.beta.var",  label = "Internal: NHANES 03-04"),
  list(beta = "calib.NHANES.beta",    var = "calib.NHANES.var",        label = "NHANES pooled calibration"),
  list(beta = "M-Based.NHANES.beta",  var = "M-Based.NHANES.beta.var", label = "NHANES external calibration"),
  list(beta = "calib.NHIS.beta",      var = "calib.NHIS.var",          label = "NHIS pooled calibration"),
  list(beta = "M-Based.NHIS.beta",    var = "M-Based.NHIS.beta.var",   label = "NHIS external calibration")
)

# Variable names as they appear in rownames(table4), plus display label and scale factor
var_map <- list(
  list(coef = "RIDAGEYR",                   label = "Age",              scale = 100),
  list(coef = "RIAGENDRFemale",             label = "Female",           scale = 10),
  list(coef = "as.factor(race4)NH-Black",   label = "NH Black",         scale = 10),
  list(coef = "as.factor(race4)Hispanic",   label = "Hispanic",         scale = 10),
  list(coef = "as.factor(race4)NH-Other",   label = "NH Other",         scale = 10),
  list(coef = "phys.all",                   label = "Physical activity", scale = 100),
  list(coef = "as.factor(smk3)1",           label = "Former smoker",    scale = 10),
  list(coef = "as.factor(smk3)2",           label = "Current smoker",   scale = 10),
  list(coef = "as.factor(hvy.drk)1",        label = "Some-day drinker", scale = 10),
  list(coef = "as.factor(hvy.drk)2",        label = "Heavy drinker",    scale = 100),
  list(coef = "DXDTOFAT",                   label = "Total fat",        scale = 1e6)
)

ci_rows <- list()
for (vm in var_map) {
  rn  <- rownames(table4)
  idx <- which(rn == vm$coef)
  if (length(idx) == 0) next
  for (cm in ci_map) {
    b   <- table4[idx, cm$beta]
    se  <- sqrt(table4[idx, cm$var])
    ci_rows[[length(ci_rows) + 1]] <- data.frame(
      variable = vm$label,
      method   = cm$label,
      estimate = b  * vm$scale,
      lower    = (b - z * se) * vm$scale,
      upper    = (b + z * se) * vm$scale,
      stringsAsFactors = FALSE
    )
  }
}
CI_table_tidy <- do.call(rbind, ci_rows)

# ==============================================================================
# 11. Write results
# ==============================================================================
saveRDS(table4,        "results/table4.rds")
saveRDS(CI_table_tidy, "results/CI_table_tidy.rds")

write.csv(table4, "results/model_results_table4.csv", row.names = TRUE)
write.csv(CI_table_tidy, "results/CI_table_tidy.csv", row.names = FALSE)



cat("\nResults written to results/:\n")
cat("  table4.rds, CI_table_tidy.rds\n")
# cat("  model_results_table4.csv, CI_table_tidy.csv, CI_table_wide.csv\n")
# cat("  Model_Bias_Final.csv (from step 6)\n")
cat("Script 02 complete.\n")
