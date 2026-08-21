# =============================================================================
# 04_rerun_table1_simulation_optional.R  [SLOW ROUTE -- OPTIONAL]
#
# Reruns the Monte Carlo simulation that produced the canonical results for
# Table 1 (tab:bias_vr_full). This is the SLOW ROUTE (NSIMU = 1000 takes
# several hours on a typical workstation).
#
# Rerun results are written to results/rerun_outputs/table1/  and will NOT
# overwrite the canonical results in results/canonical/table1/.
#
# IMPORTANT: These rerun scripts use sanitized integer seeds to make reruns
# deterministic. Because historical canonical results were produced under
# earlier seed/path handling and post-processing, rerun summaries may not
# exactly match the manuscript canonical tables. For exact reproduction of
# manuscript Table 1, use Script 01 (Fast Route) instead.
#
# Run from:  simulation/
#   setwd("path/to/simulation")
#   source("code/04_rerun_table1_simulation_optional.R")
# =============================================================================

rm(list = ls())
library(survey)
library(MASS)
library(writexl)
library(numDeriv)
library(ICC)

# ---------------------------------------------------------------------------
# User-editable parameters
# ---------------------------------------------------------------------------
NSIMU <- 1000    # set to a smaller value (e.g. 10) for a quick test run

# ---------------------------------------------------------------------------
# Seed sanitization
#
# seed.txt stores floating-point values including negative entries.
# Passing such values directly to set.seed() yields non-deterministic
# or environment-dependent behaviour. We sanitize to valid positive integers.
# ---------------------------------------------------------------------------
normalize_seed <- function(x) {
  s <- suppressWarnings(as.integer(trunc(x)))
  if (is.na(s)) stop(paste("Invalid seed value:", x))
  s <- abs(s)
  if (s == 0L) s <- 1L
  s
}

seed_raw <- read.table("code/helper_functions/seed.txt", header = TRUE)
seed_int <- as.data.frame(lapply(seed_raw, function(col) sapply(col, normalize_seed)))

# Save the sanitized seed matrix for traceability
OUT_DIR <- "results/rerun_outputs/table1"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
write.csv(seed_int, file.path(OUT_DIR, "sanitized_seed_used.csv"), row.names = FALSE)
cat("Sanitized seed matrix saved to:", file.path(OUT_DIR, "sanitized_seed_used.csv"), "\n")

# ---------------------------------------------------------------------------
# Paths (relative to simulation/ working directory)
# ---------------------------------------------------------------------------
source("code/helper_functions/subfunctions20250403_V1_0.R")

# ---------------------------------------------------------------------------
# Population generation (same design as manuscript Section 6)
# n2 = 6000 to match srs_pps_6000_0.1 canonical results
# ---------------------------------------------------------------------------
set.seed(94833)

N     <- 1000000
M.psu <- 10000
psu.size <- N / M.psu

n1 <- 2000; f1 <- n1 / N
n2 <- 6000; f2 <- n2 / N   # changed from 2000 to match Table 1 canonical
m.psu <- 60
n.psu <- m.psu * round(n2 / n1, 0)

x1 <- as.numeric(runif(N) < 0.3)

x_mu    <- c(0, 0)
sd_x    <- c(1, 1)
x_rho   <- matrix(c(1.0, 0.1, 0.1, 1.0), 2, 2)
x_sigma <- x_rho * outer(sd_x, sd_x)
x_mtrx  <- as.data.frame(mvrnorm(N, mu = x_mu, Sigma = x_sigma))
x2 <- x_mtrx[, 1]
x3 <- x_mtrx[, 2]

# Surrogate covariates: a[1]=0.67 gives rho~0.80; a[2]=1.55 gives rho~0.50
a  <- c(0.67, 1.55)
e  <- rnorm(N, mean = 0, sd = 0.5)
x3.s <- as.data.frame(sapply(seq_along(a),
                             function(i) x3 + a[i] * log(abs(x3 + e))))
names(x3.s) <- paste0("x3.s", seq_along(a))

betas.y <- c(-3.5, 1, 0.5, 0.5, 0.2)
odds.y  <- exp(cbind(1, x1, x2, x3, x1 * x3) %*% betas.y)
p.y     <- odds.y / (1 + odds.y)
y <- get.y(betas = betas.y,
           design.x = cbind(1, x1, x2, x3, x1 * x3),
           ICCy = 0.1, M.psu = M.psu, psu.size = psu.size)

gamma1   <- c(0, 0.05, 0.1, 0, 0.01)
odds.s1  <- exp(cbind(x1, x2, y, y * x1, y * x2) %*% gamma1)

gamma2   <- c(0.05, 0, 0.1, 0.01, 0)   # gamma2[3]=0.1 => folder suffix "_0.1"
odds.s2  <- exp(cbind(x1, x2, y, y * x1, y * x2) %*% gamma2)

pop <- data.frame(id = 1:N, x1, x2, x3, y, x3.s, odds.s1, odds.s2)

base_psu_size <- floor(N / M.psu)
remaining     <- N %% M.psu
psu_assignments <- rep(1:M.psu, each = base_psu_size)
if (remaining > 0) {
  extra_psus      <- sample(1:M.psu, remaining, replace = FALSE)
  psu_assignments <- c(psu_assignments, extra_psus)
}
pop$psu <- sample(psu_assignments, N, replace = FALSE)

# ---------------------------------------------------------------------------
# Method labels (16 methods)
# ---------------------------------------------------------------------------
method <- c(
  "s1.wt.r",    "s1.wt.f",    "s1.wt.m1",    "s1.wt.m2",
  "s2.wt.r",    "s2.wt.f",    "s2.wt.m1",    "s2.wt.m2",
  "s1.cwt.i.r", "s1.cwt.i.f", "s1.cwt.i.m1", "s1.cwt.i.m2",
  "s1.cwts.r",  "s1.cwts.f",  "s1.cwts.m1",  "s1.cwts.m2"
)

fit.y <- rbind(
  c("y~x1+x2+x3+x1:x3", "y~x1+x2+x3.s+x1:x3.s", "y~x1+x2"),
  c("y~x1+x2+x3",        "y~x1+x2+x3.s",          "y~x1+x2")
)

# True FP regression coefficients (fit on full population)
beta.pop.t <- glm(fit.y[1, 1], family = "binomial", pop)$coef
beta.pop.m <- glm(fit.y[2, 1], family = "binomial", pop)$coef
beta.pop.m[setdiff(names(beta.pop.t), names(beta.pop.m))] <- NA
beta.pop   <- rbind(beta.pop.t, beta.pop.m)

# ---------------------------------------------------------------------------
# Monte Carlo loop
# ---------------------------------------------------------------------------
for (i in 1:nrow(fit.y)) {
  beta.pop.fit <- beta.pop[i, ]
  beta.pop.fit <- beta.pop.fit[!is.na(beta.pop.fit)]
  n_beta.fit   <- length(beta.pop.fit)

  beta.est  <- array(0, c(NSIMU, n_beta.fit, length(method)))
  beta.var1 <- array(0, c(NSIMU, n_beta.fit, length(method)))

  for (simu in 1:NSIMU) {
    samp1 <- samp.slct(seed = seed_int[simu, 1], fnt.pop = pop,
                       n = n1, psu.name = "psu", m.psu = m.psu,
                       dsgn = "srs-pps", size = "odds.s1")
    ds1   <- svydesign(ids = ~psu, strata = NULL, weights = ~wt, data = samp1)

    samp2 <- samp.slct(seed = seed_int[simu, 2], fnt.pop = pop,
                       n = n2, psu.name = "psu", m.psu = n.psu,
                       dsgn = "srs-pps", size = "odds.s2")
    ds2   <- svydesign(ids = ~psu, strata = NULL, weights = ~wt, data = samp2)

    S    <- rbind(samp1, samp2); n.S <- nrow(S)
    S$wt <- c(samp1$wt * length(samp1$id) / n.S,
              samp2$wt * length(samp2$id) / n.S)
    ds.S <- svydesign(ids = ~psu, strata = NULL, weights = ~wt, data = S)

    calib.beta <- function(fm.idx) {
      s1.fit    <- svyglm(fit.y[i, fm.idx], ds1, family = "binomial")
      beta.1    <- s1.fit$coeff
      beta.var.1 <- diag(vcov(s1.fit))

      s2.fit    <- svyglm(fit.y[i, fm.idx], ds2, family = "binomial")
      beta.2    <- s2.fit$coeff
      beta.var.2 <- diag(vcov(s2.fit))

      S.fit   <- svyglm(fit.y[i, fm.idx], ds.S, family = "binomial")
      X.S     <- model.matrix(S.fit)
      ui      <- (samp1$y - S.fit$fitted.values[1:length(samp1$id)]) *
                   X.S[1:length(samp1$id), ]
      samp1$cwt.i <- samp1$wt *
        greg_f(wt0 = samp1$wt, v.mtx0 = ui,
               VS.hat = rep(0, ncol(X.S)))$f
      samp1$cwt.i[samp1$cwt.i < 0] <- 0
      ds1.cwt.i <- svydesign(ids = ~psu, strata = NULL,
                              weights = ~cwt.i, data = samp1)
      beta.i <- svyglm(fit.y[i, 1], ds1.cwt.i, family = "binomial")$coeff

      beta.var.i.fit <- calibrate_beta_variance(
        ds_s = ds1.cwt.i, ds_S = ds.S,
        fit.y[i, fm.idx], fit.y[i, 1], ui, samp1$wt,
        structure = "combined", var_beta_star = NULL,
        data = NULL, beta = NULL, ds2 = ds2,
        n1 = length(samp1$id), n2 = length(samp2$id), ds1 = ds1)
      beta.var.i <- diag(beta.var.i.fit$var_beta_calib)

      x.s1    <- model.matrix(lm(fit.y[i, fm.idx], data = samp1))
      odds.hat <- exp(x.s1 %*% s2.fit$coeff)
      p.hat   <- odds.hat / (1 + odds.hat)
      ui      <- c(samp1$y - p.hat) * x.s1
      samp1$cwt.s <- samp1$wt *
        greg_f(wt0 = samp1$wt, v.mtx0 = ui,
               VS.hat = rep(0, ncol(x.s1)))$f
      samp1$cwt.s[samp1$cwt.s < 0] <- 0
      ds1.cwt.s <- svydesign(ids = ~psu, strata = NULL,
                              weights = ~cwt.s, data = samp1)
      beta.s <- svyglm(fit.y[i, 1], ds1.cwt.s, family = "binomial")$coeff

      beta.var.s1.fit <- calibrate_beta_variance(
        ds_s = ds1.cwt.s, ds_S = ds1,
        fit.y[i, fm.idx], fit.y[i, 1], ui, samp1$wt,
        structure = "model-based", var_beta_star = vcov(s2.fit),
        data = samp1, beta = s2.fit$coefficients,
        ds2 = NULL, n1 = NULL, n2 = NULL, ds1 = NULL)
      beta.var.s1 <- diag(beta.var.s1.fit$var_beta_calib)

      if (length(beta.i) != length(beta.1)) {
        beta.1[setdiff(names(beta.i), names(beta.1))]    <- NA
        beta.2[setdiff(names(beta.i), names(beta.2))]    <- NA
        beta.var.1[setdiff(names(beta.var.i), names(beta.var.1))] <- NA
        beta.var.2[setdiff(names(beta.var.i), names(beta.var.2))] <- NA
      }

      list(beta      = cbind(beta.1, beta.2, beta.i, beta.s),
           beta.var1 = cbind(beta.var.1, beta.var.2, beta.var.i, beta.var.s1))
    }

    calib.r <- calib.beta(fm.idx = 3)
    beta.est [simu, , c(1, 5,  9, 13)] <- calib.r$beta
    beta.var1[simu, , c(1, 5,  9, 13)] <- calib.r$beta.var1

    calib.f <- calib.beta(fm.idx = 1)
    beta.est [simu, , c(2, 6, 10, 14)] <- calib.f$beta
    beta.var1[simu, , c(2, 6, 10, 14)] <- calib.f$beta.var1

    for (fit.k in seq_along(a)) {
      x3.s.lab    <- paste0("x3.s", fit.k)
      samp1$x3.s  <- samp1[, x3.s.lab]
      ds1  <- update(ds1,  x3.s = samp1[, x3.s.lab])
      ds2  <- update(ds2,  x3.s = samp2[, x3.s.lab])
      ds.S <- update(ds.S, x3.s = c(samp1[, x3.s.lab], samp2[, x3.s.lab]))
      calib.m <- calib.beta(fm.idx = 2)
      beta.est [simu, , c(1, 5, 9, 13) + (1 + fit.k)] <- calib.m$beta
      beta.var1[simu, , c(1, 5, 9, 13) + (1 + fit.k)] <- calib.m$beta.var1
    }
    cat("Simulation iteration:", simu, "\n")
  }

  # --- Summarise results ---
  beta.coverage <- matrix(0, nrow = length(method), ncol = n_beta.fit)
  rownames(beta.coverage) <- method
  colnames(beta.coverage) <- paste0("beta.", names(beta.pop.fit))

  for (meth.k in seq_along(method)) {
    for (beta.k in seq_len(n_beta.fit)) {
      ests  <- beta.est[, beta.k, meth.k]
      vars  <- beta.var1[, beta.k, meth.k]
      se    <- sqrt(vars)
      lower <- ests - 1.96 * se
      upper <- ests + 1.96 * se
      beta.coverage[meth.k, beta.k] <-
        mean((beta.pop.fit[beta.k] >= lower) & (beta.pop.fit[beta.k] <= upper))
    }
  }

  beta.means     <- sapply(1:n_beta.fit, function(k) apply(beta.est[, k, ],  2, mean))
  beta.var1.mean <- sapply(1:n_beta.fit, function(k) apply(beta.var1[, k, ], 2, mean))
  beta.evars     <- sapply(1:n_beta.fit, function(k) apply(beta.est[, k, ],  2, var))
  beta.bias      <- sweep(beta.means, 2, beta.pop.fit, FUN = "-") /
                    matrix(rep(betas.y, each = nrow(beta.means)), ncol = length(betas.y))
  beta.bias.true <- sweep(beta.means, 2, beta.pop.fit, FUN = "-")
  ratio1         <- beta.var1.mean / beta.evars

  coef_nm <- paste0("beta.", names(beta.pop.fit))
  colnames(beta.means) <- colnames(beta.evars) <- colnames(beta.var1.mean) <-
    colnames(ratio1) <- coef_nm
  rownames(beta.means) <- rownames(beta.evars) <- rownames(beta.var1.mean) <-
    rownames(ratio1)    <- method

  # --- Save ---
  fnames <- c(
    paste0("beta_means_fit_",                i, "_n2_", n2, ".csv"),
    paste0("beta_empirical_variances_fit_",  i, "_n2_", n2, ".csv"),
    paste0("beta_variance_candidate1_fit_",  i, "_n2_", n2, ".csv"),
    paste0("ratio_candidate1_fit_",          i, "_n2_", n2, ".csv"),
    paste0("beta_bias_fit_",                 i, "_n2_", n2, ".csv"),
    paste0("beta_bias_true_fit_",            i, "_n2_", n2, ".csv"),
    paste0("beta_coverage_fit_",             i, "_n2_", n2, ".csv")
  )
  mats <- list(beta.means, beta.evars, beta.var1.mean, ratio1,
               beta.bias,  beta.bias.true, beta.coverage)
  for (j in seq_along(mats)) {
    write.csv(mats[[j]], file.path(OUT_DIR, fnames[j]), row.names = TRUE)
  }
  cat("Saved rerun results for fit =", i, "\n")
}
