# =============================================================================
# 05_rerun_table2_simulation_optional.R  [SLOW ROUTE -- OPTIONAL]
#
# Reruns the Monte Carlo simulation that produced the canonical results for
# Table 2 (tab:simulation_theta). This is the SLOW ROUTE (NSIMU = 1000 is
# very slow: each iteration fits glm() on N = 1,000,000 observations).
#
# Rerun results are written to results/rerun_outputs/table2/ and will NOT
# overwrite the canonical results in results/canonical/table2/.
#
# IMPORTANT: These rerun scripts use sanitized integer seeds to make reruns
# deterministic. Because historical canonical results were produced under
# earlier seed/path handling and post-processing, rerun summaries may not
# exactly match the manuscript canonical tables. For exact reproduction of
# manuscript Table 2, use Script 02 (Fast Route) instead.
#
# Run from:  simulation/
#   setwd("path/to/simulation")
#   source("code/05_rerun_table2_simulation_optional.R")
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
NSIMU <- 1000    # set to a smaller value (e.g. 5) for a quick test run

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
OUT_DIR <- "results/rerun_outputs/table2"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
write.csv(seed_int, file.path(OUT_DIR, "sanitized_seed_used.csv"), row.names = FALSE)
cat("Sanitized seed matrix saved to:", file.path(OUT_DIR, "sanitized_seed_used.csv"), "\n")

# ---------------------------------------------------------------------------
# Paths (relative to simulation/ working directory)
# ---------------------------------------------------------------------------
source("code/helper_functions/subfunctions20250417_V1_3.R")

# ---------------------------------------------------------------------------
# Population generation (same design as manuscript; gamma2[3] = 0.1)
# ---------------------------------------------------------------------------
set.seed(94833)

N     <- 1000000
M.psu <- 10000
psu.size <- N / M.psu

n1 <- 2000; f1 <- n1 / N
n2 <- 6000; f2 <- n2 / N
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

a  <- c(0.67, 1.55)
e  <- rnorm(N, mean = 0, sd = 0.5)
x3.s <- as.data.frame(sapply(seq_along(a),
                             function(i) x3 + a[i] * log(abs(x3 + e))))
names(x3.s) <- paste0("x3.s", seq_along(a))
cat("Surrogate correlations:\n"); print(cor(x3, x3.s))

betas.y <- c(-3.5, 1, 0.5, 0.5, 0.2)
odds.y  <- exp(cbind(1, x1, x2, x3, x1 * x3) %*% betas.y)
p.y     <- odds.y / (1 + odds.y)
y <- get.y(betas = betas.y,
           design.x = cbind(1, x1, x2, x3, x1 * x3),
           ICCy = 0.1, M.psu = M.psu, psu.size = psu.size)

gamma1  <- c(0, 0.05, 0.1, 0, 0.01)
odds.s1 <- exp(cbind(x1, x2, y, y * x1, y * x2) %*% gamma1)

gamma2  <- c(0.05, 0, 0.1, 0.01, 0)
odds.s2 <- exp(cbind(x1, x2, y, y * x1, y * x2) %*% gamma2)

pop <- data.frame(id = 1:N, x1, x2, x3, y, x3.s, odds.s1, odds.s2)

base_psu_size   <- floor(N / M.psu)
remaining       <- N %% M.psu
psu_assignments <- rep(1:M.psu, each = base_psu_size)
if (remaining > 0) {
  extra_psus      <- sample(1:M.psu, remaining, replace = FALSE)
  psu_assignments <- c(psu_assignments, extra_psus)
}
pop$psu <- sample(psu_assignments, N, replace = FALSE)

# ---------------------------------------------------------------------------
# Method labels (32 methods: 4 scenarios x 8 estimators)
# ---------------------------------------------------------------------------
method_block <- c(
  "beta.1", "beta.2", "beta.i",
  "beta.if.S", "beta.if.S.pop",
  "beta.s",
  "beta.if.s2", "beta.if.s2.pop"
)
scenario <- c("r", "f", "m1", "m2")
method   <- as.vector(sapply(scenario, function(sc) paste0(sc, ".", method_block)))

fit.y <- rbind(
  c("y~x1+x2+x3+x1:x3", "y~x1+x2+x3.s+x1:x3.s", "y~x1+x2"),
  c("y~x1+x2+x3",        "y~x1+x2+x3.s",          "y~x1+x2")
)

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

  beta.est <- array(0, c(NSIMU, n_beta.fit, length(method)))

  idx_block <- function(sc) which(startsWith(method, paste0(sc, ".")))

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
      s1.fit <- svyglm(fit.y[i, fm.idx], ds1, family = "binomial")
      beta.1 <- s1.fit$coeff
      s2.fit <- svyglm(fit.y[i, fm.idx], ds2, family = "binomial")
      beta.2 <- s2.fit$coeff

      S.fit   <- svyglm(fit.y[i, fm.idx], ds.S, family = "binomial")
      X.S     <- model.matrix(S.fit)
      ui      <- (samp1$y - S.fit$fitted.values[1:length(samp1$id)]) *
                   X.S[1:length(samp1$id), ]
      samp1$cwt.i <- samp1$wt *
        greg_f(wt0 = samp1$wt, v.mtx0 = ui, VS.hat = rep(0, ncol(X.S)))$f
      samp1$cwt.i[samp1$cwt.i < 0] <- 0
      ds1.cwt.i <- svydesign(ids = ~psu, strata = NULL,
                              weights = ~cwt.i, data = samp1)
      beta.i <- svyglm(fit.y[i, 1], ds1.cwt.i, family = "binomial")$coeff

      # Population-based IF (Theta using full population glm)
      X.star.s1   <- model.matrix(lm(fit.y[i, fm.idx], data = samp1))
      fit_star_pop <- glm(fit.y[i, fm.idx], pop, family = "binomial")
      X.pop        <- model.matrix(fit_star_pop)
      p.hat.pop    <- plogis(X.pop %*% fit_star_pop$coeff)
      eta.u        <- drop(X.star.s1 %*% fit_star_pop$coeff)
      p.hat.u      <- plogis(eta.u)
      ui.pop       <- (samp1$y - p.hat.u) * X.star.s1

      odds.hat.S   <- exp(X.S %*% S.fit$coeff)
      p.hat.S      <- odds.hat.S / (1 + odds.hat.S)

      if.out.S <- calib_if_logit(wt0 = S$wt, X = X.S, y = S$y,
                                  p_hat = p.hat.S, u = ui,
                                  VS.hat = rep(0, ncol(X.star.s1)))
      Delta.star.S <- if.out.S$v.mtx0

      s1.fit1   <- svyglm(fit.y[i, 1], ds1, family = "binomial")
      x.s1.out  <- model.matrix(lm(fit.y[i, 1], data = samp1))
      odds.s1.f <- exp(x.s1.out %*% s1.fit1$coeff)
      p.s1.f    <- odds.s1.f / (1 + odds.s1.f)
      if.nostar <- calib_if_logit(wt0 = samp1$wt, X = x.s1.out, y = samp1$y,
                                   p_hat = p.s1.f, VS.hat = rep(0, ncol(x.s1.out)))
      Delta <- if.nostar$v.mtx0

      X.nonstar.s1 <- model.matrix(lm(fit.y[i, 1], data = samp1))
      fit_pop      <- glm(fit.y[i, 1], pop, family = "binomial")
      X.pop.ns     <- model.matrix(fit_pop)
      p.ns.pop     <- plogis(X.pop.ns %*% fit_pop$coeff)
      p.ns.u       <- plogis(drop(X.nonstar.s1 %*% fit_pop$coeff))
      u.pop        <- (samp1$y - p.ns.u) * X.nonstar.s1
      if.nostar.pop <- calib_if_logit(wt0 = 1, X = X.pop.ns, y = pop$y,
                                       p_hat = p.ns.pop, u = u.pop,
                                       VS.hat = rep(0, ncol(X.nonstar.s1)))
      Delta.pop <- if.nostar.pop$v.mtx

      if.out.pop <- calib_if_logit(wt0 = 1, X = X.pop, y = pop$y,
                                    p_hat = p.hat.pop, u = ui.pop,
                                    VS.hat = rep(0, ncol(X.star.s1)))
      Delta.star.pop <- if.out.pop$v.mtx0

      w <- samp1$wt
      A <- crossprod(Delta.star.S   * w, Delta.star.S);   B <- crossprod(Delta.star.S   * w, Delta)
      Theta.hat.S   <- solve(A, B)
      A <- crossprod(Delta.star.pop * w, Delta.star.pop); B <- crossprod(Delta.star.pop * w, Delta.pop)
      Theta.hat.pop <- solve(A, B)

      beta.if.S     <- s1.fit1$coeff + t(Theta.hat.S)   %*% (S.fit$coeff - beta.1)
      beta.if.S.pop <- s1.fit1$coeff + t(Theta.hat.pop) %*% (S.fit$coeff - beta.1)

      x.s1.star <- model.matrix(lm(fit.y[i, fm.idx], data = samp1))
      odds.hat  <- exp(x.s1.star %*% s2.fit$coeff)
      p.hat     <- odds.hat / (1 + odds.hat)
      ui        <- c(samp1$y - p.hat) * x.s1.star
      samp1$cwt.s <- samp1$wt *
        greg_f(wt0 = samp1$wt, v.mtx0 = ui, VS.hat = rep(0, ncol(x.s1.star)))$f
      samp1$cwt.s[samp1$cwt.s < 0] <- 0
      ds1.cwt.s <- svydesign(ids = ~psu, strata = NULL,
                              weights = ~cwt.s, data = samp1)
      beta.s <- svyglm(fit.y[i, 1], ds1.cwt.s, family = "binomial")$coeff

      if.out    <- calib_if_logit(wt0 = samp1$wt, X = x.s1.star, y = samp1$y,
                                   p_hat = p.hat, VS.hat = rep(0, ncol(x.s1.star)))
      Delta.star <- if.out$v.mtx0
      A <- crossprod(Delta.star * w, Delta.star); B <- crossprod(Delta.star * w, Delta)
      Theta.hat.s2 <- solve(A, B)

      beta.if.s2     <- s1.fit1$coeff + t(Theta.hat.s2)  %*% (beta.2 - beta.1)
      beta.if.s2.pop <- s1.fit1$coeff + t(Theta.hat.pop) %*% (beta.2 - beta.1)

      if (length(beta.i) != length(beta.1)) {
        beta.1[setdiff(names(beta.i), names(beta.1))] <- NA
        beta.2[setdiff(names(beta.i), names(beta.2))] <- NA
      }

      list(beta = cbind(beta.1, beta.2, beta.i, beta.if.S, beta.if.S.pop,
                        beta.s, beta.if.s2, beta.if.s2.pop))
    }

    calib.r <- calib.beta(fm.idx = 3)
    beta.est[simu, , idx_block("r")] <- calib.r$beta

    calib.f <- calib.beta(fm.idx = 1)
    beta.est[simu, , idx_block("f")] <- calib.f$beta

    for (fit.k in seq_along(a)) {
      x3.s.lab   <- paste0("x3.s", fit.k)
      samp1$x3.s <- samp1[, x3.s.lab]
      ds1  <- update(ds1,  x3.s = samp1[, x3.s.lab])
      ds2  <- update(ds2,  x3.s = samp2[, x3.s.lab])
      ds.S <- update(ds.S, x3.s = c(samp1[, x3.s.lab], samp2[, x3.s.lab]))
      pop$x3.s <- pop[, x3.s.lab]
      calib.m <- calib.beta(fm.idx = 2)
      sc <- if (fit.k == 1) "m1" else "m2"
      beta.est[simu, , idx_block(sc)] <- calib.m$beta
    }
    cat("Simulation iteration:", simu, "\n")
  }

  # --- Summarise ---
  beta.means     <- sapply(1:n_beta.fit, function(k) apply(beta.est[, k, ], 2, mean))
  beta.evars     <- sapply(1:n_beta.fit, function(k) apply(beta.est[, k, ], 2, var))
  beta.bias.true <- sweep(beta.means, 2, beta.pop.fit, FUN = "-")

  coef_nm <- paste0("beta.", names(beta.pop.fit))
  colnames(beta.means) <- colnames(beta.evars) <- colnames(beta.bias.true) <- coef_nm
  rownames(beta.means) <- rownames(beta.evars) <- method

  fnames <- c(
    paste0("beta_means_fit_",               i, "_n2_", n2, ".csv"),
    paste0("beta_empirical_variances_fit_",  i, "_n2_", n2, ".csv"),
    paste0("beta_bias_true_fit_",            i, "_n2_", n2, ".csv")
  )
  mats <- list(beta.means, beta.evars, beta.bias.true)
  for (j in seq_along(mats)) {
    write.csv(mats[[j]], file.path(OUT_DIR, fnames[j]), row.names = TRUE)
  }
  cat("Saved rerun results for fit =", i, "\n")
}
