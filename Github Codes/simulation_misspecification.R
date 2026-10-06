############################################################
## simulation_misspecification.R
##
## MR-CCC SIMULATION UNDER MODEL MISSPECIFICATION  (scenarios S4-S7)
##
## The simulations in simulation_mrccc.R (S1-S3) generate data from exactly
## the model that MR-CCC fits, with strong instruments. They therefore test
## estimation, not robustness. This script adds the four departures from the
## assumed model that matter most for the real-data application:
##
##   S4  WEAK INSTRUMENTS
##       Instrument strength is calibrated to a TARGET first-stage F rather
##       than to a fixed effect size, so results map directly onto the F
##       values observed in the OneK1K analysis.
##
##   S5  CORRELATED INSTRUMENTS
##       Instruments are given an AR(1) correlation structure mimicking
##       linkage disequilibrium. This probes the correlation that remains
##       among cis-SNPs after LD clumping (up to r^2 = 0.8, i.e. |r| = 0.9).
##
##   S6  INVALID / PLEIOTROPIC INSTRUMENTS
##       A fraction of instruments is given a DIRECT effect on the outcome,
##       violating the exclusion restriction. Both the fraction and the
##       magnitude are varied to locate the breakdown point.
##
##   S7  NONLINEAR LIGAND EFFECT
##       The outcome is generated with a quadratic ligand effect, testing
##       whether the linear-plus-interaction working model spuriously
##       activates the inclusion indicator. A threshold form is also
##       available in generate_data_mis() but is not used in the paper.
##
## Every scenario is run under THREE arms, mirroring scenarios S1-S3 of the
## correctly specified study in the same order and with the same effect
## sizes:
##
##   null    beta_X = 0,   beta_XZ = 0    false positive rate
##   signal  beta_X = 0.3, beta_XZ = 0.3  power
##   main    beta_X = 0.3, beta_XZ = 0    power, plus the question the
##           interaction model must answer for itself: with a real main
##           effect but NO true interaction, does the misspecification
##           manufacture a spurious one? That appears as a systematic
##           beta_XZ bias against a truth of zero.
##
## The four methods compared, the scoring conventions, the rejection rules,
## the number of replicates and the prior settings are identical to
## simulation_mrccc.R, so S4-S7 are directly comparable with S1-S3. The
## MR-CCC chain is longer here (one chain of 400,000 sweeps against 20,000
## in S1-S3); see the note on chain length under USER SETTINGS.
##
## PREREQUISITES (in this order):
##   Rcpp::sourceCpp("Github Codes/mr_ccc_gibbs.cpp")
##   source("Github Codes/simulation_mrccc.R")
##       # provides run_methods(), mr_bma_simple(), generate_data(), etc.
##   source("Github Codes/simulation_misspecification.R")
##
## OUTPUT (written to Results/):
##   simulation_S4_weak.csv          per-replicate results, weak instruments
##   simulation_S5_correlated.csv    per-replicate results, correlated IVs
##   simulation_S6_pleiotropy.csv    per-replicate results, invalid IVs
##   simulation_S7_nonlinear.csv     per-replicate results, nonlinear effect
##   simulation_S4_S7_summary.csv    aggregated summary across all scenarios
##
## FIGURES (written to Plots/ via save_fig() from figure_style.R):
##   sim_S4_weak.pdf, sim_S5_correlated.pdf, sim_S6_pleiotropy.pdf,
##   sim_S7_nonlinear.pdf          rejection rate against the swept setting,
##                                 faceted by arm and sample size
##   sim_S4_S7_score_box.pdf,      per-replicate estimates at representative
##   sim_S4_S7_betaX_box.pdf,      settings of each scenario, with the truth
##   sim_S4_S7_betaXZ_box.pdf      marked per arm
##   All are rebuildable from the CSVs above with regenerate_figures.R.
############################################################
library(parallel)
library(dplyr)
library(tidyr)

## The helper functions defined in simulation_mrccc.R are reused rather than
## duplicated, so that the two scripts cannot drift apart.
if (!exists("run_methods") || !exists("generate_data")) {
  stop("run_methods()/generate_data() not found.\n",
       "  Source the main simulation first, from the project root:\n",
       "    Rcpp::sourceCpp(\"Github Codes/mr_ccc_gibbs.cpp\")\n",
       "    source(\"Github Codes/simulation_mrccc.R\")")
}

############################################################
## USER SETTINGS
############################################################

SEED_MIS <- 20260921          # reproducibility

## Sample sizes. The default is deliberately smaller than the S1-S3 grid:
## misspecification behaviour is of interest at sample sizes comparable to
## the real analysis (n = 651 donors), and the full grid multiplied by the
## number of misspecification settings is expensive. Extend if required.
mis_sample_sizes <- c(500, 1000)

## Replicates per cell; matches R_reps in simulation_mrccc.R, where the
## Monte Carlo standard error of a rejection rate is worked through. It
## matters more here than there: S4-S7 are read as calibration statements
## ("is the false-positive rate still 0.05?"), and at 20 replicates the
## error on such a rate is +/-0.049, which is the nominal level itself.
R_reps_mis  <- 100            # replicates per cell (matches S1-S3)
pip_thr_mis <- 0.5            # discovery threshold (matches S1-S3)
alpha_sig_mis <- 0.05

## ---- Chain length ----------------------------------------------------------
## Scenarios S1-S3 use 20,000 iterations, which is sufficient there because
## the simulated posterior inclusion probabilities mostly sit far from the
## 0.5 threshold (well below it under the null, near 1 under signal), so
## chain noise rarely changes a rejection decision.
##
## That argument does not carry over to misspecification. Departures from
## the assumed model make the posterior less decisive and pull the inclusion
## probabilities toward the middle of the unit interval, which is exactly
## where Monte Carlo error decides rejections.
##
## All four scenarios are therefore run on the same NUMBER OF DRAWS as the
## real data, approximately 400,000. Applying one setting across the study is
## simpler to state than a scenario-by-scenario rule, and it removes any
## question of whether a difference between scenarios reflects the
## misspecification or the sampler. At this length the Monte Carlo standard
## error of the PIP is small relative to its distance from 0.5: in the
## saved results it is at most about 0.03 in every replicate where it is
## defined (it is undefined when the chain never changes state).
## MR_pip_mcse is written for every replicate so this can be checked from
## the output.
##
## The draws are taken as ONE chain of 400,000 here, whereas the real data
## pools FOUR chains of 100,000. The distinction is deliberate and is a
## matter of what each analysis reports, not of precision: combining chains in
## quadrature reproduces the single-chain standard error exactly once the
## total is matched, so the two configurations are equally precise. Several
## chains are used for the real data because the Gelman-Rubin statistic is
## reported there for every triplet; a single chain is used here because
## between-chain diagnostics on 10,800 independent fits would be reported
## nowhere, and quadrupling the number of fits to compute a statistic that
## never appears would buy nothing.
n_iter_mis  <- 400000
burn_in_mis <- 2000

pG_mis <- 5; pH_mis <- 5; pV_mis <- 3
conf_mis <- 0.7               # confounder loading (matches S1-S3)

## Signal used for the "signal" arms; matches scenario S2 of the main study
## so that the weak-instrument results can be read directly against it.
beta_X_sig  <- 0.3
beta_XZ_sig <- 0.3
beta_Z_all  <- 0.5

## Target first-stage F values for S4. The real-data analysis has median
## F ~ 1.5 (ligand) and ~1.8 (receptor), with first-stage F below 10
## throughout, while the S1-S3 design (pi = 0.5) corresponds to a
## calibration F of about 70 at n = 500, rising in proportion to n at the
## larger S1-S3 sample sizes. The grid therefore spans the observed range and
## the regime of scenarios S1-S3 at n = 500.
target_F_grid <- c(1, 2, 5, 10, 25, 70)

## S5: AR(1) correlation among instrument columns.
iv_corr_grid <- c(0, 0.3, 0.6, 0.9)

## S6: fraction of instruments with a direct effect on Y, and its magnitude.
pleio_frac_grid  <- c(0.2, 0.4)
pleio_delta_grid <- c(0.1, 0.3)

## S7: nonlinear forms for the ligand effect.
## Quadratic only. A threshold form is implemented in generate_data_mis()
## and can be added here, but a single clear departure from linearity
## isolates the effect of misspecification; piling on harder forms would
## obscure rather than sharpen the point.
nl_form_grid <- c("quadratic")

## STRENGTH of the nonlinearity, swept rather than fixed. The quadratic term
## enters as  nl_coef * (X^2 - E[X^2]), so nl_coef is on the same scale as
## beta_X, which is 0.3 throughout. The grid therefore runs from a curvature
## one sixth the size of the linear effect -- a mild departure from linearity
## of the kind a real ligand-response curve might show, and the setting the
## headline statement is made at -- up to a curvature as large as the linear
## effect itself, which is severe.
##
## Sweeping the strength rather than fixing one value matters for how the
## result can be read. A single severe setting describes only an extreme
## departure, and a single mild setting says nothing about stronger ones. A
## dose-response covers both: it shows whether the working model is affected
## while the departure is mild, identifies where it begins to manufacture an
## interaction, and reports how large that spurious interaction is at each
## level, so a reader can judge which regime their own application
## resembles. The model itself is not changed in response to this scenario;
## the scenario characterises the range over which the linear working model
## remains adequate.
nl_coef_grid <- c(0.05, 0.1, 0.2, 0.3)

############################################################
## INSTRUMENT STRENGTH CALIBRATION
##
## With G ~ N(0, I) and pi_X = pi * 1, the sender equation is
##   X = G pi_X + V alpha_X + conf * U + e,   e ~ N(0,1)
## so that, writing rest = Var(V alpha_X) + conf^2 + Var(e),
##   R^2 / (1 - R^2) = pG * pi^2 / rest
## and the calibration uses
##   F = (R^2 / (1 - R^2)) * (n - pG - pV - 1) / pG
##     = pi^2 * (n - pG - pV - 1) / rest.
##
## This is not exactly the population partial F of G given V: its
## denominator `rest` includes the covariate variance Var(V alpha_X), which
## the partial F (V already in the model) excludes. The realised partial F
## recorded in F_observed is therefore larger than the target -- a median of
## about 85 at a target of 70 -- and the realised values are reported
## alongside the targets.
##
## Inverting gives the pi that attains a target F at a given n. Note that F
## is held FIXED as n grows, which is the weak-instrument asymptotic regime;
## holding pi fixed instead would make every design strong as n increases and
## would miss the behaviour of interest.
############################################################
pi_for_target_F <- function(target_F, n,
                            pG = pG_mis, pV = pV_mis,
                            alpha_val = 0.3, conf = conf_mis, resid_var = 1) {
  rest <- pV * alpha_val^2 + conf^2 + resid_var
  df2  <- n - pG - pV - 1
  if (df2 <= 0) return(NA_real_)
  sqrt(target_F * rest / df2)
}

# Calibration F (same formula as above) for a given pi. The realised
# partial F of each replicate is computed separately, as F_observed.
F_for_pi <- function(pi_val, n,
                     pG = pG_mis, pV = pV_mis,
                     alpha_val = 0.3, conf = conf_mis, resid_var = 1) {
  rest <- pV * alpha_val^2 + conf^2 + resid_var
  pi_val^2 * (n - pG - pV - 1) / rest
}

############################################################
## DATA GENERATION UNDER MISSPECIFICATION
##
## Extends generate_data() from simulation_mrccc.R with four switches. With
## all switches at their defaults the generator reduces EXACTLY to the S1-S3
## data-generating model, which is checked by the self-test that runs before
## the study.
##
##   pi_val      first-stage effect size (set via pi_for_target_F for S4)
##   iv_corr     AR(1) correlation between adjacent instrument columns (S5)
##   pleio_frac  fraction of sender instruments with a direct Y effect (S6)
##   pleio_delta magnitude of that direct effect (S6)
##   nl_form     "none", "quadratic" or "threshold" ligand effect (S7)
############################################################
generate_data_mis <- function(n,
                              pG = pG_mis, pH = pH_mis, pV = pV_mis,
                              beta_X = 0, beta_XZ = 0, beta_Z = beta_Z_all,
                              conf_strength = conf_mis,
                              pi_val = 0.5,
                              iv_corr = 0,
                              pleio_frac = 0, pleio_delta = 0,
                              nl_form = "none", nl_coef = 0.3) {

  # ---- Instruments, optionally with AR(1) correlation (S5) ----
  ar1_chol <- function(p, rho) {
    if (rho == 0) return(diag(p))
    S <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
    chol(S)
  }
  G <- matrix(rnorm(n * pG), n, pG) %*% ar1_chol(pG, iv_corr)
  H <- matrix(rnorm(n * pH), n, pH) %*% ar1_chol(pH, iv_corr)

  V <- matrix(rnorm(n * pV), n, pV)
  U <- rnorm(n)

  pi_X <- rep(pi_val, pG)
  pi_Z <- rep(pi_val, pH)
  alpha_X <- rep(0.3, pV); alpha_Z <- rep(0.3, pV); alpha_Y <- rep(0.3, pV)

  X <- as.numeric(G %*% pi_X + V %*% alpha_X + conf_strength * U + rnorm(n))
  Z <- as.numeric(H %*% pi_Z + V %*% alpha_Z + conf_strength * U + rnorm(n))

  # ---- Ligand contribution to Y: linear, quadratic or threshold (S7) ----
  # The nonlinear forms are centred so that they do not also shift the mean,
  # which would otherwise confound a nonlinearity effect with an intercept
  # effect.
  lig_term <- switch(
    nl_form,
    none      = beta_X * X,
    quadratic = beta_X * X + nl_coef * (X^2 - mean(X^2)),
    threshold = beta_X * X + nl_coef * (pmax(X, 0) - mean(pmax(X, 0))),
    stop("unknown nl_form: ", nl_form)
  )

  Y <- as.numeric(lig_term +
                    beta_Z  * Z +
                    beta_XZ * (X * Z) +
                    V %*% alpha_Y +
                    conf_strength * U +
                    rnorm(n))

  # ---- Invalid instruments: direct G -> Y path (S6) ----
  # Violates the exclusion restriction for the affected instruments. The
  # first m sender instruments are used, so the violated set is well defined.
  n_bad <- floor(pleio_frac * pG)
  if (n_bad > 0 && pleio_delta != 0) {
    Y <- Y + as.numeric(G[, seq_len(n_bad), drop = FALSE] %*%
                          rep(pleio_delta, n_bad))
  }

  list(X = matrix(X, ncol = 1),
       Z = matrix(Z, ncol = 1),
       Y = matrix(Y, ncol = 1),
       G = G, H = H, V = V,
       n_invalid_iv = n_bad)
}

############################################################
## SCENARIO GRID
##
## Each row is one simulation cell. `arm` distinguishes the null
## (beta_X = beta_XZ = 0) from the signal arm, because the two answer
## different questions: calibration versus power.
############################################################
build_grid <- function() {

  # ---- S4: weak instruments, null and signal arms ----
  s4 <- expand.grid(target_F = target_F_grid,
                    arm      = c("null", "signal", "main"),
                    n        = mis_sample_sizes,
                    stringsAsFactors = FALSE) %>%
    mutate(scenario = "S4", iv_corr = 0,
           pleio_frac = 0, pleio_delta = 0, nl_form = "none")

  # ---- S5: correlated instruments (strong IVs, to isolate the effect) ----
  s5 <- expand.grid(iv_corr = iv_corr_grid,
                    arm     = c("null", "signal", "main"),
                    n       = mis_sample_sizes,
                    stringsAsFactors = FALSE) %>%
    mutate(scenario = "S5", target_F = 70,
           pleio_frac = 0, pleio_delta = 0, nl_form = "none")

  # ---- S6: invalid / pleiotropic instruments ----
  s6 <- expand.grid(pleio_frac  = pleio_frac_grid,
                    pleio_delta = pleio_delta_grid,
                    arm         = c("null", "signal", "main"),
                    n           = mis_sample_sizes,
                    stringsAsFactors = FALSE) %>%
    mutate(scenario = "S6", target_F = 70, iv_corr = 0, nl_form = "none")

  # ---- S7: nonlinear ligand effect ----
  s7 <- expand.grid(nl_form = nl_form_grid,
                    nl_coef = nl_coef_grid,
                    arm     = c("null", "signal", "main"),
                    n       = mis_sample_sizes,
                    stringsAsFactors = FALSE) %>%
    mutate(scenario = "S7", target_F = 70, iv_corr = 0,
           pleio_frac = 0, pleio_delta = 0)

  # S4-S6 are linear, so their nonlinearity strength is zero; the column is
  # added explicitly so that bind_rows() aligns rather than filling with NA.
  s4$nl_coef <- 0
  s5$nl_coef <- 0
  s6$nl_coef <- 0

  bind_rows(s4, s5, s6, s7) %>%
    mutate(cell_id = row_number())
}

############################################################
## ONE REPLICATE
############################################################
run_one_mis <- function(cell, rep_id) {

  n       <- cell$n
  pi_val  <- pi_for_target_F(cell$target_F, n)

  # Three arms, mirroring scenarios S1-S3 of the correctly specified study
  # in the same order and with the same effect sizes:
  #   null    beta_X = 0,   beta_XZ = 0    (S1 analogue: false positive rate)
  #   signal  beta_X = 0.3, beta_XZ = 0.3  (S2 analogue: power)
  #   main    beta_X = 0.3, beta_XZ = 0    (S3 analogue: power, and whether
  #                                         misspecification manufactures a
  #                                         SPURIOUS interaction, visible in
  #                                         the beta_XZ bias against truth 0)
  b_X  <- if (cell$arm %in% c("signal", "main")) beta_X_sig  else 0
  b_XZ <- if (identical(cell$arm, "signal"))     beta_XZ_sig else 0

  dat <- generate_data_mis(
    n           = n,
    beta_X      = b_X,
    beta_XZ     = b_XZ,
    beta_Z      = beta_Z_all,
    pi_val      = pi_val,
    iv_corr     = cell$iv_corr,
    pleio_frac  = cell$pleio_frac,
    pleio_delta = cell$pleio_delta,
    nl_form     = cell$nl_form,
    nl_coef     = cell$nl_coef
  )

  # Identical method suite, scoring and rejection rules as S1-S3.
  fit <- run_methods(dat, pip_thr = pip_thr_mis, alpha_sig = alpha_sig_mis,
                     n_iter = n_iter_mis, burn_in = burn_in_mis)

  # Realised first-stage partial F in this replicate, for reporting
  # alongside the target (they differ by sampling variation and by the
  # covariate term noted under INSTRUMENT STRENGTH CALIBRATION).
  f_obs <- tryCatch({
    d  <- data.frame(x = as.numeric(dat$X), dat$V, dat$G)
    m0 <- lm(x ~ ., data = data.frame(x = d$x, dat$V))
    m1 <- lm(x ~ ., data = d)
    a  <- anova(m0, m1); a[2, "F"]
  }, error = function(e) NA_real_)

  tibble::as_tibble(fit) %>%
    mutate(scenario   = cell$scenario,
           arm        = cell$arm,
           n          = n,
           target_F   = cell$target_F,
           pi_val     = pi_val,
           F_observed = f_obs,
           iv_corr    = cell$iv_corr,
           pleio_frac = cell$pleio_frac,
           pleio_delta= cell$pleio_delta,
           nl_form    = cell$nl_form,
           nl_coef    = cell$nl_coef,
           rep_id     = rep_id,
           # gamma indexes the causal block (beta_X, beta_XZ) jointly, so it
           # is truly 1 whenever either coefficient is nonzero -- including
           # the main-only arm. Matches the rule used for S1-S3.
           gamma_true = as.integer(b_X != 0 || b_XZ != 0))
}

############################################################
## DRIVER
##
## mclapply uses forking and is not available on Windows; set n_cores = 1
## there, or substitute parallel::parLapply with a PSOCK cluster.
##
## The core count comes from .default_sim_cores(), defined in
## simulation_mrccc.R and shared by both studies so that parallelisation is
## configured in one place. It claims every core but one by default, which
## is right when a simulation is the only thing running; assigning
## .mrccc_sim_cores before sourcing caps it when the machine is also busy
## with the twenty-pair real-data sweep or the sensitivity analyses. The
## value actually in use is printed either way.
############################################################
run_misspecification <- function(n_cores = .default_sim_cores()) {

  # Parallel-safe RNG: mclapply then gives each worker its own stream, so
  # the parallel runs are reproducible from the seed.
  RNGkind("L'Ecuyer-CMRG")
  set.seed(SEED_MIS)
  grid <- build_grid()

  jobs <- expand.grid(cell_id = grid$cell_id,
                      rep_id  = seq_len(R_reps_mis)) %>%
    arrange(cell_id, rep_id)

  cat("Misspecification study:", nrow(grid), "cells x", R_reps_mis,
      "replicates =", nrow(jobs), "fits\n")
  cat("  scenarios:", paste(unique(grid$scenario), collapse = ", "), "\n")
  cat("  sample sizes:", paste(mis_sample_sizes, collapse = ", "), "\n")
  cat("  cores:", n_cores, "\n\n")

  t0 <- Sys.time()
  res <- mclapply(seq_len(nrow(jobs)), function(k) {
    cell <- grid[grid$cell_id == jobs$cell_id[k], , drop = FALSE]
    tryCatch(run_one_mis(cell, jobs$rep_id[k]),
             error = function(e) NULL)
  }, mc.cores = n_cores)

  out <- bind_rows(res[!vapply(res, is.null, logical(1))])
  cat("Elapsed:",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
      "minutes;", nrow(out), "of", nrow(jobs), "fits returned\n")
  out
}

############################################################
## SUMMARY
##
## Reports, per cell and method: mean communication score, rejection rate,
## and bias / mean absolute deviation for beta_X and beta_XZ.
##
## Under the NULL arm the rejection rate is the false positive rate, which
## is the quantity that decides whether weak instruments (S4), correlated
## instruments (S5) or pleiotropy (S6) compromise calibration. Under the
## SIGNAL arm it is power.
##
## beta_XZ is reported separately and deserves particular attention: its
## design column is a product of two first-stage projections, so it degrades
## before beta_X as instruments weaken.
############################################################
summarise_mis <- function(res) {
  res %>%
    pivot_longer(
      cols = matches("^(OLS|MVMR|MRBMA|MR)_(score|reject|beta_X|beta_XZ)$"),
      names_to  = c("method", ".value"),
      names_pattern = "^(OLS|MVMR|MRBMA|MR)_(.*)$"
    ) %>%
    mutate(method = recode(method,
                           OLS = "OLS", MVMR = "MVMR",
                           MRBMA = "MR-BMA", MR = "MR-CCC")) %>%
    group_by(scenario, arm, n, target_F, iv_corr, pleio_frac, pleio_delta,
             nl_form, nl_coef, method) %>%
    summarise(
      F_obs_median = median(F_observed, na.rm = TRUE),
      score_mean   = mean(score,  na.rm = TRUE),
      reject_rate  = mean(reject, na.rm = TRUE),
      # Monte Carlo standard error of the rejection rate, sqrt(p(1-p)/R).
      # Reported per cell so that a rate can be read against the nominal
      # 0.05 on its own scale: a departure smaller than this is noise.
      reject_mcse  = sd(reject, na.rm = TRUE) /
                       sqrt(sum(is.finite(reject))),
      # True values per arm: beta_X is nonzero in the signal AND main arms;
      # beta_XZ only in the signal arm. In the main arm the beta_XZ bias is
      # therefore measured against zero, which is exactly the spurious-
      # interaction question that arm exists to answer.
      bias_beta_X  = mean(beta_X  - if (first(arm) %in% c("signal", "main"))
                                      beta_X_sig  else 0, na.rm = TRUE),
      mad_beta_X   = mean(abs(beta_X - if (first(arm) %in% c("signal", "main"))
                                         beta_X_sig else 0), na.rm = TRUE),
      bias_beta_XZ = mean(beta_XZ - if (first(arm) == "signal")
                                      beta_XZ_sig else 0, na.rm = TRUE),
      mad_beta_XZ  = mean(abs(beta_XZ - if (first(arm) == "signal")
                                          beta_XZ_sig else 0), na.rm = TRUE),
      .groups = "drop"
    )
}

############################################################
## SELF-TEST
##
## With every misspecification switch off, generate_data_mis() must reduce to
## the S1-S3 generator. This guards against the extended generator silently
## changing the baseline against which S4-S7 are compared. It runs before the
## study, so that a failure is visible before the long run starts. It is
## independent of the study's results: run_misspecification() selects its
## generator and sets its seed itself, so the draws made here do not reach it.
############################################################
local({
  set.seed(99)
  a <- generate_data_mis(n = 400, beta_X = 0.3, beta_XZ = 0.3, pi_val = 0.5)
  set.seed(99)
  b <- generate_data(n = 400, pG = pG_mis, pH = pH_mis, pV = pV_mis,
                     beta_X = 0.3, beta_XZ = 0.3, beta_Z = beta_Z_all,
                     conf_strength = conf_mis, pi_val = 0.5)
  ok <- isTRUE(all.equal(as.numeric(a$Y), as.numeric(b$Y), tolerance = 1e-10))
  cat("\nSelf-test (defaults reduce to the S1-S3 generator):",
      if (ok) "PASS\n" else
        "FAIL -- generate_data_mis() at default settings does not reduce to generate_data(); investigate before relying on S4-S7\n")
})

############################################################
## RUN AND SAVE
############################################################
res_mis <- run_misspecification()

dir.create("Results", showWarnings = FALSE)
write.csv(res_mis %>% filter(scenario == "S4"),
          "Results/simulation_S4_weak.csv", row.names = FALSE)
write.csv(res_mis %>% filter(scenario == "S5"),
          "Results/simulation_S5_correlated.csv", row.names = FALSE)
write.csv(res_mis %>% filter(scenario == "S6"),
          "Results/simulation_S6_pleiotropy.csv", row.names = FALSE)
write.csv(res_mis %>% filter(scenario == "S7"),
          "Results/simulation_S7_nonlinear.csv", row.names = FALSE)

sum_mis <- summarise_mis(res_mis)
write.csv(sum_mis, "Results/simulation_S4_S7_summary.csv", row.names = FALSE)

cat("\n================ S4: WEAK INSTRUMENTS ================\n")
cat("Null arm -- rejection rate is the FALSE POSITIVE rate.\n")
cat("This is the key table: it shows whether MR-CCC stays calibrated at the\n")
cat("instrument strengths actually observed in the OneK1K analysis.\n\n")
print(as.data.frame(
  sum_mis %>%
    filter(scenario == "S4", arm == "null") %>%
    select(n, target_F, F_obs_median, method, score_mean, reject_rate) %>%
    arrange(n, target_F, method)
))

cat("\n---------------- S4: power (signal arm) ----------------\n")
print(as.data.frame(
  sum_mis %>%
    filter(scenario == "S4", arm == "signal") %>%
    select(n, target_F, method, reject_rate, bias_beta_X, bias_beta_XZ) %>%
    arrange(n, target_F, method)
))

cat("\n------------- S4: main-effect arm (no true interaction) -------------\n")
cat("beta_X = 0.3, beta_XZ = 0. Rejection rate is power for the causal\n")
cat("block; bias_beta_XZ is measured against 0, so a systematic departure\n")
cat("from 0 is a SPURIOUS interaction manufactured by the misspecification.\n\n")
print(as.data.frame(
  sum_mis %>%
    filter(scenario == "S4", arm == "main") %>%
    select(n, target_F, method, reject_rate, bias_beta_X, bias_beta_XZ) %>%
    arrange(n, target_F, method)
))

cat("\n-------- Spurious interaction across ALL scenarios (main arm) --------\n")
cat("OLS and MR-CCC only; MVMR and MR-BMA do not model the interaction.\n\n")
print(as.data.frame(
  sum_mis %>%
    filter(arm == "main", method %in% c("OLS", "MR-CCC")) %>%
    select(scenario, n, target_F, iv_corr, pleio_frac, pleio_delta, nl_form,
           nl_coef, method, bias_beta_XZ, mad_beta_XZ) %>%
    arrange(scenario, n, method)
))

cat("\nWritten to Results/simulation_S4_S7_summary.csv\n")

############################################################
## FIGURES
##
## One figure per scenario, each plotting the quantity that scenario is
## designed to interrogate against the departure being varied. The NULL arm
## and the SIGNAL arm are shown together because they answer different
## questions from the same run: under the null the rejection rate is the
## false positive rate, so calibration is the issue; under signal it is
## power.
##
## Each figure is built by a function in figure_functions.R (sourced
## through figure_style.R), which is the single definition of that figure;
## regenerate_figures.R calls the same functions on the saved summary CSV.
## Theme, the Okabe-Ito method palette, point shapes and export settings
## come from figure_style.R, so these figures match S1-S3 and the real-data
## figures exactly.
############################################################
## figure_style.R is sourced by simulation_mrccc.R, which is a prerequisite
## of this script; sourcing it again here is harmless and makes this
## section runnable on its own once sum_mis exists.
source("Github Codes/figure_style.R")

## The summary is written above; recomputed here only if that object is not
## in scope, so this section can also be run on its own.
if (!exists("sum_mis")) sum_mis <- summarise_mis(res_mis)

## S4: weak instruments -- the decisive figure, instrument strength swept
## from the regime observed in the real analysis to that of S1-S3.
p_s4 <- plot_s4_weak(sum_mis)

## S5: correlated instruments.
p_s5 <- plot_s5_correlated(sum_mis)

## S6: invalid / pleiotropic instruments, faceted by invalid fraction.
p_s6 <- plot_s6_pleiotropy(sum_mis)

## S7: nonlinear ligand effect.
p_s7 <- plot_s7_nonlinear(sum_mis)

## Estimate-level boxplots at two representative settings per scenario
## (see mis_representative() in figure_functions.R), in the same visual
## language as the S1-S3 figures. Built from the raw per-replicate results,
## so regenerate_figures.R can rebuild them from the saved CSVs.
p_box_score  <- plot_mis_estimate_box(res_mis, "score")
p_box_betaX  <- plot_mis_estimate_box(res_mis, "beta_X")
p_box_betaXZ <- plot_mis_estimate_box(res_mis, "beta_XZ")

save_fig(p_s4, "sim_S4_weak",       15, 10)
save_fig(p_s5, "sim_S5_correlated", 15, 10)
save_fig(p_s6, "sim_S6_pleiotropy", 18, 10)
save_fig(p_s7, "sim_S7_nonlinear",  14, 10)
save_fig(p_box_score,  "sim_S4_S7_score_box",  16, 11)
save_fig(p_box_betaX,  "sim_S4_S7_betaX_box",  16, 11)
save_fig(p_box_betaXZ, "sim_S4_S7_betaXZ_box", 16, 11)

## Package versions used for this run. The study has no reduced-run
## override, so the file name carries no "smoke_" prefix.
writeLines(utils::capture.output(sessionInfo()),
           "Results/sessionInfo_simulation_misspecification.txt")
