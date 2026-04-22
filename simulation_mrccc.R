############################################################
## simulation_mrccc.R
##
## MR-CCC SIMULATION FRAMEWORK
## 3 scenarios x 4 methods x 4 sample sizes
##
## USAGE:
##   1. Compile the Gibbs sampler:
##        Rcpp::sourceCpp("mr_ccc_gibbs.cpp")
##   2. Source this file:
##        source("simulation_mrccc.R")
##
## Methods compared:
##   OLS    -- naive OLS with interaction; no IV correction
##   MVMR   -- Multivariable MR (frequentist 2SLS)
##   MR-BMA -- Bayesian Model Averaging MR (Zuber et al., 2020)
##   MR-CCC -- Bayesian 2SLS with spike-and-slab prior (this paper)
##
## Scoring convention:
##   OLS   : score = 1 - p_t  (Wald t-test, H0: beta_X = 0;
##            reject if p_t <= alpha_sig)
##   MVMR  : score = 1 - p_t  (Wald t-test, H0: beta_X = 0;
##            no interaction term; reject if p_t <= alpha_sig)
##   MR-BMA: score = MIP_X    (marginal inclusion probability
##            for ligand X; reject if MIP_X > pip_thr)
##   MR-CCC: score = PIP      (posterior inclusion probability
##            P(gamma=1|data); reject if PIP >= pip_thr)
############################################################
# rm(list = ls())
# gc()
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)


############################################################
## MR-BMA: SELF-CONTAINED IMPLEMENTATION
##
## MR.BMA is not available on CRAN for newer R versions, so
## we implement the core algorithm directly. The method is
## described in Zuber et al. (2020, PLOS Genetics).
##
## Model:
##   For K = 2 candidate exposures (X, Z), enumerate all
##   2^K = 4 models M_k (subsets of active exposures).
##   Standardise data: y_std = betaY / se.betaY,
##                     W_std = betaX / se.betaY (row-wise)
##   so the outcome has unit error variance.
##
## G-prior (Zellner, g = J instruments by default):
##   log BF_{k,0} = -p_k/2 * log(1+g)
##                  + g/(2(1+g)) * y_std' W_k (W_k'W_k)^{-1} W_k' y_std
##
## Prior on models:  P(M_k) proportional to pi^{p_k} (1-pi)^{K - p_k}
##
## Outputs:
##   MIP  -- marginal inclusion probability per exposure
##   MACE -- model-averaged causal effect per exposure
##           (weighted average of within-model MAP estimates)
############################################################
mr_bma_simple <- function(betaX, betaY, se.betaX, se.betaY,
                          g = NULL, prior_pi = 0.1) {
  J <- nrow(betaX)
  K <- ncol(betaX)
  if (is.null(g)) g <- J  # Zellner's unit-information g = J
  
  # Standardise by se.betaY (makes outcome variance = 1)
  y_std <- betaY / se.betaY
  W_std <- sweep(betaX, 1, se.betaY, "/")  # J x K
  
  # Enumerate 2^K models (binary inclusion matrix)
  model_mat <- as.matrix(expand.grid(rep(list(0:1), K)))
  colnames(model_mat) <- colnames(betaX)
  n_models <- nrow(model_mat)
  
  log_BF     <- numeric(n_models)
  log_prior  <- numeric(n_models)
  theta_list <- vector("list", n_models)
  
  for (m in seq_len(n_models)) {
    inc <- which(model_mat[m, ] == 1)
    p_k <- length(inc)
    log_prior[m] <- p_k * log(prior_pi / (1 - prior_pi))
    
    if (p_k == 0L) {
      log_BF[m]       <- 0
      theta_list[[m]] <- numeric(K)
    } else {
      W_k <- W_std[, inc, drop = FALSE]
      XtX <- crossprod(W_k)
      Xty <- crossprod(W_k, y_std)
      
      XtX_inv <- tryCatch(solve(XtX), error = function(e) NULL)
      if (is.null(XtX_inv)) {
        log_BF[m]       <- -Inf
        theta_list[[m]] <- numeric(K)
        next
      }
      
      SS_model_k <- as.numeric(t(Xty) %*% XtX_inv %*% Xty)
      log_BF[m]  <- -p_k / 2 * log(1 + g) +
        g / (2 * (1 + g)) * SS_model_k
      
      # MAP estimate under g-prior (shrinkage factor g/(1+g))
      th_full      <- numeric(K)
      th_full[inc] <- as.numeric(g / (1 + g) * XtX_inv %*% Xty)
      theta_list[[m]] <- th_full
    }
  }
  
  # Posterior model probabilities (numerically stable)
  log_unnorm <- log_BF + log_prior
  log_unnorm <- log_unnorm - max(log_unnorm)
  PP <- exp(log_unnorm) / sum(exp(log_unnorm))
  
  # Marginal inclusion probabilities and model-averaged effects
  MIP  <- numeric(K)
  MACE <- numeric(K)
  for (k in seq_len(K)) {
    idx_inc <- which(model_mat[, k] == 1L)
    MIP[k]  <- sum(PP[idx_inc])
    MACE[k] <- sum(PP[idx_inc] *
                     vapply(idx_inc,
                            function(m) theta_list[[m]][k],
                            numeric(1L)))
  }
  
  names(MIP)  <- colnames(betaX)
  names(MACE) <- colnames(betaX)
  
  list(MIP = MIP, MACE = MACE, PP = PP,
       model = model_mat, BF = exp(log_BF))
}


############################################################
## USER SETTINGS
############################################################
sample_sizes  <- c(500, 1000, 10000, 30000)  # donor sample sizes to evaluate
R_reps        <- 20                            # Number of repetitions
conf_strength <- 0.7                           # confounder loading on X, Z, Y
pip_thr       <- 0.5                           # PIP / MIP threshold for rejection
alpha_sig     <- 0.05                          # frequentist significance level
pG_val        <- 5                             # number of sender cis-eQTL instruments
pH_val        <- 5                             # number of receiver cis-eQTL instruments
pV_val        <- 3                             # number of shared covariates
pi_val        <- 0.5                           # true first-stage eQTL effect size

method_levels <- c("OLS", "MVMR", "MR-BMA", "MR-CCC")

# Simulation scenarios:
#   S1 -- null (no causal communication; beta_X = beta_XZ = 0)
#   S2 -- signal with receptor-modulated interaction (beta_X = beta_XZ = 0.3)
#   S3 -- signal without interaction (beta_X = 0.3, beta_XZ = 0)
scenario_list <- list(
  S1 = list(beta_X = 0.0, beta_XZ = 0.0, beta_Z = 0.5),
  S2 = list(beta_X = 0.3, beta_XZ = 0.3, beta_Z = 0.5),
  S3 = list(beta_X = 0.3, beta_XZ = 0.0, beta_Z = 0.5)
)


############################################################
## DATA GENERATION
##
## Structural model:
##   X_i = G_i' pi_X + V_i' alpha_X + conf_strength * U_i + eps_Xi
##   Z_i = H_i' pi_Z + V_i' alpha_Z + conf_strength * U_i + eps_Zi
##   Y_i = beta_X * X_i + beta_Z * Z_i + beta_XZ * X_i * Z_i
##         + V_i' alpha_Y + conf_strength * U_i + eps_Yi
##
## U_i is a scalar unmeasured confounder shared by X, Z, Y.
## No intercept in the DGM for Y -- the working model intercept
## mu emerges only after IV substitution (Proposition 1 of the paper).
## All components are zero-mean by construction.
############################################################
generate_data <- function(n             = 300,
                          pG            = 5,
                          pH            = 5,
                          pV            = 3,
                          beta_X        = 0,
                          beta_XZ       = 0,
                          beta_Z        = 0.5,
                          conf_strength = 0.7,
                          pi_val        = 0.5,
                          seed          = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  G <- matrix(rnorm(n * pG), n, pG)  # sender cis-eQTL genotypes
  H <- matrix(rnorm(n * pH), n, pH)  # receiver cis-eQTL genotypes
  V <- matrix(rnorm(n * pV), n, pV)  # shared observed covariates
  U <- rnorm(n)                       # shared unmeasured confounder
  
  pi_X    <- rep(pi_val, pG)  # true first-stage eQTL effects for sender
  pi_Z    <- rep(pi_val, pH)  # true first-stage eQTL effects for receiver
  alpha_X <- rep(0.3, pV)     # covariate effects on ligand
  alpha_Z <- rep(0.3, pV)     # covariate effects on receptor
  alpha_Y <- rep(0.3, pV)     # covariate effects on pathway
  
  X <- as.numeric(G %*% pi_X + V %*% alpha_X + conf_strength * U + rnorm(n))
  Z <- as.numeric(H %*% pi_Z + V %*% alpha_Z + conf_strength * U + rnorm(n))
  Y <- as.numeric(beta_X  * X +
                    beta_Z  * Z +
                    beta_XZ * (X * Z) +
                    V %*% alpha_Y +
                    conf_strength * U +
                    rnorm(n))
  
  list(
    X = matrix(X, ncol = 1),
    Z = matrix(Z, ncol = 1),
    Y = matrix(Y, ncol = 1),
    G = G, H = H, V = V
  )
}


############################################################
## HELPERS
############################################################

# Safely extract a named coefficient from a fitted lm object.
# Returns NA if the term is absent (e.g., dropped due to collinearity).
safe_coef <- function(fit, term) {
  cf <- coef(fit)
  if (!term %in% names(cf)) return(NA_real_)
  unname(cf[term])
}

# Format mean (SD) of a numeric vector for summary tables.
# Finite values only; returns "---" if nothing is finite.
fmt_mean_sd <- function(x, digits = 3) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("---")
  paste0(
    format(round(mean(x), digits), nsmall = digits), " (",
    format(round(sd(x),   digits), nsmall = digits), ")"
  )
}


############################################################
## RUN ALL 4 METHODS ON ONE DATASET
############################################################
run_methods <- function(dat, pip_thr = 0.5, alpha_sig = 0.05) {
  
  X <- dat$X;  Z <- dat$Z;  Y <- dat$Y
  G <- dat$G;  H <- dat$H;  V <- dat$V
  n  <- nrow(X)
  pV <- ncol(V)
  
  # ----------------------------------------------------------
  # CENTERING: subtract sample means from x, z, y.
  # Enforces the paper's centering assumption exactly in each
  # sample. Slopes beta_X and beta_XZ are invariant to this
  # location shift; only the intercept mu changes.
  # ----------------------------------------------------------
  x <- as.numeric(X) - mean(X)
  z <- as.numeric(Z) - mean(Z)
  y <- as.numeric(Y) - mean(Y)
  
  V_names <- paste0("V", seq_len(pV))
  df <- data.frame(
    Y = y, X = x, Z = z,
    setNames(as.data.frame(V), V_names)
  )
  
  all_IVs  <- cbind(G, H)         # combined instrument matrix [G | H]
  n_all_IV <- ncol(all_IVs)       # pG + pH
  IV_names <- paste0("IV", seq_len(n_all_IV))
  
  ##########################################################
  ## (1) OLS with interaction
  ##
  ## Model:  Y_i = mu + beta_X X_i + beta_Z Z_i +
  ##               beta_XZ X_i Z_i + V_i' alpha_Y + eps_i
  ## Fit via OLS (no IV correction; confounding ignored).
  ##
  ## Communication score: 1 - p_F
  ##   Joint Wald F-test  H0: beta_X = beta_XZ = 0
  ##   Contrast matrix C selects (beta_X, beta_XZ) from coef.
  ##   F  = (1/2)(C beta)' [C vcov C']^{-1} (C beta)
  ##   p_F = P(F(2, df.res) > F)
  ##   score = 1 - p_F;   reject if p_F <= alpha_sig
  ##########################################################
  ols_fit     <- lm(Y ~ X * Z + ., data = df)
  OLS_beta_X  <- safe_coef(ols_fit, "X")
  OLS_beta_XZ <- safe_coef(ols_fit, "X:Z")
  
  cf_names <- names(coef(ols_fit))
  idx_X    <- which(cf_names == "X")
  idx_XZ   <- which(cf_names == "X:Z")
  
  beta_ols  <- coef(ols_fit)
  V_ols     <- vcov(ols_fit)
  C_ols     <- matrix(0.0, 2L, length(beta_ols))
  C_ols[1L, idx_X]  <- 1.0
  C_ols[2L, idx_XZ] <- 1.0
  
  Cb_ols  <- C_ols %*% beta_ols
  CVC_ols <- C_ols %*% V_ols %*% t(C_ols)
  F_ols   <- as.numeric(0.5 * t(Cb_ols) %*% solve(CVC_ols) %*% Cb_ols)
  p_ols   <- pf(F_ols, 2L, df.residual(ols_fit), lower.tail = FALSE)
  
  OLS_score  <- 1 - p_ols
  OLS_reject <- as.integer(p_ols <= alpha_sig)
  
  ##########################################################
  ## (2) MVMR  (Multivariable Mendelian Randomization)
  ##
  ## Standard MVMR is an additive model: it estimates the
  ## causal effect of each exposure adjusting for the other,
  ## but does NOT include an interaction term. The interaction
  ## beta_XZ is therefore not estimable from this method.
  ##
  ## All instruments [G, H] are used jointly in the first
  ## stage for BOTH exposures X and Z (standard MVMR practice).
  ##
  ## First stage (OLS per exposure):
  ##   X_i = mu_X + [G_i, H_i]' pi_X + V_i' alpha_X + eps_Xi
  ##   Z_i = mu_Z + [G_i, H_i]' pi_Z + V_i' alpha_Z + eps_Zi
  ##
  ## Second-stage design matrix (additive, no interaction):
  ##   W = [1_n | Xhat_c | Zhat_c | V]
  ##   theta = (W'W)^{-1} W'y;   theta[2] = beta_X
  ##
  ## Communication score: 1 - p_t
  ##   Wald t-test  H0: beta_X = 0
  ##   t = beta_X_hat / se(beta_X_hat)
  ##   p_t = 2 * P(t(n - p) > |t|)
  ##   score = 1 - p_t;   reject if p_t <= alpha_sig
  ##########################################################
  df_fsX_mv <- data.frame(
    X = x,
    setNames(as.data.frame(all_IVs), IV_names),
    setNames(as.data.frame(V),       V_names)
  )
  df_fsZ_mv <- data.frame(
    Z = z,
    setNames(as.data.frame(all_IVs), IV_names),
    setNames(as.data.frame(V),       V_names)
  )
  
  fs_X_mv   <- lm(X ~ ., data = df_fsX_mv)
  fs_Z_mv   <- lm(Z ~ ., data = df_fsZ_mv)
  Xhat_mv_c <- fitted(fs_X_mv) - mean(fitted(fs_X_mv))
  Zhat_mv_c <- fitted(fs_Z_mv) - mean(fitted(fs_Z_mv))
  
  # Additive second stage: no interaction term
  p_mv     <- 3L + pV   # intercept + Xhat + Zhat + V columns
  W_mv     <- cbind(1, Xhat_mv_c, Zhat_mv_c, V)
  WtW_mv   <- crossprod(W_mv)
  WtW_mv_i <- solve(WtW_mv)
  theta_mv <- WtW_mv_i %*% crossprod(W_mv, y)
  e_mv     <- y - W_mv %*% theta_mv
  s2_mv    <- sum(e_mv^2) / (n - p_mv)
  
  MVMR_beta_X  <- theta_mv[2L, 1L]
  MVMR_beta_XZ <- NA_real_   # interaction not modeled in standard MVMR
  
  # 1-df Wald t-test: H0: beta_X = 0
  se_betaX_mv <- sqrt(s2_mv * WtW_mv_i[2L, 2L])
  t_mv        <- MVMR_beta_X / se_betaX_mv
  p_mv_F      <- 2 * pt(abs(t_mv), df = n - p_mv, lower.tail = FALSE)
  
  MVMR_score  <- 1 - p_mv_F
  MVMR_reject <- as.integer(p_mv_F <= alpha_sig)
  
  ##########################################################
  ## (3) MR-BMA  (Bayesian Model Averaging MR)
  ##
  ## Candidate exposures: X (ligand) and Z (receptor).
  ## Instruments: all SNPs in [G, H].
  ##
  ## For each instrument g_j, one-sample summary statistics
  ## are obtained via simple OLS (intercept included):
  ##   beta_{X,j}  <- coef of g_j in: lm(x ~ g_j)
  ##   beta_{Z,j}  <- coef of g_j in: lm(z ~ g_j)
  ##   beta_{Y,j}  <- coef of g_j in: lm(y ~ g_j)
  ## with corresponding standard errors.
  ##
  ## mr_bma_simple() evaluates all 2^2 = 4 exposure models
  ## {null, {X}, {Z}, {X,Z}} via a Bayesian g-prior and returns
  ##   MIP  -- marginal inclusion probability per exposure
  ##   MACE -- model-averaged causal effect per exposure
  ##
  ## NOTE: MR-BMA does not model the ligand-receptor
  ## interaction beta_XZ; that parameter is not estimable
  ## from this method (beta_XZ_hat is set to NA).
  ##
  ## Communication score: MIP_X
  ##   = posterior probability that ligand X has a nonzero
  ##     causal effect on the outcome Y.
  ##   reject if MIP_X > pip_thr
  ##########################################################
  betaX_bma    <- matrix(NA_real_, n_all_IV, 2L,
                         dimnames = list(NULL, c("X", "Z")))
  se_betaX_bma <- matrix(NA_real_, n_all_IV, 2L,
                         dimnames = list(NULL, c("X", "Z")))
  betaY_bma    <- numeric(n_all_IV)
  se_betaY_bma <- numeric(n_all_IV)
  
  for (j in seq_len(n_all_IV)) {
    iv_j <- all_IVs[, j]
    
    sXj <- summary(lm(x ~ iv_j))$coefficients
    betaX_bma[j, 1]    <- sXj[2L, 1L]
    se_betaX_bma[j, 1] <- sXj[2L, 2L]
    
    sZj <- summary(lm(z ~ iv_j))$coefficients
    betaX_bma[j, 2]    <- sZj[2L, 1L]
    se_betaX_bma[j, 2] <- sZj[2L, 2L]
    
    sYj <- summary(lm(y ~ iv_j))$coefficients
    betaY_bma[j]    <- sYj[2L, 1L]
    se_betaY_bma[j] <- sYj[2L, 2L]
  }
  
  bma_res <- tryCatch(
    mr_bma_simple(betaX    = betaX_bma,
                  betaY    = betaY_bma,
                  se.betaX = se_betaX_bma,
                  se.betaY = se_betaY_bma),
    error = function(e) NULL
  )
  
  if (is.null(bma_res)) {
    MRBMA_score   <- NA_real_
    MRBMA_reject  <- NA_integer_
    MRBMA_beta_X  <- NA_real_
    MRBMA_beta_XZ <- NA_real_
  } else {
    MRBMA_score   <- as.numeric(bma_res$MIP["X"])
    MRBMA_reject  <- as.integer(MRBMA_score > pip_thr)
    MRBMA_beta_X  <- as.numeric(bma_res$MACE["X"])
    MRBMA_beta_XZ <- NA_real_  # interaction not modeled by MR-BMA
  }
  
  ##########################################################
  ## (4) MR-CCC  (Bayesian 2SLS with spike-and-slab prior)
  ##
  ## Separate first-stage models: X ~ G + V;  Z ~ H + V
  ## (each exposure instrumented only by its own gene's IVs,
  ##  consistent with the structural CCC model).
  ##
  ## Spike-and-slab prior on beta = (beta_X, beta_XZ)':
  ##   gamma = 1  =>  beta ~ N(0, gBeta * sigma_Y^2 * (W'W)^{-1})
  ##   gamma = 0  =>  beta ~ N(0, nu1 * gBeta * sigma_Y^2 * (W'W)^{-1})
  ##                  (near-zero spike; nu1 = 1e-4)
  ##
  ## g-prior scales follow the paper's recommendation:
  ##   gG = gH = gV = gZ = gBeta = min(n, 100)
  ##
  ## Communication score: P(gamma=1 | data)  (PIP)
  ##   Estimated as the posterior mean of gamma across
  ##   Gibbs iterations (after burn-in, with thinning).
  ##   reject if PIP >= pip_thr
  ##########################################################
  res_mr <- mr_ccc_gibbs(
    matrix(x, ncol = 1),
    matrix(z, ncol = 1),
    matrix(y, ncol = 1),
    G, H, V,
    n_iter  = 20000,
    burn_in = 2000,
    thin    = 5,
    gG      = min(n, 100),
    gH      = min(n, 100),
    gV      = min(n, 100),
    gZ      = min(n, 100),
    gBeta   = min(n, 100)
  )
  
  MR_beta_X  <- res_mr$Beta_X_mean
  MR_beta_XZ <- res_mr$Beta_XZ_mean
  MR_score   <- res_mr$gamma_mean
  
  list(
    OLS_beta_X  = OLS_beta_X,   OLS_beta_XZ  = OLS_beta_XZ,
    OLS_score   = OLS_score,    OLS_reject   = OLS_reject,
    
    MVMR_beta_X  = MVMR_beta_X,  MVMR_beta_XZ  = MVMR_beta_XZ,
    MVMR_score   = MVMR_score,   MVMR_reject   = MVMR_reject,
    
    MRBMA_beta_X  = MRBMA_beta_X,  MRBMA_beta_XZ  = MRBMA_beta_XZ,
    MRBMA_score   = MRBMA_score,   MRBMA_reject   = MRBMA_reject,
    
    MR_beta_X  = MR_beta_X,   MR_beta_XZ  = MR_beta_XZ,
    MR_score   = MR_score
  )
}


############################################################
## ONE REPLICATE
##
## Generates one dataset for a given (scenario, n, seed)
## combination, runs all four methods, and returns a
## one-row data frame with all results.
##
## NOTE: scenario_list is accessed from the global environment.
############################################################
run_one_replicate <- function(scenario_name,
                              n,
                              seed,
                              pG            = 5,
                              pH            = 5,
                              pV            = 3,
                              conf_strength = 0.7,
                              pi_val        = 0.5,
                              pip_thr       = 0.5,
                              alpha_sig     = 0.05) {
  pars <- scenario_list[[scenario_name]]
  
  dat <- generate_data(
    n             = n,
    pG            = pG,
    pH            = pH,
    pV            = pV,
    beta_X        = pars$beta_X,
    beta_XZ       = pars$beta_XZ,
    beta_Z        = pars$beta_Z,
    conf_strength = conf_strength,
    pi_val        = pi_val,
    seed          = seed
  )
  
  fit <- run_methods(dat, pip_thr = pip_thr, alpha_sig = alpha_sig)
  
  # gamma_true = 1 if at least one causal effect is nonzero
  data.frame(
    scenario     = scenario_name,
    n            = n,
    gamma_true   = as.integer((pars$beta_X != 0) || (pars$beta_XZ != 0)),
    beta_X_true  = pars$beta_X,
    beta_XZ_true = pars$beta_XZ,
    
    OLS_score   = fit$OLS_score,    OLS_reject  = fit$OLS_reject,
    OLS_beta_X  = fit$OLS_beta_X,   OLS_beta_XZ = fit$OLS_beta_XZ,
    
    MVMR_score   = fit$MVMR_score,   MVMR_reject  = fit$MVMR_reject,
    MVMR_beta_X  = fit$MVMR_beta_X,  MVMR_beta_XZ = fit$MVMR_beta_XZ,
    
    MRBMA_score   = fit$MRBMA_score,   MRBMA_reject  = fit$MRBMA_reject,
    MRBMA_beta_X  = fit$MRBMA_beta_X,  MRBMA_beta_XZ = fit$MRBMA_beta_XZ,
    
    MR_score   = fit$MR_score,
    MR_reject  = as.integer(fit$MR_score >= pip_thr),
    MR_beta_X  = fit$MR_beta_X,
    MR_beta_XZ = fit$MR_beta_XZ
  )
}


############################################################
## RUN FULL SIMULATION GRID
##
## Builds the complete (scenario x n x replicate) job table,
## then dispatches each job in parallel via mclapply.
##
## NOTE: mclapply uses forking and is not supported on Windows.
## On Windows, replace with:
##   lapply(seq_len(nrow(jobs)), function(k) { ... })
## or use parallel::parLapply with a PSOCK cluster instead.
##
## OMP/BLAS thread counts are set to 1 to prevent nested
## parallelism conflicts with mclapply.
############################################################
run_simulation_grid <- function(sample_sizes  = c(500, 1000, 10000, 30000),
                                R_reps        = 20,
                                pG            = 5,
                                pH            = 5,
                                pV            = 3,
                                conf_strength = 0.7,
                                pi_val        = 0.5,
                                pip_thr       = 0.5,
                                alpha_sig     = 0.05,
                                mc_cores      = max(1, detectCores() - 1)) {
  # Build job table: one row per (scenario, n, replicate)
  # Sorted by descending n so large jobs start first (better load balancing)
  jobs <- expand.grid(
    n        = sample_sizes,
    scenario = names(scenario_list),
    rep      = seq_len(R_reps),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(n)) %>%
    mutate(seed = seq_len(n()))
  
  # Prevent BLAS/OpenMP from spawning threads inside each forked worker
  Sys.setenv(OMP_NUM_THREADS      = "1",
             OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS      = "1")
  
  res_list <- mclapply(seq_len(nrow(jobs)), function(k) {
    run_one_replicate(
      scenario_name = jobs$scenario[k],
      n             = jobs$n[k],
      seed          = jobs$seed[k],
      pG            = pG,
      pH            = pH,
      pV            = pV,
      conf_strength = conf_strength,
      pi_val        = pi_val,
      pip_thr       = pip_thr,
      alpha_sig     = alpha_sig
    )
  }, mc.cores = mc_cores)
  
  bind_rows(res_list)
}


############################################################
## TIDY LONG FORMAT
##
## Reshapes the wide results data frame into long format
## (one row per method per replicate) for plotting and
## summarisation. pip_thr is used to re-derive MR-CCC
## rejection here to keep the long format self-consistent.
############################################################
make_long_results <- function(df_raw, pip_thr = 0.5) {
  bind_rows(
    df_raw %>% transmute(
      scenario, n, gamma_true, beta_X_true, beta_XZ_true,
      method      = "OLS",
      score       = OLS_score,
      reject      = OLS_reject,
      beta_X_hat  = OLS_beta_X,
      beta_XZ_hat = OLS_beta_XZ
    ),
    df_raw %>% transmute(
      scenario, n, gamma_true, beta_X_true, beta_XZ_true,
      method      = "MVMR",
      score       = MVMR_score,
      reject      = MVMR_reject,
      beta_X_hat  = MVMR_beta_X,
      beta_XZ_hat = MVMR_beta_XZ
    ),
    df_raw %>% transmute(
      scenario, n, gamma_true, beta_X_true, beta_XZ_true,
      method      = "MR-BMA",
      score       = MRBMA_score,
      reject      = MRBMA_reject,
      beta_X_hat  = MRBMA_beta_X,
      beta_XZ_hat = MRBMA_beta_XZ  # NA: interaction not modeled
    ),
    df_raw %>% transmute(
      scenario, n, gamma_true, beta_X_true, beta_XZ_true,
      method      = "MR-CCC",
      score       = MR_score,
      reject      = as.integer(MR_score >= pip_thr),
      beta_X_hat  = MR_beta_X,
      beta_XZ_hat = MR_beta_XZ
    )
  ) %>%
    mutate(method = factor(method, levels = method_levels))
}


############################################################
## TABLE SUMMARIES
############################################################

# Produce a summary table with mean (SD) for each metric
# grouped by (scenario, n, method).
make_table_summary <- function(df_long) {
  df_long %>%
    group_by(scenario, n, method) %>%
    summarise(
      Score          = fmt_mean_sd(score),
      Rejection_rate = fmt_mean_sd(reject),
      Bias_betaX     = fmt_mean_sd(beta_X_hat - beta_X_true),
      MAD_betaX      = fmt_mean_sd(abs(beta_X_hat - beta_X_true)),
      Bias_betaXZ    = fmt_mean_sd(beta_XZ_hat - beta_XZ_true),
      MAD_betaXZ     = fmt_mean_sd(abs(beta_XZ_hat - beta_XZ_true)),
      .groups = "drop"
    ) %>%
    arrange(scenario, n, method)
}

# Split summary table by scenario for convenient printing.
make_scenario_tables <- function(df_long) {
  tab <- make_table_summary(df_long)
  list(
    table_s1 = tab %>% filter(scenario == "S1") %>% arrange(n, method),
    table_s2 = tab %>% filter(scenario == "S2") %>% arrange(n, method),
    table_s3 = tab %>% filter(scenario == "S3") %>% arrange(n, method)
  )
}


############################################################
## RUN SIMULATION
############################################################
cat("Running simulation (4 methods, 3 scenarios, 4 sample sizes)...\n")

sim_raw  <- run_simulation_grid(
  sample_sizes  = sample_sizes,
  R_reps        = R_reps,
  pG            = pG_val,
  pH            = pH_val,
  pV            = pV_val,
  pi_val        = pi_val,
  conf_strength = conf_strength,
  pip_thr       = pip_thr,
  alpha_sig     = alpha_sig
)

sim_long <- make_long_results(sim_raw, pip_thr = pip_thr)

# Print scenario-level summary tables
tabs <- make_scenario_tables(sim_long)
cat("\n--- Scenario 1 (null) ---\n")
print(tabs$table_s1)
cat("\n--- Scenario 2 (signal + interaction) ---\n")
print(tabs$table_s2)
cat("\n--- Scenario 3 (signal, no interaction) ---\n")
print(tabs$table_s3)


############################################################
## PUBLICATION THEME
############################################################
theme_pub <- theme_bw(base_size = 14) +
  theme(
    axis.title    = element_text(size = 16, face = "bold"),
    axis.text     = element_text(size = 12),
    strip.text    = element_text(size = 13, face = "bold"),
    plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    legend.title  = element_text(size = 14, face = "bold"),
    legend.text   = element_text(size = 12),
    axis.text.x   = element_text(angle = 25, hjust = 1)
  )

# Reference lines for true parameter values
gamma_df  <- sim_long %>%
  distinct(scenario) %>%
  mutate(gamma_true = case_when(scenario == "S1" ~ 0, TRUE ~ 1))

betaX_df  <- sim_long %>% distinct(scenario, n, beta_X_true)
betaXZ_df <- sim_long %>% distinct(scenario, n, beta_XZ_true)


############################################################
## PLOT 1: COMMUNICATION SCORE
##
## Score definition per method:
##   OLS, MVMR : 1 - p_F  (joint Wald F-test)
##   MR-BMA    : MIP_X    (marginal inclusion probability for X)
##   MR-CCC    : PIP      (posterior inclusion probability P(gamma=1|data))
## Red dashed line: true communication state (0 for S1, 1 for S2/S3)
############################################################
p_score <- ggplot(sim_long, aes(x = method, y = score)) +
  geom_boxplot(outlier.alpha = 0.35) +
  geom_hline(
    data      = gamma_df,
    aes(yintercept = gamma_true),
    color     = "red",
    linetype  = "dashed",
    linewidth = 0.9
  ) +
  facet_grid(scenario ~ n, scales = "free_y") +
  theme_pub +
  labs(
    title    = "Communication scores across methods and scenarios",
    subtitle = "Rows = scenario; columns = sample size",
    x        = "Method",
    y        = "Communication score"
  )
print(p_score)


############################################################
## PLOT 2: beta_X ESTIMATES
##
## OLS / MVMR / MR-CCC: point estimate of main causal effect
## MR-BMA: model-averaged causal effect (MACE) for X
## Red dashed line: true beta_X value
############################################################
p_betaX <- ggplot(sim_long, aes(x = method, y = beta_X_hat)) +
  geom_boxplot(outlier.alpha = 0.35) +
  geom_hline(
    data      = betaX_df,
    aes(yintercept = beta_X_true),
    colour    = "red",
    linetype  = "dashed",
    linewidth = 0.9
  ) +
  facet_grid(scenario ~ n, scales = "free_y") +
  theme_pub +
  labs(
    title    = expression(paste("Estimated ", beta[X], " across methods and scenarios")),
    subtitle = "Red dashed line = true value",
    x        = "Method",
    y        = expression(hat(beta)[X])
  )
print(p_betaX)


############################################################
## PLOT 3: beta_XZ ESTIMATES
##
## Only OLS and MR-CCC estimate the interaction term.
## MR-BMA and MVMR are excluded (neither models beta_XZ).
## Red dashed line: true beta_XZ value.
############################################################
sim_long_xz <- sim_long %>%
  filter(method %in% c("OLS", "MR-CCC")) %>%
  mutate(method = factor(method, levels = c("OLS", "MR-CCC")))

betaXZ_df_xz <- sim_long_xz %>%
  distinct(scenario, n, beta_XZ_true)

p_betaXZ <- ggplot(sim_long_xz, aes(x = method, y = beta_XZ_hat)) +
  geom_boxplot(outlier.alpha = 0.35) +
  geom_hline(
    data      = betaXZ_df_xz,
    aes(yintercept = beta_XZ_true),
    colour    = "red",
    linetype  = "dashed",
    linewidth = 0.9
  ) +
  facet_grid(scenario ~ n, scales = "free_y") +
  theme_pub +
  labs(
    title    = expression(paste("Estimated ", beta[XZ], " across methods and scenarios")),
    subtitle = "Red dashed line = true value  (MR-BMA and MVMR excluded: interaction not modeled)",
    x        = "Method",
    y        = expression(hat(beta)[XZ])
  )
print(p_betaXZ)


############################################################
## SAVE (uncomment to write outputs)
############################################################
# write.csv(tabs$table_s1, "table_s1.csv", row.names = FALSE)
# write.csv(tabs$table_s2, "table_s2.csv", row.names = FALSE)
# write.csv(tabs$table_s3, "table_s3.csv", row.names = FALSE)
# ggsave("Plots/sim_score.pdf",  p_score,  width = 14, height = 9)
# ggsave("Plots/sim_betaX.pdf",  p_betaX,  width = 14, height = 9)
# ggsave("Plots/sim_betaXZ.pdf", p_betaXZ, width = 14, height = 9)
