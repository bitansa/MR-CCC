############################################################
## posterior_sign_probabilities.R
##
## Refit the declared discoveries and report two quantities
## that the primary output cannot provide, because it stores
## summaries rather than draws:
##
##   1. POSTERIOR SIGN PROBABILITIES
##      P(beta_XZ > 0) and P(beta_XZ < 0).
##
##      A 95% credible interval answers only "does it exclude
##      zero", which collapses everything from a near-certain
##      sign to a coin flip into one bit. A posterior
##      probability of 0.92 and one of 0.51 both "include
##      zero" and are not remotely the same statement. For a
##      Bayesian method the probability is the natural report.
##
##   2. INTERVALS CONDITIONAL ON gamma = 1
##      The marginal posterior of beta_XZ is a MIXTURE: the
##      slab when the interaction is included, and a tight
##      spike near zero when it is not. The marginal interval
##      therefore mixes "how big is the effect" with "is
##      there an effect at all". Conditioning on gamma = 1
##      separates them and answers the magnitude question on
##      its own.
##
##      Both are reported. The conditional interval is NOT a
##      replacement for the marginal one -- it answers a
##      different question, and presenting only the tighter
##      of the two would misrepresent the evidence.
##
## PREREQUISITES (run in this order, from the project root):
##   Rcpp::sourceCpp("Github Codes/mr_ccc_gibbs.cpp")
##   source("Github Codes/build_lr_database.R")
##   source("Github Codes/real_data_analysis.R")
##   source("Github Codes/posterior_sign_probabilities.R")
##
## COST. Only the declared discoveries are refitted, and the
## sampler itself is fast (~14 s per chain at n = 651); the
## time goes on rebuilding pathway scores, roughly 2-3
## minutes per triplet, so a few minutes per discovery.
##
## OUTPUT (written to Results/):
##   posterior_sign_probabilities_<Cell1>_<Cell2>.csv
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})

## ---- Inherited objects ---------------------------------------------------
.needed <- c("Output_MR_CCC", "lr_filtered", "key_cols",
             "RNA.count.adj_Cell1", "RNA.count.adj_Cell2", "donor", "V_donor", "centre_cols", "orient_pc1",
             "get_SNP_matrix", "gsva_scores", "mr_ccc_gibbs",
             "N_ITER", "BURN_IN", "THIN", "N_CHAINS", "INIT_SCALE",
             "pip_thresh", "Cell1", "Cell2")
.missing <- .needed[!vapply(.needed, exists, logical(1))]
if (length(.missing) > 0) {
  stop("Objects missing -- source real_data_analysis.R first.\n  Missing: ",
       paste(.missing, collapse = ", "))
}

source("Github Codes/build_triplet_inputs.R")

## Which triplets to refit. Default: the declared discoveries.
TARGETS <- Output_MR_CCC %>%
  filter(gamma_mean > pip_thresh) %>%
  arrange(desc(gamma_mean))

cat("\nRefitting", nrow(TARGETS), "triplets with", N_CHAINS,
    "dispersed chains each.\n")

rows <- vector("list", nrow(TARGETS))

for (i in seq_len(nrow(TARGETS))) {

  lab <- paste0(TARGETS$ligand_col_name[i], "-", TARGETS$receptor_col_name[i])
  cat("  [", i, "/", nrow(TARGETS), "] ", lab, "\n", sep = "")

  inp <- build_triplet_inputs(TARGETS$ligand_col_name[i],
                              TARGETS$receptor_col_name[i],
                              TARGETS$pathway_name[i])
  if (is.null(inp)) {
    warning("Could not rebuild ", lab, " -- skipped."); next
  }

  n <- nrow(inp$X); g_scale <- min(n, 100.0)

  ## Identical sampler settings to the primary analysis.
  fits <- lapply(seq_len(N_CHAINS), function(cc) {
    mr_ccc_gibbs(
      inp$X, inp$Z, inp$Y, inp$G, inp$H, inp$V,
      n_iter = N_ITER, burn_in = BURN_IN, thin = THIN,
      a_sigma = 3.0, b_sigma = 2.0, a_rho = 3.0, b_rho = 1.0,
      nu1 = 1e-4, gG = g_scale, gH = g_scale, gV = g_scale,
      gZ = g_scale, gBeta = g_scale, ridge = 1e-8,
      init_gamma = if (cc %% 2L == 0L) 0L else 1L,
      init_scale = INIT_SCALE
    )
  })

  ## Pooled draws, rescaled to the standardized scale used throughout.
  sxz_y  <- inp$sd_X * inp$sd_Z / inp$sd_Y
  bXZ    <- unlist(lapply(fits, function(f) as.numeric(f$Beta_XZ_draws)),
                   use.names = FALSE) * sxz_y
  g_dr   <- unlist(lapply(fits, function(f) as.numeric(f$gamma_draws)),
                   use.names = FALSE)

  ## Draws in which the interaction is INCLUDED. When gamma = 0 the
  ## coefficient is drawn from the spike (variance nu1 = 1e-4), so those
  ## draws describe absence, not magnitude.
  inc <- g_dr == 1
  bXZ_cond <- bXZ[inc]

  q_marg <- stats::quantile(bXZ, c(0.025, 0.5, 0.975), names = FALSE)
  q_cond <- if (sum(inc) >= 100L) {
    stats::quantile(bXZ_cond, c(0.025, 0.5, 0.975), names = FALSE)
  } else rep(NA_real_, 3)

  rows[[i]] <- tibble(
    ligand   = TARGETS$ligand_col_name[i],
    receptor = TARGETS$receptor_col_name[i],
    pathway  = TARGETS$pathway_name[i],
    representation = inp$representation,
    PIP      = TARGETS$gamma_mean[i],

    # --- posterior sign probabilities (marginal) ---
    P_XZ_gt0 = mean(bXZ > 0),
    P_XZ_lt0 = mean(bXZ < 0),
    # The larger of the two: "the sign is X with probability P".
    P_sign   = max(mean(bXZ > 0), mean(bXZ < 0)),

    # --- sign probability CONDITIONAL on inclusion ---
    P_XZ_gt0_cond = if (sum(inc) >= 100L) mean(bXZ_cond > 0) else NA_real_,
    P_XZ_lt0_cond = if (sum(inc) >= 100L) mean(bXZ_cond < 0) else NA_real_,

    # --- intervals ---
    XZ_marg_lo = q_marg[1], XZ_marg_med = q_marg[2], XZ_marg_hi = q_marg[3],
    XZ_cond_lo = q_cond[1], XZ_cond_med = q_cond[2], XZ_cond_hi = q_cond[3],

    n_draws_total    = length(bXZ),
    n_draws_included = sum(inc)
  )
}

sign_tab <- bind_rows(rows[!vapply(rows, is.null, logical(1))])

dir.create("Results", showWarnings = FALSE)
out_csv <- file.path("Results",
             paste0("posterior_sign_probabilities_", Cell1, "_", Cell2, ".csv"))
write.csv(sign_tab, out_csv, row.names = FALSE)

cat("\n=========== POSTERIOR SIGN PROBABILITIES ===========\n")
print(as.data.frame(
  sign_tab %>%
    select(ligand, receptor, PIP, P_sign, P_XZ_lt0,
           XZ_marg_lo, XZ_marg_hi, XZ_cond_lo, XZ_cond_hi)),
  digits = 3)

cat("\nRead as: P_sign is the posterior probability that beta_XZ has the sign",
    "\nof its point estimate. XZ_marg_* mixes spike and slab; XZ_cond_* is",
    "\nconditional on the interaction being included in the model.\n")
cat("\nWritten to", out_csv, "\n")
