############################################################
## leave_one_instrument_out.R
##
## LEAVE-ONE-INSTRUMENT-OUT SENSITIVITY FOR THE DECLARED DISCOVERIES
##
## Assesses whether any single instrument drives a finding. For
## each declared discovery, each ligand instrument and each
## receptor instrument is dropped in turn, the model is refit
## under exactly the primary sampler settings, and the PIP is
## recorded. A finding that depends on one SNP will show a
## large drop in PIP -- typically below the discovery
## threshold -- when that SNP is removed.
##
## WHAT IS REPORTED, per discovery:
##   PIP_full        the primary PIP (refit here, not copied,
##                   so every number in the table comes from
##                   the same code path)
##   PIP_min_dropL   smallest PIP over all single ligand-IV
##                   drops, and which SNP produced it
##   PIP_min_dropR   likewise for receptor IVs
##   n_drops_below   how many single drops push the PIP below
##                   the threshold (0 = robust to any one SNP)
##   PIP_range       full range of PIPs across all drops
##
## INTERPRETATION. This is a robustness check on the
## INSTRUMENT SET, not a test of the exclusion restriction.
## An instrument with a direct effect on the pathway would
## bias the estimate whether or not it dominates the PIP, so
## a clean leave-one-out result rules out "one bad SNP" but
## not "several mildly invalid SNPs". Simulation scenario S6
## addresses the latter.
##
## Note on weak instruments. With first-stage F around 1.5,
## the analysis measures how much each PIP depends on any
## single instrument. Its results should be read alongside
## the F values, not as evidence of instrument strength.
##
## PREREQUISITES (run in this order, from the project root):
##   Rcpp::sourceCpp("Github Codes/mr_ccc_gibbs.cpp")
##   source("Github Codes/build_lr_database.R")
##   source("Github Codes/real_data_analysis.R")
##   source("Github Codes/leave_one_instrument_out.R")
##
## COST. Each discovery needs (n_iv_L + n_iv_R + 1) refits, up
## to 21 at ten instruments per gene. Pathway scores are
## rebuilt ONCE per triplet and reused across its drops, so
## the cost is ~2-3 min of scoring plus ~56 s of sampling per
## refit at the primary configuration of four chains of
## 100,000 iterations: roughly 20-25 min per discovery.
##
## OUTPUT (written to Results/):
##   leave_one_out_summary_<Cell1>_<Cell2>.csv   one row per discovery
##   leave_one_out_detail_<Cell1>_<Cell2>.csv    one row per drop
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})

## ---- Inherited objects ---------------------------------------------------
.needed <- c("Output_MR_CCC", "lr_filtered", "key_cols",
             "RNA.count.adj_Cell1", "RNA.count.adj_Cell2", "donor", "V_donor", "centre_cols", "orient_pc1",
             "get_SNP_matrix", "gsva_scores", "mr_ccc_gibbs",
             "N_ITER", "BURN_IN", "THIN", "N_CHAINS", "INIT_SCALE",
             "pip_thresh", "Cell1", "Cell2",
             # Instrument selection is LD-clumped; a session holding an
             # older get_SNP_matrix() lacks these and must re-source.
             "LD_R2_MAX", "clump_by_ld")
.missing <- .needed[!vapply(.needed, exists, logical(1))]
if (length(.missing) > 0) {
  stop("Objects missing -- source real_data_analysis.R first.\n  Missing: ",
       paste(.missing, collapse = ", "))
}

source("Github Codes/build_triplet_inputs.R")

## Null-coalescing helper: instrument matrices from get_SNP_matrix() may or
## may not carry SNP identifiers as column names. Must be defined BEFORE the
## loop that uses it.
`%||%` <- function(a, b) if (is.null(a)) b else a

## Which triplets to test. Default: the declared discoveries.
TARGETS <- Output_MR_CCC %>%
  filter(gamma_mean > pip_thresh) %>%
  arrange(desc(gamma_mean))

## ---- One refit under the primary settings ---------------------------------
## Identical sampler call to real_data_analysis.R Step 7. Returns the PIP.
fit_pip <- function(inp, G, H) {
  n <- nrow(inp$X); g_scale <- min(n, 100.0)
  fits <- lapply(seq_len(N_CHAINS), function(cc) {
    mr_ccc_gibbs(
      inp$X, inp$Z, inp$Y, G, H, inp$V,
      n_iter = N_ITER, burn_in = BURN_IN, thin = THIN,
      a_sigma = 3.0, b_sigma = 2.0, a_rho = 3.0, b_rho = 1.0,
      nu1 = 1e-4, gG = g_scale, gH = g_scale, gV = g_scale,
      gZ = g_scale, gBeta = g_scale, ridge = 1e-8,
      init_gamma = if (cc %% 2L == 0L) 0L else 1L,
      init_scale = INIT_SCALE
    )
  })
  mean(unlist(lapply(fits, function(f) as.numeric(f$gamma_draws)),
              use.names = FALSE))
}

cat("\nLeave-one-instrument-out on", nrow(TARGETS), "discoveries.\n")

summary_rows <- vector("list", nrow(TARGETS))
detail_rows  <- list()

for (i in seq_len(nrow(TARGETS))) {

  lab <- paste0(TARGETS$ligand_col_name[i], "-", TARGETS$receptor_col_name[i])
  cat("  [", i, "/", nrow(TARGETS), "] ", lab, "\n", sep = "")

  ## Pathway scores are the expensive part; build them once per triplet.
  inp <- build_triplet_inputs(TARGETS$ligand_col_name[i],
                              TARGETS$receptor_col_name[i],
                              TARGETS$pathway_name[i])
  if (is.null(inp)) {
    warning("Could not rebuild ", lab, " -- skipped."); next
  }

  G <- inp$G; H <- inp$H
  nL <- ncol(G); nR <- ncol(H)

  ## Full-set refit, so PIP_full is on the same footing as the drops.
  pip_full <- fit_pip(inp, G, H)
  cat("      full set: PIP =", round(pip_full, 3),
      " (", nL, "L +", nR, "R instruments)\n")

  ## Drop each ligand instrument in turn. A gene with a single instrument
  ## cannot have it dropped (the model needs >= 1), so it is skipped and
  ## flagged; this happens for genes where selection returned only one SNP.
  drops <- list()
  if (nL >= 2L) {
    for (j in seq_len(nL)) {
      p <- fit_pip(inp, G[, -j, drop = FALSE], H)
      drops[[length(drops) + 1L]] <- tibble(
        ligand = TARGETS$ligand_col_name[i],
        receptor = TARGETS$receptor_col_name[i],
        side = "ligand", dropped = colnames(G)[j] %||% paste0("L", j),
        PIP = p)
    }
  }
  if (nR >= 2L) {
    for (j in seq_len(nR)) {
      p <- fit_pip(inp, G, H[, -j, drop = FALSE])
      drops[[length(drops) + 1L]] <- tibble(
        ligand = TARGETS$ligand_col_name[i],
        receptor = TARGETS$receptor_col_name[i],
        side = "receptor", dropped = colnames(H)[j] %||% paste0("R", j),
        PIP = p)
    }
  }
  dd <- bind_rows(drops)
  detail_rows[[i]] <- dd

  minL <- dd %>% filter(side == "ligand")   %>% slice_min(PIP, n = 1, with_ties = FALSE)
  minR <- dd %>% filter(side == "receptor") %>% slice_min(PIP, n = 1, with_ties = FALSE)

  summary_rows[[i]] <- tibble(
    ligand   = TARGETS$ligand_col_name[i],
    receptor = TARGETS$receptor_col_name[i],
    pathway  = TARGETS$pathway_name[i],
    PIP_primary_run = TARGETS$gamma_mean[i],   # from Output_MR_CCC
    PIP_full        = pip_full,                # refit here
    n_iv_L = nL, n_iv_R = nR,
    PIP_min_dropL   = if (nrow(minL)) minL$PIP else NA_real_,
    worst_dropL     = if (nrow(minL)) minL$dropped else NA_character_,
    PIP_min_dropR   = if (nrow(minR)) minR$PIP else NA_real_,
    worst_dropR     = if (nrow(minR)) minR$dropped else NA_character_,
    PIP_min_any     = if (nrow(dd)) min(dd$PIP) else NA_real_,
    PIP_max_any     = if (nrow(dd)) max(dd$PIP) else NA_real_,
    n_drops         = nrow(dd),
    n_drops_below   = sum(dd$PIP < pip_thresh),
    robust_to_any_single_drop = all(dd$PIP > pip_thresh)
  )

  cat("      min PIP over drops:", round(min(dd$PIP), 3),
      " | drops below", pip_thresh, ":", sum(dd$PIP < pip_thresh),
      "of", nrow(dd), "\n")
}

loo_summary <- bind_rows(summary_rows[!vapply(summary_rows, is.null, logical(1))])
loo_detail  <- bind_rows(detail_rows[!vapply(detail_rows, is.null, logical(1))])

dir.create("Results", showWarnings = FALSE)
f_sum <- file.path("Results", paste0("leave_one_out_summary_", Cell1, "_", Cell2, ".csv"))
f_det <- file.path("Results", paste0("leave_one_out_detail_",  Cell1, "_", Cell2, ".csv"))
write.csv(loo_summary, f_sum, row.names = FALSE)
write.csv(loo_detail,  f_det, row.names = FALSE)

cat("\n=========== LEAVE-ONE-INSTRUMENT-OUT ===========\n")
print(as.data.frame(
  loo_summary %>%
    select(ligand, receptor, PIP_full, PIP_min_any, PIP_max_any,
           n_drops, n_drops_below, robust_to_any_single_drop)),
  digits = 3)
cat("\nDiscoveries robust to every single-instrument drop:",
    sum(loo_summary$robust_to_any_single_drop), "of", nrow(loo_summary), "\n")
cat("Written to", f_sum, "and", f_det, "\n")
