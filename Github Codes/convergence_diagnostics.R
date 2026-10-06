############################################################
## convergence_diagnostics.R
##
## MCMC TRACE PLOTS AND RUNTIME BENCHMARK
##
## The PRIMARY analysis (real_data_analysis.R) runs N_CHAINS chains per triplet
## from overdispersed starting values, reports the pooled posterior, and
## reports the Gelman-Rubin statistic for every triplet. This script refits
## selected triplets with the same settings, recomputes R-hat (via
## gelman_rubin()), ESS and MCSE for each monitored parameter from those
## refits, and adds the visual companion to that table, plus the timing
## figures:
##
##   1. TRACE PLOTS for representative triplets, showing all chains on the
##      same axes. Two series are drawn by default: the joint log-likelihood,
##      which is the conventional single-number summary of model fit, and the
##      running mean of the inclusion indicator, which is the posterior
##      inclusion probability as it accumulates. See TRACE_PARS below for why
##      both are shown rather than the log-likelihood alone.
##
##   2. A DIAGNOSTIC TABLE for the same triplets covering every monitored
##      quantity -- the two causal coefficients, the inclusion indicator and
##      the log-likelihood -- with R-hat, pooled effective sample size and
##      pooled Monte Carlo standard error, computed from the refitted
##      chains with the same estimators as the primary analysis.
##
##   3. RUNTIME as a function of sample size, together with a statement of
##      the per-sweep computational complexity.
##
## Why the inclusion indicator gets special attention. The posterior
## inclusion probability is a chain average of a binary indicator, and that
## chain is strongly autocorrelated, so its effective sample size is a
## small fraction of the number of retained draws. Since discovery
## is decided by comparing the PIP with a fixed threshold, the Monte Carlo
## error of that average is the quantity that determines whether a reported
## discovery is reproducible. R-hat and ESS are therefore reported for gamma
## as well as for the coefficients.
##
## PREREQUISITES (in this order):
##   Rcpp::sourceCpp("Github Codes/mr_ccc_gibbs.cpp")
##   source("Github Codes/build_lr_database.R")
##   source("Github Codes/real_data_analysis.R")     # supplies the preprocessed objects
##   source("Github Codes/convergence_diagnostics.R")
##
## OUTPUT:
##   Plots/convergence_trace_<Cell1>_<Cell2>.pdf
##   Results/convergence_trace_<Cell1>_<Cell2>.rds   (display series behind
##                                                    the trace plot, so the
##                                                    figure can be rebuilt
##                                                    by regenerate_figures.R)
##   Results/convergence_diagnostics_<Cell1>_<Cell2>.csv
##   Results/runtime_benchmark.csv
##   Results/sessionInfo_convergence_diagnostics.txt
## Each name is prefixed "smoke_" when the primary run was a check.
############################################################
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

## Theme, chain palette, axis helpers, save_fig() and plot_trace(), shared
## with every other figure in the repository.
source("Github Codes/figure_style.R")

.needed <- c("RNA.count.adj_Cell1", "RNA.count.adj_Cell2", "donor", "V_donor", "centre_cols", "orient_pc1",
             "get_SNP_matrix", "lr_filtered", "key_cols", "Output_MR_CCC",
             "Cell1", "Cell2", "N_ITER", "BURN_IN", "THIN", "pip_thresh",
             "gsva_scores", "mcse_batch", "mcse_pooled",
             "N_CHAINS", "INIT_SCALE")
.missing <- .needed[!vapply(.needed, exists, logical(1))]
if (length(.missing) > 0) {
  stop("Objects missing -- source real_data_analysis.R first.\n  Missing: ",
       paste(.missing, collapse = ", "))
}

## ---- Output names ---------------------------------------------------------
## A session whose primary run was a check (shortened chains, a truncated
## triplet universe or a non-canonical clumping threshold) writes to names
## prefixed "smoke_", by the rule used in real_data_analysis.R, so that a
## check never overwrites a reported result. The prefix set by the primary
## script is reused; the rule is re-evaluated on the inherited settings.
.out_prefix <- if ((exists(".out_prefix") && identical(.out_prefix, "smoke_")) ||
                   isTRUE(get0(".mrccc_truncated", ifnotfound = FALSE)) ||
                   N_ITER != 100000L ||
                   !isTRUE(all.equal(get0("LD_R2_MAX", ifnotfound = 0.8), 0.8))) {
  "smoke_"
} else ""

############################################################
## USER SETTINGS
############################################################
SEED_DIAG <- 20260921
set.seed(SEED_DIAG)

## N_CHAINS, INIT_SCALE, N_ITER, BURN_IN and THIN are all INHERITED from the
## primary script, so the chains drawn here use exactly the configuration the
## paper reports (each chain is seeded separately below, so the draws are a
## fresh replicate of that configuration). Nothing about the sampler
## configuration is overridden. The
## assertion below states that requirement rather than leaving it implicit: a
## trace figure produced under a different configuration from the reported
## results would misrepresent them.
stopifnot(
  "convergence_diagnostics.R requires N_CHAINS >= 2 (inherited from the primary)" =
    N_CHAINS >= 2,
  "convergence_diagnostics.R requires dispersed starts: INIT_SCALE > 0" =
    INIT_SCALE > 0
)

## Which series are drawn in the trace figure.
##
##   "loglik"   the joint log-likelihood of the three equations at each
##              retained draw. This is the conventional quantity to show when
##              a model has too many parameters to trace one by one: it is a
##              single scalar that every parameter feeds into, so drift
##              anywhere in the model tends to show up in it.
##
##   "gamma"    the running mean of the inclusion indicator, which IS the
##              posterior inclusion probability as it accumulates.
##
## Both are drawn by default, and the reason for not showing the
## log-likelihood alone is specific to this model. The reported quantity is
## the PIP, a chain average of a binary indicator, and the log-likelihood is
## dominated by the three residual variances and the first-stage blocks. It is
## therefore entirely possible for the log-likelihood to appear well mixed
## while the indicator is still moving between the spike and the slab. Showing
## the indicator alongside it means the figure speaks directly to the
## quantity that decides a discovery.
##
## "Beta_X" and "Beta_XZ" may be added; they are always present in the
## diagnostic table regardless of what is plotted.
TRACE_PARS <- c("loglik", "gamma")

## Quantities that appear in the diagnostic table, independently of TRACE_PARS.
DIAG_PARS  <- c("Beta_X", "Beta_XZ", "gamma", "loglik")

stopifnot(
  "TRACE_PARS must be a subset of DIAG_PARS" = all(TRACE_PARS %in% DIAG_PARS)
)

## Triplets to diagnose.
##
##   SPANNING SET: the strongest discovery, the triplet closest to the
##   decision threshold, and the weakest non-discovery. These show sampler
##   behaviour across the range of the posterior rather than only where it
##   is best behaved.
##
##   FLAGGED SET: every triplet whose R-hat exceeds 1.01 on any monitored
##   parameter in the primary output. These are included DELIBERATELY.
##   Presenting trace plots only for well-behaved triplets while others carry
##   flags would be selective; a reader assessing convergence should see
##   precisely the cases the diagnostic singled out.
##
## TRACE_SET controls which triplets are drawn:
##
##   "discoveries"  every triplet above the discovery threshold. This is the
##                  set a reader cares about, and it fits a journal figure.
##
##   "diagnostic"   the spanning set plus every R-hat-flagged triplet. Wider,
##                  and useful while assessing the sampler, because it shows a
##                  borderline case and a clear non-discovery for contrast and
##                  does not hide the cases the diagnostic singled out.
##
## Note that these differ: the flagged set and the discovery set may overlap
## only partially, so "discoveries" can omit a flagged triplet and vice versa.
## Whichever is used for the figure, R-hat is reported for ALL triplets by the
## primary script, so nothing is concealed by the choice made here.
TRACE_SET <- "discoveries"

RHAT_FLAG <- 1.01

pick_triplets <- function(out, n_each = 1L) {
  o <- out %>% arrange(desc(gamma_mean))

  # Spanning set: strongest, closest to the threshold, weakest.
  idx_span <- unique(c(
    head(seq_len(nrow(o)), n_each),
    order(abs(o$gamma_mean - pip_thresh))[seq_len(n_each)],
    tail(seq_len(nrow(o)), n_each)
  ))

  idx_disc <- if ("disc_pip" %in% names(o)) which(o$disc_pip) else integer(0)

  idx_flag <- if ("rhat_max" %in% names(o)) {
    which(is.finite(o$rhat_max) & o$rhat_max > RHAT_FLAG)
  } else integer(0)

  idx <- switch(
    TRACE_SET,
    discoveries = if (length(idx_disc)) idx_disc else idx_span,
    diagnostic  = unique(c(idx_span, idx_flag)),
    stop("TRACE_SET must be \"discoveries\" or \"diagnostic\".")
  )

  keep <- intersect(c("ligand_col_name", "receptor_col_name", "pathway_name",
                      "gamma_mean", "disc_pip", "rhat_max"), names(o))
  out_sel <- o[idx, keep, drop = FALSE]
  out_sel$rhat_flagged <- idx %in% idx_flag
  out_sel
}

############################################################
## HELPERS
############################################################

## Gelman-Rubin R-hat for a list of equal-length chains (Gelman & Rubin,
## 1992). Values near 1 indicate that between-chain and within-chain
## variance agree; a common working rule is R-hat < 1.01.
gelman_rubin <- function(chains) {
  m <- length(chains); n <- length(chains[[1]])
  if (m < 2 || n < 2) return(NA_real_)
  means <- vapply(chains, mean, numeric(1))
  vars  <- vapply(chains, stats::var, numeric(1))
  W <- mean(vars)
  B <- n * stats::var(means)
  if (!is.finite(W)) return(NA_real_)
  # W == 0 means every chain is constant. If the chains are constant at
  # DIFFERENT values (B > 0) that is total non-convergence -- two chains
  # stuck at gamma = 1 and two at gamma = 0 -- and must be reported as
  # Inf, not NA, or na.rm = TRUE downstream would count it as converged.
  # Constant AND identical (B == 0) is genuinely undefined.
  if (W <= 0) return(if (is.finite(B) && B > 0) Inf else NA_real_)
  var_hat <- ((n - 1) / n) * W + B / n
  sqrt(var_hat / W)
}

## Assemble one triplet exactly as the primary pipeline does, so the chains
## being diagnosed are the chains actually used for inference.
##
## The assembly logic lives in build_triplet_inputs.R and is SHARED with the
## other post-hoc refit scripts. A single shared definition ensures that all
## six pathway representations are offered to the adaptive rule, as in the
## primary analysis, and the representation chosen is checked against the
## one recorded in Output_MR_CCC, with a warning if they differ. One
## definition, one place.
source("Github Codes/build_triplet_inputs.R")

build_inputs <- function(lig_sym, rec_sym, pth) {
  build_triplet_inputs(lig_sym, rec_sym, pth)
}

## Run N_CHAINS chains from dispersed starts.
run_chains <- function(inp) {
  n <- nrow(inp$X); g_scale <- min(n, 100.0)
  lapply(seq_len(N_CHAINS), function(cc) {
    set.seed(SEED_DIAG + 1000L * cc)
    mr_ccc_gibbs(
      inp$X, inp$Z, inp$Y, inp$G, inp$H, inp$V,
      n_iter = N_ITER, burn_in = BURN_IN, thin = THIN,
      a_sigma = 3.0, b_sigma = 2.0, a_rho = 3.0, b_rho = 1.0,
      nu1 = 1e-4, gG = g_scale, gH = g_scale, gV = g_scale,
      gZ = g_scale, gBeta = g_scale, ridge = 1e-8,
      init_gamma = if (cc %% 2 == 0) 0L else 1L,   # alternate spike/slab start
      init_scale = INIT_SCALE                       # dispersed beta starts
    )
  })
}

############################################################
## PART 1: TRACE PLOTS AND R-hat
############################################################
sel <- pick_triplets(Output_MR_CCC)
cat("\nDiagnosing", nrow(sel), "triplets with", N_CHAINS, "chains each.\n")
print(as.data.frame(sel))

diag_rows <- list(); trace_rows <- list(); d <- 0
for (i in seq_len(nrow(sel))) {

  inp <- build_inputs(sel$ligand_col_name[i], sel$receptor_col_name[i],
                      sel$pathway_name[i])
  if (is.null(inp)) next
  fits <- run_chains(inp)
  lab  <- paste0(sel$ligand_col_name[i], "–", sel$receptor_col_name[i])
  # Facet column label for the trace plot: ligand on the first line and
  # receptor on the second, so that long gene symbols do not widen the
  # panels. The table label `lab` keeps the single-line form.
  lab_facet <- paste0(sel$ligand_col_name[i], "\n", sel$receptor_col_name[i])

  for (par in DIAG_PARS) {
    ch <- lapply(fits, function(f) as.numeric(f[[paste0(par, "_draws")]]))
    # Per-chain batch means, combined in quadrature -- the SAME estimator as
    # real_data_analysis.R. Applying mcse_batch() to the concatenated chains
    # would place batches across the joins between chains that started at
    # different values and inflate the variance, so the ESS/MCSE written here
    # would not certify the PIP_ess/PIP_mcse of the primary analysis.
    mc <- mcse_pooled(lapply(ch, mcse_batch))
    d <- d + 1
    diag_rows[[d]] <- tibble(
      triplet = lab, parameter = par,
      rhat = gelman_rubin(ch),
      ess  = unname(mc["ess"]),
      mcse = unname(mc["mcse"]),
      post_mean = mean(unlist(ch))
    )

    # ---- Display series ------------------------------------------------
    # Only the series named in TRACE_PARS are carried into the figure; the
    # diagnostics above cover every parameter regardless.
    if (!(par %in% TRACE_PARS)) next

    # The coefficients and the log-likelihood are plotted as ordinary traces.
    # The inclusion indicator is NOT: it takes only the values 0 and 1, so a
    # raw trace of several hundred thousand draws is a solid block of vertical
    # lines that conveys nothing. It is plotted instead as a RUNNING MEAN,
    #
    #   running_mean(t) = (1/t) * sum_{s<=t} gamma_s,
    #
    # which is the posterior inclusion probability as it accumulates. A
    # converged chain's running mean settles to a horizontal asymptote, and
    # several chains settling to the SAME level is the visual counterpart of
    # the Gelman-Rubin statistic. The running mean is computed on the FULL
    # chain and thinned only afterwards, so thinning affects which points are
    # drawn but not the value of the curve.
    #
    # Diagnostics above are always computed on the full, unthinned chains.
    par_label <- switch(par,
      gamma  = "gamma (running mean)",
      loglik = "log-likelihood",
      par)
    for (cc in seq_along(ch)) {
      series <- if (par == "gamma") {
        cumsum(ch[[cc]]) / seq_along(ch[[cc]])
      } else ch[[cc]]
      keep <- seq(1, length(series), by = max(1, floor(length(series) / 2000)))
      trace_rows[[length(trace_rows) + 1L]] <- tibble(
        triplet = lab_facet, parameter = par_label, chain = factor(cc),
        iter = keep, value = series[keep]
      )
    }
  }
  cat("  done:", lab, "\n")
}
diag_tab  <- bind_rows(diag_rows)
trace_tab <- bind_rows(trace_rows)

# The run settings travel with the display series as attributes, so that
# plot_trace() can state them in the subtitle when the figure is rebuilt
# from the saved object.
attr(trace_tab, "n_chains") <- N_CHAINS
attr(trace_tab, "n_iter")   <- N_ITER
attr(trace_tab, "burn_in")  <- BURN_IN

dir.create("Results", showWarnings = FALSE)
write.csv(diag_tab,
          paste0("Results/", .out_prefix, "convergence_diagnostics_",
                 Cell1, "_", Cell2, ".csv"),
          row.names = FALSE)
saveRDS(trace_tab,
        paste0("Results/", .out_prefix, "convergence_trace_",
               Cell1, "_", Cell2, ".rds"))

cat("\n---- Convergence diagnostics ----\n")
print(as.data.frame(diag_tab %>% mutate(across(where(is.numeric), ~round(.x, 4)))))

############################################################
## TRACE PLOT
##
## Built by plot_trace() in figure_functions.R (sourced through
## figure_style.R), which is the single definition of this figure;
## regenerate_figures.R calls the same function on the saved trace_tab.
## Theme, chain colours (Okabe-Ito) and the iteration axis (at most four
## breaks, in thousands) come from figure_style.R.
############################################################
p_trace <- plot_trace(trace_tab)
save_fig(p_trace, paste0(.out_prefix, "convergence_trace_", Cell1, "_", Cell2),
         14, 12)

############################################################
## PART 2: RUNTIME BENCHMARK AND COMPLEXITY
##
## Per sweep the cost is dominated by the two first-stage blocks and the
## outcome block, each of which forms cross-products of an n x p design and
## inverts a p x p matrix. With p_G, p_H and p_V small and fixed, the cost is
## therefore LINEAR in the number of donors n and linear in the number of
## iterations; triplets are independent and embarrassingly parallel.
############################################################
bench_n <- c(250, 500, 1000, 2000)
bench <- lapply(bench_n, function(nn) {
  set.seed(SEED_DIAG)
  G <- matrix(rnorm(nn * 5), nn); H <- matrix(rnorm(nn * 5), nn)
  V <- matrix(rnorm(nn * 3), nn)
  X <- matrix(G %*% rep(.5, 5) + rnorm(nn), nn)
  Z <- matrix(H %*% rep(.5, 5) + rnorm(nn), nn)
  Y <- matrix(.3 * X + .3 * X * Z + rnorm(nn), nn)
  tt <- system.time(
    mr_ccc_gibbs(X, Z, Y, G, H, V, n_iter = 2000, burn_in = 200, thin = 1)
  )[["elapsed"]]
  tibble(n = nn, seconds_per_2000_iter = tt,
         seconds_per_100k_iter = tt * 50)
})
bench_tab <- bind_rows(bench)
write.csv(bench_tab, paste0("Results/", .out_prefix, "runtime_benchmark.csv"),
          row.names = FALSE)

cat("\n---- Runtime benchmark (single triplet) ----\n")
print(as.data.frame(bench_tab %>% mutate(across(where(is.numeric), ~round(.x, 2)))))
cat("\nCost is linear in n and in the number of iterations; triplets are\n",
    "independent, so genome-scale screens parallelise across cores.\n", sep = "")

## Package versions used for this run.
writeLines(utils::capture.output(sessionInfo()),
           paste0("Results/", .out_prefix,
                  "sessionInfo_convergence_diagnostics.txt"))
