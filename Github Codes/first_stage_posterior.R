############################################################
## first_stage_posterior.R
##
## Posterior summaries of the first-stage instrument effects
## (Pi_X for the ligand, Pi_Z for the receptor), triplet by
## triplet.
##
## WHY. Median first-stage F statistics are low in these data
## (about 1-2), so it is natural to ask whether the instruments
## carry signal inside the model itself. The sampler draws
## Pi_X and Pi_Z at every iteration but by default keeps only
## their means. With .mrccc_first_stage <- TRUE,
## real_data_analysis.R asks it to keep the draws and calls
## the two functions below.
##
## USAGE (from the project root, in a fresh R session):
##   .mrccc_cell1 <- "NKCells"; .mrccc_cell2 <- "MonocytesCells"
##   .mrccc_first_stage <- TRUE
##   Rcpp::sourceCpp("Github Codes/mr_ccc_gibbs.cpp")
##   source("Github Codes/build_lr_database.R")
##   source("Github Codes/real_data_analysis.R")
##
## REPRODUCIBILITY. Keeping the draws and summarising them
## uses no random numbers, so such a run reproduces the
## primary run exactly (same seed, same data, same package
## versions). Its outputs are written with the prefix
## "firststage_" and its PIPs are compared with those of the
## primary output at the end.
##
## SUMMARIES, for each gene of each triplet (suffix _L for
## the ligand instruments G, _R for the receptor instruments H):
##   k        number of instruments
##   n_excl0  number of instruments whose 95% equal-tailed
##            credible interval excludes 0
##   d2       squared Mahalanobis distance of 0 from the
##            posterior mean of the instrument-effect vector,
##            in the metric of its posterior covariance
##   p_joint  upper chi-square(k) tail probability of d2. Under
##            the approximately normal posterior, 0 lies outside
##            the joint 95% credible ellipsoid exactly when
##            p_joint < 0.05; this is the joint counterpart of
##            n_excl0, which looks at one instrument at a time
##   r2_med, r2_lo, r2_hi
##            posterior median and 95% interval of the share of
##            expression variance explained by the instruments,
##            var(G Pi_X) / var(X) (and likewise for Z)
##
## OUTPUT (written to Results/):
##   first_stage_posterior_<Cell1>_<Cell2>.csv
## prefixed "smoke_" when the run was a check.
############################################################

## Number of draws used for the variance-share summary. A systematic
## subsample, as for the effect-curve band in real_data_analysis.R, so that
## no random numbers are used.
FS_DRAWS <- 2000L

first_stage_summary <- function(fits, field, D, y_c, side) {
  dr <- do.call(rbind, lapply(fits, function(f) f[[field]]))
  if (is.null(dr) || nrow(dr) < 2L) {
    stop("No first-stage draws found in '", field,
         "'; was the sampler called with save_first_stage = TRUE?")
  }
  k <- ncol(dr)

  ## One instrument at a time.
  q <- apply(dr, 2, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
  q <- matrix(q, nrow = 2)
  excl0 <- q[1, ] > 0 | q[2, ] < 0

  ## All instruments jointly.
  m  <- colMeans(dr)
  S  <- stats::cov(dr)
  d2 <- tryCatch(as.numeric(crossprod(m, solve(S, m))),
                 error = function(e) NA_real_)
  p_joint <- if (is.finite(d2)) {
    stats::pchisq(d2, df = k, lower.tail = FALSE)
  } else NA_real_

  ## Share of expression variance explained by the instruments.
  idx <- if (nrow(dr) > FS_DRAWS) {
    round(seq(1, nrow(dr), length.out = FS_DRAWS))
  } else seq_len(nrow(dr))
  fitted <- D %*% t(dr[idx, , drop = FALSE])
  share  <- apply(fitted, 2, stats::var) / stats::var(as.numeric(y_c))
  qs <- stats::quantile(share, c(0.025, 0.5, 0.975), names = FALSE)

  out <- tibble::tibble(k = k, n_excl0 = sum(excl0), d2 = d2,
                        p_joint = p_joint,
                        r2_med = qs[2], r2_lo = qs[1], r2_hi = qs[3])
  names(out) <- paste0(names(out), "_", side)
  out
}

write_first_stage_posterior <- function(rows, output, Cell1, Cell2, prefix) {
  tab <- dplyr::bind_rows(rows)
  base_prefix <- sub("firststage_$", "", prefix)
  path <- file.path("Results",
                    paste0(base_prefix, "first_stage_posterior_",
                           Cell1, "_", Cell2, ".csv"))
  utils::write.csv(tab, path, row.names = FALSE)
  cat("First-stage posterior summaries written:", path, "\n")

  ## Overview for the triplets and for the declared discoveries.
  disc <- tab$PIP > 0.5
  overview <- function(t, lab) {
    cat(sprintf(paste0("  %-12s n = %2d | ligand: >=1 instrument excludes 0",
                       " in %2d, joint p < 0.05 in %2d | receptor: %2d, %2d\n"),
                lab, nrow(t),
                sum(t$n_excl0_L >= 1), sum(t$p_joint_L < 0.05, na.rm = TRUE),
                sum(t$n_excl0_R >= 1), sum(t$p_joint_R < 0.05, na.rm = TRUE)))
  }
  cat("\n=========== FIRST-STAGE POSTERIOR ===========\n")
  overview(tab, "all triplets")
  overview(tab[disc, , drop = FALSE], "discoveries")

  ## Reproduction check against the primary output. The comparison is made
  ## only for a canonical run, since a check run is not expected to agree.
  primary <- file.path("Results", paste0(Cell1, "_", Cell2, "_MR_CCC.rds"))
  if (nzchar(base_prefix)) {
    cat("\nCheck run: comparison with the primary output skipped.\n")
  } else if (!file.exists(primary)) {
    cat("\nPrimary output", primary, "not found: comparison skipped.\n")
  } else {
    key <- c("ligand_col_name", "receptor_col_name", "pathway_name")
    p0  <- readRDS(primary)[, c(key, "gamma_mean", "Beta_XZ_mean")]
    p1  <- output[, c(key, "gamma_mean", "Beta_XZ_mean")]
    j   <- dplyr::inner_join(p0, p1, by = key, suffix = c("_primary", "_rerun"))
    d_pip  <- max(abs(j$gamma_mean_primary - j$gamma_mean_rerun))
    d_beta <- max(abs(j$Beta_XZ_mean_primary - j$Beta_XZ_mean_rerun))
    cat(sprintf(paste0("\nReproduction check (%d of %d triplets matched): ",
                       "max |PIP difference| = %.3g, ",
                       "max |beta_XZ difference| = %.3g\n"),
                nrow(j), nrow(p0), d_pip, d_beta))
    if (nrow(j) == nrow(p0) && d_pip == 0 && d_beta == 0) {
      cat("REPRODUCED EXACTLY: the rerun matches the primary analysis.\n")
    } else {
      warning("The rerun does not reproduce the primary analysis exactly; ",
              "check package versions against the primary sessionInfo.",
              call. = FALSE)
    }
  }
  invisible(tab)
}
