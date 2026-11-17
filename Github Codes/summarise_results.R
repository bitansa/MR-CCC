############################################################
## MR-CCC: summaries that combine several saved result files
##
## Part 1 reads the twenty Results/<Cell1>_<Cell2>_MR_CCC.rds files written
## by real_data_analysis.R and reports the screen-wide numbers of the paper:
## the numbers of triplets and discoveries, the discoveries at pair-level
## Bayesian FDR <= 0.10, the largest Gelman-Rubin statistic, the triplets
## whose PIP lies within two Monte Carlo standard errors of 0.5, and the
## discoveries that stay above 0.5 under a neutral prior inclusion
## probability.
##
## Part 2 reads Results/simulation_S4_weak.csv written by
## simulation_misspecification.R and reports, for each sample size and
## instrument strength, the area under the ROC curve (AUC) of the three MR
## methods, and the type I error rate and power of MR-CCC under the neutral
## prior.
##
## Run from the repository root, after the analyses:
##   source("Github Codes/summarise_results.R")
## Nothing is resampled, so the output is exact.
############################################################

## Neutral prior. The prior inclusion probability rho ~ Beta(3, 1) has mean
## 0.75. The Bayes factor for gamma = 1 does not depend on this prior, so a
## prior with mean 0.5 divides the posterior odds by three:
##   PIP' = PIP / (PIP + 3 (1 - PIP)),
## and PIP' > 0.5 is the same as PIP > 0.75.
neutral_pip <- function(pip) pip / (pip + 3 * (1 - pip))

## Largest PIP-ranked set whose Bayesian FDR does not exceed alpha, the rule
## of real_data_analysis.R.
bfdr_set <- function(pip, alpha = 0.10) {
  o <- order(pip, decreasing = TRUE)
  f <- cumsum(1 - pip[o]) / seq_along(pip)
  k <- suppressWarnings(max(which(f <= alpha)))
  keep <- rep(FALSE, length(pip))
  if (is.finite(k) && k >= 1) keep[o[seq_len(k)]] <- TRUE
  keep
}

## Area under the ROC curve: the probability that a signal replicate scores
## above a null replicate, with ties counted as one half.
auc <- function(signal, null) {
  signal <- signal[is.finite(signal)]
  null <- null[is.finite(null)]
  if (length(signal) == 0L || length(null) == 0L) return(NA_real_)
  wins <- outer(signal, null, ">") + 0.5 * outer(signal, null, "==")
  mean(wins)
}

############################################################
## Part 1: real data, all twenty ordered cell-type pairs
############################################################
files <- list.files("Results",
                    pattern = "^[A-Za-z0-9]+_[A-Za-z0-9]+_MR_CCC\\.rds$",
                    full.names = TRUE)
if (length(files) != 20L) {
  stop("Expected the 20 pair files in Results/, found ", length(files),
       ". Run real_data_analysis.R for every pair first.", call. = FALSE)
}

need <- c("ligand_col_name", "receptor_col_name", "pathway_name",
          "gamma_mean", "PIP_mcse", "rhat_max")
screen <- do.call(rbind, lapply(files, function(f) {
  d <- as.data.frame(readRDS(f))
  miss <- setdiff(need, names(d))
  if (length(miss)) {
    stop(basename(f), " lacks the column(s) ", toString(miss), ".",
         call. = FALSE)
  }
  data.frame(pair      = sub("_MR_CCC\\.rds$", "", basename(f)),
             ligand    = d$ligand_col_name,
             receptor  = d$receptor_col_name,
             pathway   = d$pathway_name,
             pip       = d$gamma_mean,
             mcse      = d$PIP_mcse,
             rhat      = d$rhat_max,
             disc_bfdr = bfdr_set(d$gamma_mean, 0.10),
             stringsAsFactors = FALSE)
}))
screen$discovery <- screen$pip > 0.5
screen$pip_neutral <- neutral_pip(screen$pip)
screen$borderline <- abs(screen$pip - 0.5) < 2 * screen$mcse

cat("\n==================== Real data: all 20 pairs ====================\n")
cat("Triplets:", nrow(screen), "\n")
cat("Discoveries (PIP > 0.5):", sum(screen$discovery), "\n")
cat("Discoveries at pair-level Bayesian FDR <= 0.10:",
    sum(screen$disc_bfdr), "\n")
cat("Pairs with at least one discovery:",
    length(unique(screen$pair[screen$discovery])), "of",
    length(unique(screen$pair)), "\n")
cat("Largest R-hat:", round(max(screen$rhat, na.rm = TRUE), 3),
    "\nTriplets with R-hat above 1.01:",
    sum(screen$rhat > 1.01, na.rm = TRUE), "\n")
bad <- which(screen$rhat > 1.01)
if (length(bad)) {
  cat("Their distance from 0.5, in Monte Carlo standard errors:",
      round(abs(screen$pip[bad] - 0.5) / screen$mcse[bad], 1), "\n")
}
cat("Triplets whose PIP is within 2 MCSE of 0.5:",
    sum(screen$borderline, na.rm = TRUE),
    "\nDiscoveries whose PIP is within 2 MCSE of 0.5 (borderline):",
    sum(screen$borderline & screen$discovery, na.rm = TRUE),
    "\nDiscoveries that stay above 0.5 under the neutral prior:",
    sum(screen$discovery & screen$pip_neutral > 0.5), "of",
    sum(screen$discovery), "\n\n")

disc <- screen[screen$discovery, c("pair", "ligand", "receptor", "pathway",
                                   "pip", "pip_neutral", "mcse", "rhat",
                                   "disc_bfdr", "borderline")]
disc <- disc[order(disc$pair, -disc$pip), ]
rownames(disc) <- NULL
print(disc, digits = 3)

############################################################
## Part 2: scenario S4, weak instruments
############################################################
s4 <- utils::read.csv("Results/simulation_S4_weak.csv")
cells <- unique(s4[, c("n", "target_F")])
cells <- cells[order(cells$n, cells$target_F), ]

s4_summary <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  d <- s4[s4$n == cells$n[i] & s4$target_F == cells$target_F[i], ]
  sig <- d[d$arm == "signal", ]
  nul <- d[d$arm == "null", ]
  data.frame(n               = cells$n[i],
             target_F        = cells$target_F[i],
             median_F        = stats::median(nul$F_observed, na.rm = TRUE),
             AUC_MRCCC       = auc(sig$MR_score, nul$MR_score),
             AUC_MVMR        = auc(sig$MVMR_score, nul$MVMR_score),
             AUC_MRBMA       = auc(sig$MRBMA_score, nul$MRBMA_score),
             typeI_neutral   = mean(neutral_pip(nul$MR_score) > 0.5),
             power_neutral   = mean(neutral_pip(sig$MR_score) > 0.5))
}))

cat("\n============== Scenario S4: weak instruments ==============\n")
cat("AUC compares the signal and null arms. MR_score is the PIP of\n",
    "MR-CCC. The last two columns use the neutral prior.\n\n", sep = "")
print(s4_summary, digits = 2, row.names = FALSE)
