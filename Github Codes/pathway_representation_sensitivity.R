############################################################
## pathway_representation_sensitivity.R
##
## SENSITIVITY OF THE RESULTS TO THE CHOICE OF PATHWAY REPRESENTATION
##
## The primary analysis computes six donor-level pathway activity scores for
## every triplet -- first principal component, mean expression, AUCell,
## UCell, ssGSEA and GSVA -- and retains the one whose absolute Pearson
## correlation with the OBSERVED ligand expression X is largest.
##
## That rule selects the outcome using its association with the exposure,
## which raises two distinct questions:
##
##   (1) ROBUSTNESS.  Would the same triplets be discovered if a single
##       representation were fixed in advance instead of being chosen?
##
##   (2) CALIBRATION. Does the selection step itself inflate posterior
##       inclusion probabilities, by allowing the most favourable of six
##       candidate outcomes to be picked?
##
## Question (1) is answered by re-running the analysis six times with each
## representation held fixed, and comparing both the discovery sets and the
## full PIP rankings.
##
## Question (2) is answered by a permutation experiment. The six Y vectors
## are computed as usual and then their DONOR ORDER is permuted -- the same
## permutation applied to all six -- which preserves each representation's
## marginal distribution while destroying any association with X. Under this
## null, a well-calibrated procedure should return low PIPs. Running the
## adaptive rule and a fixed representation through identical permutations
## isolates the inflation attributable to selection alone.
##
## A mitigating argument, which the permutation experiment tests rather than
## assumes: selection uses the OBSERVED X, whereas inference uses the
## instrument-projected X*, and the confounded component of X that drives
## much of the observed correlation is precisely what the MR projection
## removes. The two are nevertheless correlated, so the question is empirical.
##
## PREREQUISITES (in this order):
##   Rcpp::sourceCpp("Github Codes/mr_ccc_gibbs.cpp")
##   source("Github Codes/build_lr_database.R")
##   source("Github Codes/real_data_analysis.R")      # supplies the preprocessed objects
##   source("Github Codes/pathway_representation_sensitivity.R")
##
## Re-using the objects built by real_data_analysis.R guarantees that the QC,
## normalisation, instrument selection and MCMC settings are identical to the
## primary analysis, so the only thing that varies here is the outcome
## representation.
##
## OUTPUT (written to Results/):
##   representation_fixed_<Cell1>_<Cell2>.csv   PIPs under each fixed choice
##   representation_concord_<Cell1>_<Cell2>.csv rank concordance summary
##   representation_null_<Cell1>_<Cell2>.csv    permutation null calibration
##   sessionInfo_pathway_representation_sensitivity.txt
## Each name is prefixed "smoke_" when the primary run was a check. The fixed
## and null tables carry PIP_rhat, the Gelman-Rubin statistic of the
## inclusion indicator across the inherited chains (NA for a single chain).
############################################################
library(dplyr)
library(tidyr)
library(tibble)
library(AUCell)
library(UCell)
library(GSVA)

## ---- Objects inherited from real_data_analysis.R -------------------------
.needed <- c("RNA.count.adj_Cell1", "RNA.count.adj_Cell2", "donor", "V_donor", "centre_cols", "orient_pc1",
             "genotype.mat", "vcf.ref", "gene.ref", "gene.promoter.ref",
             "lr_filtered", "key_cols", "get_SNP_matrix", "gsva_scores",
             "Cell1", "Cell2", "N_ITER", "BURN_IN", "THIN",
             # The chain settings and R-hat function are inherited too, so
             # that the "adaptive" arm here reproduces the primary analysis
             # exactly rather than approximating it with a single chain.
             "N_CHAINS", "INIT_SCALE", "gelman_rubin",
             "pip_thresh", "bayes_fdr",
             # Instrument selection is LD-clumped; both objects are defined
             # by real_data_analysis.R alongside get_SNP_matrix().
             "LD_R2_MAX", "clump_by_ld")
.missing <- .needed[!vapply(.needed, exists, logical(1))]
if (length(.missing) > 0) {
  stop("Objects missing -- source real_data_analysis.R first.\n  Missing: ",
       paste(.missing, collapse = ", "))
}

## check_primary_representation(), shared with the other refit scripts,
## compares the representation chosen here with the one recorded in
## Output_MR_CCC.
source("Github Codes/build_triplet_inputs.R")

## ---- Output names ---------------------------------------------------------
## Prefixed "smoke_" when the primary run in this session was a check, by the
## rule used in real_data_analysis.R, so that a check never overwrites a
## reported result.
.out_prefix <- if ((exists(".out_prefix") && identical(.out_prefix, "smoke_")) ||
                   isTRUE(get0(".mrccc_truncated", ifnotfound = FALSE)) ||
                   N_ITER != 100000L ||
                   !isTRUE(all.equal(LD_R2_MAX, 0.8))) {
  "smoke_"
} else ""

############################################################
## USER SETTINGS
############################################################
SEED_REP <- 20260921
set.seed(SEED_REP)

REPRESENTATIONS <- c("PC1", "AUCell", "Mean", "UCell", "ssGSEA", "GSVA")

## Representation to promote if the primary analysis is ever switched from
## the adaptive rule to a pre-specified one. PC1 is the natural default: it
## is defined without reference to X and requires no external gene ranking.
PRESPECIFIED <- "PC1"

## Number of permutations for the null-calibration experiment. Each
## permutation re-runs every triplet under two rules (adaptive and fixed),
## so cost is  2 x N_PERM x n_triplets  Gibbs runs. Start small.
N_PERM <- 10

## Optionally restrict the permutation experiment to a subset of triplets:
## the first PERM_MAX_TRIPLETS triplets in processing order (the order of
## lr_filtered, skipping triplets that cannot be built). Set to NULL to use
## all triplets.
PERM_MAX_TRIPLETS <- 20

############################################################
## HELPER: the six pathway representations
##
## Reproduces exactly the construction used in real_data_analysis.R, but
## returns ALL six rather than only the selected one.
############################################################
compute_representations <- function(pathway_gene_ids, RNA_receiver) {

  expr      <- RNA_receiver[pathway_gene_ids, , drop = FALSE]
  expr_full <- RNA_receiver

  # Oriented against the mean score exactly as in the primary pipeline.
  Y_PC1 <- orient_pc1(prcomp(t(expr), center = TRUE, scale. = TRUE)$x[, 1],
                      colMeans(expr))

  geneRanks <- AUCell_buildRankings(expr_full, plotStats = FALSE)
  path_genes <- intersect(rownames(expr_full), pathway_gene_ids)
  geneSets   <- list(pathway = path_genes)
  auc <- AUCell_calcAUC(geneSets, geneRanks,
                        aucMaxRank = ceiling(0.05 * nrow(expr_full)))
  Y_AUCell <- as.numeric(getAUC(auc)["pathway", ])

  Y_Mean  <- colMeans(expr)
  Y_UCell <- as.numeric(
    ScoreSignatures_UCell(expr_full, features = geneSets)[, "pathway_UCell"])
  # gsva_scores() is inherited from real_data_analysis.R and tolerates both
  # the legacy and the parameter-object GSVA interfaces.
  Y_ssGSEA <- as.numeric(
    gsva_scores(expr_full, geneSets, method = "ssgsea")["pathway", ])
  Y_GSVA <- as.numeric(
    gsva_scores(expr_full, geneSets, method = "gsva",
                kcdf = "Gaussian")["pathway", ])

  list(PC1 = Y_PC1, AUCell = Y_AUCell, Mean = Y_Mean,
       UCell = Y_UCell, ssGSEA = Y_ssGSEA, GSVA = Y_GSVA)
}

# The adaptive rule used by the primary analysis.
select_adaptive <- function(Y_list, X) {
  cors <- vapply(Y_list, function(y) cor(y, X, use = "complete.obs"), numeric(1))
  names(which.max(abs(cors)))
}

############################################################
## HELPER: fit MR-CCC for one triplet given a chosen Y
##
## Centring, standardisation, covariates, priors and MCMC settings are taken
## from the primary analysis, so the fitted model differs only in Y.
############################################################

## Counter for fits that could not be completed. It lives in an environment
## rather than a plain variable so that fit_one() can update it without
## resorting to <<- into the global environment.
.fit_failures <- new.env(parent = emptyenv())
.fit_failures$n   <- 0L
.fit_failures$msg <- character(0)

## Reports and resets the counter. Called at the end of each part, so that a
## representation that fails on some triplets is visible in the output rather
## than silently producing a smaller table.
report_fit_failures <- function(part) {
  if (.fit_failures$n > 0L) {
    cat("\n  NOTE:", .fit_failures$n, "fit(s) in", part,
        "could not be completed and were skipped.\n")
    cat("        Reason(s):", paste(.fit_failures$msg, collapse = "; "), "\n")
    cat("        This occurs when a FIXED representation yields an outcome",
        "that is\n        near-constant or numerically collinear with the",
        "projected exposures.\n        It is a property of that",
        "representation on that triplet.\n")
  } else {
    cat("\n  All fits in", part, "completed.\n")
  }
  .fit_failures$n   <- 0L
  .fit_failures$msg <- character(0)
  invisible(NULL)
}

fit_one <- function(X, Z, Y, G, H, V) {
  sd_X <- sd(as.numeric(X)); sd_Z <- sd(as.numeric(Z)); sd_Y <- sd(as.numeric(Y))
  if (!is.finite(sd_Y) || sd_Y == 0) return(NULL)

  # A FIXED representation can produce an outcome the adaptive rule would
  # never have selected: a pathway score that is near-constant, or so nearly
  # collinear with the projected exposures that the causal block's
  # cross-product matrix is numerically singular. The sampler then fails in
  # inv_sympd(). That is a property of the representation on that triplet,
  # not an error in the analysis, and it must not terminate a sweep over six
  # representations and every triplet. Such fits are skipped and counted; the
  # count is reported at the end of each part so the number of skips is
  # visible rather than silent.
  if (!all(is.finite(as.numeric(Y)))) return(NULL)
  # Near-constant outcome relative to its own scale: the spike-and-slab has
  # nothing to explain and the g-prior scale degenerates.
  if (sd_Y < 1e-8 * max(1, abs(mean(as.numeric(Y))))) return(NULL)

  X_c <- matrix(as.numeric(X) - mean(X), ncol = 1)
  Z_c <- matrix(as.numeric(Z) - mean(Z), ncol = 1)
  Y_c <- matrix(as.numeric(Y) - mean(Y), ncol = 1)

  # Sampler configuration inherited from real_data_analysis.R (N_CHAINS,
  # N_ITER, BURN_IN, THIN, INIT_SCALE) so it is IDENTICAL to the primary.
  #
  # This matters for the comparison, not just for tidiness. The "adaptive"
  # arm of this sensitivity analysis is meant to REPRODUCE the primary
  # result. Any difference in chain length or start would let a gap between
  # a fixed representation and the adaptive rule be Monte Carlo noise rather
  # than the effect of the representation, and the comparison would not
  # isolate the effect of the representation choice.
  n <- nrow(X_c); g_scale <- min(n, 100.0)
  fits <- tryCatch(
    lapply(seq_len(N_CHAINS), function(cc) {
      mr_ccc_gibbs(
        X_c, Z_c, Y_c, G, H, V,
        n_iter = N_ITER, burn_in = BURN_IN, thin = THIN,
        a_sigma = 3.0, b_sigma = 2.0, a_rho = 3.0, b_rho = 1.0,
        nu1 = 1e-4, gG = g_scale, gH = g_scale, gV = g_scale,
        gZ = g_scale, gBeta = g_scale, ridge = 1e-8,
        init_gamma = if (cc %% 2L == 0L) 0L else 1L,
        init_scale = INIT_SCALE
      )
    }),
    error = function(e) {
      .fit_failures$n   <- .fit_failures$n + 1L
      .fit_failures$msg <- unique(c(.fit_failures$msg, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(fits)) return(NULL)

  g_ch <- lapply(fits, function(f) as.numeric(f$gamma_draws))
  list(PIP        = mean(unlist(g_ch, use.names = FALSE)),
       # Convergence on the same basis as the primary: the Gelman-Rubin
       # statistic across the inherited chains. NA if N_CHAINS == 1.
       PIP_rhat   = gelman_rubin(g_ch),
       Beta_X_std = mean(unlist(lapply(fits, function(f)
                      as.numeric(f$Beta_X_draws)), use.names = FALSE)) *
                      sd_X / sd_Y,
       Beta_XZ_std= mean(unlist(lapply(fits, function(f)
                      as.numeric(f$Beta_XZ_draws)), use.names = FALSE)) *
                      sd_X * sd_Z / sd_Y)
}

############################################################
## HELPER: assemble one triplet's inputs
##
## Returns NULL for any triplet the primary pipeline would skip, so the
## triplet universe here matches the primary analysis exactly.
############################################################
prepare_triplet <- function(run) {

  pathway_name       <- lr_filtered$pathway_name[run]
  ligand_gene_id     <- lr_filtered$ligand_ensembl[run]
  receptor_gene_id   <- lr_filtered$receptor_ensembl[run]

  keyword <- sub(".*by\\s+", "", pathway_name)
  if (!nzchar(keyword)) return(NULL)
  hits <- key_cols %>%
    dplyr::filter(stringr::str_detect(gs_name,
                                      stringr::regex(keyword, ignore_case = TRUE)))
  pw <- unique(hits$ensembl_gene)
  pw <- pw[pw %in% rownames(RNA.count.adj_Cell2)]
  if (length(pw) == 0) return(NULL)

  if (!ligand_gene_id   %in% rownames(RNA.count.adj_Cell1)) return(NULL)
  if (!receptor_gene_id %in% rownames(RNA.count.adj_Cell2)) return(NULL)

  X <- matrix(RNA.count.adj_Cell1[ligand_gene_id,   ], ncol = 1)
  Z <- matrix(RNA.count.adj_Cell2[receptor_gene_id, ], ncol = 1)

  # Raw dosages are centred, and the covariates are the primary pipeline's
  # V_donor, centred -- identical to the inputs of the primary fit.
  G <- get_SNP_matrix(ligand_gene_id,   RNA.count.adj_Cell1, min_snps = 1)
  if (is.null(G)) return(NULL)
  H <- get_SNP_matrix(receptor_gene_id, RNA.count.adj_Cell2, min_snps = 1)
  if (is.null(H)) return(NULL)
  G <- centre_cols(G)
  H <- centre_cols(H)
  V <- centre_cols(V_donor)

  list(ligand   = lr_filtered$ligand_symbol[run],
       receptor = lr_filtered$receptor_symbol[run],
       pathway  = pathway_name,
       X = X, Z = Z, G = G, H = H, V = V,
       Y_list = compute_representations(pw, RNA.count.adj_Cell2))
}

############################################################
## PART 1: FIXED-REPRESENTATION RE-ANALYSIS
##
## Each triplet is refitted seven times: once under each of the six
## representations held fixed, and once under the adaptive rule (which
## reproduces the primary analysis and therefore acts as a control).
############################################################
cat("\n=== Part 1: fixed-representation re-analysis ===\n")

trip_cache <- list(); rows <- list(); k <- 0
for (run in seq_len(nrow(lr_filtered))) {
  tr <- prepare_triplet(run)
  if (is.null(tr)) next
  trip_cache[[length(trip_cache) + 1L]] <- tr

  chosen <- select_adaptive(tr$Y_list, tr$X)
  # The adaptive arm is meant to reproduce the primary fit; a different
  # representation is reported (the refit proceeds with the one chosen here).
  check_primary_representation(tr$ligand, tr$receptor, tr$pathway, chosen)

  for (rep_name in c(REPRESENTATIONS, "Adaptive")) {
    Y <- if (rep_name == "Adaptive") tr$Y_list[[chosen]] else tr$Y_list[[rep_name]]
    f <- fit_one(tr$X, tr$Z, matrix(Y, ncol = 1), tr$G, tr$H, tr$V)
    if (is.null(f)) next
    k <- k + 1
    rows[[k]] <- tibble(
      ligand = tr$ligand, receptor = tr$receptor, pathway = tr$pathway,
      representation = rep_name,
      adaptive_choice = chosen,
      PIP = f$PIP, Beta_X = f$Beta_X_std, Beta_XZ = f$Beta_XZ_std,
      PIP_rhat = f$PIP_rhat
    )
  }
  if (length(trip_cache) %% 10 == 0)
    cat("  ", length(trip_cache), "triplets done\n")
}
report_fit_failures("Part 1 (fixed representations)")

fixed_res <- bind_rows(rows)

dir.create("Results", showWarnings = FALSE)
write.csv(fixed_res,
          paste0("Results/", .out_prefix, "representation_fixed_",
                 Cell1, "_", Cell2, ".csv"),
          row.names = FALSE)

## ---- Discovery sets and rank concordance --------------------------------
## Spearman correlation of the full PIP ranking is reported alongside the
## overlap of discovery sets: two representations can agree on which
## triplets are discovered while ordering the remainder quite differently.
adaptive_pip <- fixed_res %>% filter(representation == "Adaptive") %>%
  select(ligand, receptor, pathway, PIP_adaptive = PIP)

## The join below is keyed on gene symbols and pathway. A key that occurs
## more than once would pair each fixed-representation row with several
## adaptive rows, so it is reported instead of joined.
.dup_keys <- adaptive_pip %>% count(ligand, receptor, pathway) %>% filter(n > 1L)
if (nrow(.dup_keys) > 0L) {
  stop("The (ligand, receptor, pathway) key is not unique for ",
       nrow(.dup_keys), " triplet(s), e.g. ",
       paste(.dup_keys$ligand[1], .dup_keys$receptor[1], .dup_keys$pathway[1],
             sep = " / "),
       "; the concordance join would duplicate rows.", call. = FALSE)
}

concord <- fixed_res %>%
  filter(representation != "Adaptive") %>%
  inner_join(adaptive_pip, by = c("ligand", "receptor", "pathway")) %>%
  group_by(representation) %>%
  summarise(
    n_triplets       = n(),
    n_disc_fixed     = sum(PIP > pip_thresh),
    n_disc_adaptive  = sum(PIP_adaptive > pip_thresh),
    n_disc_shared    = sum(PIP > pip_thresh & PIP_adaptive > pip_thresh),
    spearman_rank    = suppressWarnings(cor(PIP, PIP_adaptive, method = "spearman")),
    median_PIP_diff  = median(PIP - PIP_adaptive),
    .groups = "drop"
  )

write.csv(concord,
          paste0("Results/", .out_prefix, "representation_concord_",
                 Cell1, "_", Cell2, ".csv"),
          row.names = FALSE)

cat("\n---- Concordance with the adaptive rule ----\n")
print(as.data.frame(concord))

cat("\n---- Discoveries under each fixed representation ----\n")
print(as.data.frame(
  fixed_res %>% group_by(representation) %>%
    summarise(n_discoveries = sum(PIP > pip_thresh),
              median_PIP = round(median(PIP), 3), .groups = "drop")
))

############################################################
## PART 2: PERMUTATION NULL CALIBRATION
##
## The donor order of the six Y vectors is permuted (one permutation shared
## across all six), which preserves each representation's marginal
## distribution while removing any association with X. Every triplet is then
## refitted twice: once under the adaptive rule and once under the
## pre-specified representation.
##
## Interpretation. Under this null both procedures should return low PIPs.
## If the adaptive rule returns systematically higher PIPs -- or declares
## more discoveries -- than the fixed rule across permutations, that excess
## is the inflation attributable to selecting the outcome using X.
############################################################
cat("\n=== Part 2: permutation null calibration ===\n")

perm_pool <- if (is.null(PERM_MAX_TRIPLETS)) seq_along(trip_cache) else
  head(seq_along(trip_cache), PERM_MAX_TRIPLETS)

perm_rows <- list(); m <- 0
for (b in seq_len(N_PERM)) {
  idx <- sample.int(nrow(donor))          # one permutation, shared across reps
  for (ti in perm_pool) {
    tr <- trip_cache[[ti]]
    Yp <- lapply(tr$Y_list, function(y) y[idx])

    chosen_p <- select_adaptive(Yp, tr$X)
    for (rule in c("Adaptive", PRESPECIFIED)) {
      Y <- if (rule == "Adaptive") Yp[[chosen_p]] else Yp[[rule]]
      f <- fit_one(tr$X, tr$Z, matrix(Y, ncol = 1), tr$G, tr$H, tr$V)
      if (is.null(f)) next
      m <- m + 1
      perm_rows[[m]] <- tibble(
        perm = b, ligand = tr$ligand, receptor = tr$receptor,
        rule = rule, PIP = f$PIP, PIP_rhat = f$PIP_rhat
      )
    }
  }
  cat("  permutation", b, "of", N_PERM, "done\n")
}
report_fit_failures("Part 2 (permutation null)")

perm_res <- bind_rows(perm_rows)

write.csv(perm_res,
          paste0("Results/", .out_prefix, "representation_null_",
                 Cell1, "_", Cell2, ".csv"),
          row.names = FALSE)

cat("\n---- Null calibration (permuted outcomes) ----\n")
print(as.data.frame(
  perm_res %>% group_by(rule) %>%
    summarise(mean_PIP        = round(mean(PIP), 4),
              median_PIP      = round(median(PIP), 4),
              pct_above_thresh= round(100 * mean(PIP > pip_thresh), 2),
              .groups = "drop")
))
cat("\nA higher false positive rate for 'Adaptive' than for '", PRESPECIFIED,
    "' is the inflation attributable to selecting the outcome using X.\n", sep = "")

cat("\nWritten to Results/", .out_prefix, "representation_{fixed,concord,null}_",
    Cell1, "_", Cell2, ".csv\n", sep = "")

## Package versions used for this run.
writeLines(utils::capture.output(sessionInfo()),
           paste0("Results/", .out_prefix,
                  "sessionInfo_pathway_representation_sensitivity.txt"))
