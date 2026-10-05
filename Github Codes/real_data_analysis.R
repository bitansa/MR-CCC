############################################################
## real_data_analysis.R
##
## MR-CCC real data analysis: one ordered cell-type pair
## (Cell1 -> Cell2), iterating over all ligand-receptor-
## pathway triplets in lr_filtered.
##
## PREREQUISITES (run before sourcing this file):
##   Rcpp::sourceCpp("Github Codes/mr_ccc_gibbs.cpp")
##   source("Github Codes/build_lr_database.R")   # produces lr_matrix, key_cols
##
## DATA:
##   Download the two .rda files from Zenodo:
##     https://doi.org/10.5281/zenodo.19675075
##   Files:
##     B_T_NK_monocytes.rda  -- pseudo-bulk expression matrices
##                              for five cell types (genes x donors)
##     donor.rda             -- donor metadata, SNP genotype matrix,
##                              and GRanges objects for SNP/gene coords
##   Save both files to a local folder and set DATA_DIR below.
##
## OUTPUT:
##   Output_MR_CCC   -- one row per triplet, columns:
##                       ligand, receptor, pathway, posterior
##                       mean effect sizes (standardized),
##                       PIP (gamma_mean), Z_vec, Effect_vec,
##                       Effect_lo / Effect_hi (pointwise 95% band),
##                       95% credible intervals for Beta_X and Beta_XZ,
##                       sign-reversal threshold tau with interval,
##                       MCMC reliability (PIP_mcse, PIP_ess, n_chains,
##                       n_keep, and rhat_Beta_X / rhat_Beta_XZ /
##                       rhat_gamma / rhat_loglik / rhat_max, which are
##                       NA unless N_CHAINS >= 2),
##                       instrument strength (n_iv_L, n_iv_R, F_L, F_R),
##                       error control (bFDR, disc_pip, disc_bfdr)
##   Three ggplot objects: p1 (bubble), p2 (lollipop), p3 (curves),
##   built by the functions in figure_functions.R and styled by
##   figure_style.R, which every figure in this repository shares.
##   PDF files saved to Plots/ subdirectory.
##   RDS file saved to Results/ subdirectory for re-use; the figures can
##   be rebuilt from it with regenerate_figures.R without re-running
##   the sampler.
############################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)
library(purrr)
library(tibble)
library(ggrepel)
library(ggtext)
library(AUCell)
library(UCell)
library(GSVA)
library(Matrix)
library(GenomicRanges)
library(stringr)

## Theme, palettes, notation helpers, save_fig() and the plot-building
## functions shared by every figure in the repository.
source("Github Codes/figure_style.R")

############################################################
## PREREQUISITE CHECK
##
## The ligand-receptor database and the pathway gene sets are built by
## build_lr_database.R and are not rebuilt here. They are first USED in
## Section 3, several minutes into the run: after donor QC, genotype
## conversion and the carrier-frequency filter. Checking for them here rather than there
## turns a late failure that discards that work into an immediate message.
##
## The compiled sampler is checked for the same reason, and specifically for
## the presence of loglik_draws in its output, so that a compiled object that
## does not match the current source is reported here rather than surfacing
## later as an unexplained failure in the diagnostics.
############################################################
.prereq_objects <- c("lr_matrix", "key_cols")
.prereq_missing <- .prereq_objects[!vapply(.prereq_objects, exists, logical(1))]
if (length(.prereq_missing) > 0) {
  stop("Prerequisite object(s) not found: ",
       paste(.prereq_missing, collapse = ", "),
       ".\n  Run first:  source(\"Github Codes/build_lr_database.R\")",
       call. = FALSE)
}

if (!exists("mr_ccc_gibbs")) {
  stop("The compiled sampler is not available.\n",
       "  Run first:  Rcpp::sourceCpp(\"Github Codes/mr_ccc_gibbs.cpp\")",
       call. = FALSE)
}
## loglik_draws is an element of the RETURN VALUE, not an argument, so it
## cannot be checked with formals(). A throwaway fit on random data is used
## instead; at 20 iterations on 40 rows the cost is not measurable.
local({
  probe <- mr_ccc_gibbs(
    matrix(rnorm(40), 40), matrix(rnorm(40), 40), matrix(rnorm(40), 40),
    matrix(rnorm(80), 40), matrix(rnorm(80), 40), matrix(rnorm(40), 40),
    n_iter = 20, burn_in = 5, thin = 1)
  if (!("loglik_draws" %in% names(probe))) {
    stop("The compiled sampler does not return 'loglik_draws', so it does ",
         "not match the current source file.\n",
         "  Recompile:  Rcpp::sourceCpp(\"Github Codes/mr_ccc_gibbs.cpp\")",
         call. = FALSE)
  }
})

############################################################
## USER SETTINGS
## Set Cell1 (sender) and Cell2 (receiver) before running.
##
## Valid values:
##   "BCells", "CD4Cells", "CD8Cells", "NKCells", "MonocytesCells"
############################################################
## Interactively, edit the two defaults below. To sweep several ordered pairs
## non-interactively, set .mrccc_cell1 / .mrccc_cell2 before sourcing this
## file (run_pair.R does exactly that); the defaults then do not apply.
##
## The override is checked rather than assigned unconditionally because
## source()-ing this script would otherwise reset Cell1 and Cell2 on every
## iteration of any wrapper loop, and every pair would silently be analysed
## as NK cells -> monocytes.
Cell1 <- if (exists(".mrccc_cell1")) .mrccc_cell1 else "NKCells"
Cell2 <- if (exists(".mrccc_cell2")) .mrccc_cell2 else "MonocytesCells"
stopifnot(Cell1 != Cell2)
cat("\n=== Ordered pair:", Cell1, "->", Cell2, "===\n")

## ---- Reproducibility -----------------------------------------------------
## The Gibbs sampler draws through R's RNG, so a seed makes every run exactly
## reproducible. Without one, posterior inclusion probabilities differ between
## runs by an amount governed by the Monte Carlo error of the gamma chain,
## which can be large enough to move a triplet across the discovery threshold.
SEED <- 20260921
set.seed(SEED)

## ---- Memory ---------------------------------------------------------------
## This analysis holds a genes x donors expression matrix per cell type and a
## SNPs x donors genotype matrix simultaneously, which is several gigabytes.
## Run it in a FRESH R session: leftover objects from the sensitivity scripts
## (in particular a full numeric genotype matrix) are large enough on their
## own to exhaust the default vector memory limit.
##
## On macOS the limit can be raised if required, e.g.
##   Sys.setenv(R_MAX_VSIZE = "32Gb")   # must be set BEFORE R starts
## or by adding R_MAX_VSIZE=32Gb to ~/.Renviron.
if (length(ls(envir = .GlobalEnv)) > 50L) {
  warning("The global environment already holds ",
          length(ls(envir = .GlobalEnv)), " objects. ",
          "If memory errors occur, restart R and run this script first.")
}

## ---- MCMC settings -------------------------------------------------------
## thin = 1 by default. The PIP is the proportion of retained iterations with
## gamma = 1, so its Monte Carlo precision depends on the number of retained
## draws; thinning discards information that directly reduces that precision.
## Thinning is useful for memory, which is not a constraint here.
##
## The script reports a Monte Carlo standard error (MCSE) for every PIP, so
## chain length can be chosen empirically: increase N_ITER until MCSE is small
## relative to the distance between the PIP and the discovery threshold.
##
## CANONICAL MCMC SETTINGS -- these values define the real-data protocol.
## Every script that refits real-data triplets inherits them from this file
## (leave_one_instrument_out.R, pathway_representation_sensitivity.R,
## convergence_diagnostics.R, posterior_sign_probabilities.R) or restates
## them under a KEEP IN SYNC comment (instrument_overlap_sensitivity.R,
## Tier 2). The simulation scripts use their own chain lengths, for the
## reasons given in the next block and beside each setting there.
##
## Two settings deserve comment, both driven by the Monte Carlo precision
## of the PIP:
##
##   thin = 1. Thinning does not shorten the chain, it only discards
##   recorded draws, so thin = 1 yields ten times as many retained draws as
##   thin = 10 at exactly the same runtime.
##
##   TOTAL DRAWS: four chains of 100000. The inclusion indicator is strongly
##   autocorrelated: its effective sample size is roughly one percent of the
##   retained draws, so a short chain leaves a Monte Carlo standard error on
##   the PIP large enough for a triplet near the 0.5 threshold to have its
##   discovery status decided by Monte Carlo noise rather than by evidence.
##   The total of approximately 400000 draws is chosen so that the MCSE of
##   every PIP is small relative to its distance from the threshold.
##
##   The stopping rule is empirical and is the thing to defend, not the
##   number itself: the script reports MCSE and ESS for every triplet, so
##   the adequacy of the chain length can be verified from the output.
##
##   WHY FOUR CHAINS OF 100000 RATHER THAN ONE OF 400000. The two are
##   equivalent in precision. The Monte Carlo standard error of a chain
##   average scales as the reciprocal square root of its length, so a chain of
##   100000 carries twice the MCSE of a chain of 400000; combining C chains in
##   quadrature divides by C, giving sqrt(C * (2m)^2) / C = 2m / sqrt(C),
##   which at C = 4 returns exactly m, the single-chain value. Four chains of
##   100000 therefore reproduce the precision of one chain of 400000 at the
##   same total cost, and they additionally permit the Gelman-Rubin statistic,
##   which no single chain can supply. Sampling takes about 14 seconds per
##   100,000 sweeps at n = 651, i.e. per chain per triplet, against roughly
##   two minutes of pathway scoring, so the cost is negligible either way.
##
## WHY THIS DIFFERS FROM THE SIMULATION STUDY (20,000 iterations there).
## The requirement follows from what each analysis reports, not from any
## ranking of their importance. Here every PIP is reported individually for
## a named triplet and decides a discovery claim, and the observed PIPs SPAN
## the 0.5 threshold, so Monte Carlo error propagates directly into the
## reported result. In the simulation nothing is reported per replicate:
## operating characteristics are averaged over 100 replicates, and the
## simulated PIPs sit far from the threshold (near 0 under the null, near 1
## under signal), five or more Monte Carlo standard errors away, so chain
## noise essentially never flips a decision. This script reports the MCSE of
## every PIP; the S1-S3 simulation output does not store it, while the S4-S7
## output does.
##
## The misspecification study (S4-S7) is the exception among the
## simulations: there the departures from the assumed model pull the
## inclusion probabilities toward the threshold, so it runs a single chain
## of 400,000, matching the total draws used here. The difference from this
## script is only in how those draws are partitioned: several chains here,
## because that is what makes the Gelman-Rubin statistic available for the
## results that are reported triplet by triplet, and a single chain there,
## because between-chain diagnostics on 10,800 independent fits would be
## reported nowhere.
##
## SMOKE-TEST OVERRIDE. N_ITER may be reduced for a quick end-to-end check by
## setting .mrccc_n_iter before sourcing this file, in the same way as
## .mrccc_cell1 / .mrccc_cell2. The override is announced on every run and the
## value is recorded in the output object, so a reduced run can never be
## mistaken for a reported one. Assigning N_ITER directly before source() has
## no effect, because the assignment below would overwrite it.
N_ITER  <- if (exists(".mrccc_n_iter")) as.integer(.mrccc_n_iter) else 100000L
THIN    <- 1

## The canonical burn-in is 2,000. Under a smoke-test override that would leave
## a short chain with no retained draws at all, so the burn-in is reduced to a
## quarter of the chain in that case only. The canonical path is untouched:
## min(2000, 100000 %/% 4) is 2000.
BURN_IN <- min(2000L, N_ITER %/% 4L)

if (N_ITER != 100000L) {
  cat("\n!!!! NON-CANONICAL CHAIN LENGTH: N_ITER =", N_ITER,
      "(canonical is 100000), BURN_IN =", BURN_IN, "\n",
      "     This run is a check, NOT a reportable result.",
      "Remove .mrccc_n_iter to restore.\n\n")
}
stopifnot(
  "N_ITER must exceed BURN_IN" = N_ITER > BURN_IN,
  "the configuration must retain at least two draws per chain" =
    (N_ITER - BURN_IN) %/% THIN >= 2L
)

## ---- Number of chains and starting values --------------------------------
## The primary analysis runs N_CHAINS chains per triplet from OVERDISPERSED
## starting values and reports the pooled posterior. Convergence is assessed
## by the Gelman-Rubin statistic across those chains, and Monte Carlo
## precision by batch-means MCSE and ESS pooled in quadrature.
##
## Why several chains rather than one long chain. The two are equivalent in
## Monte Carlo precision once the total number of draws is matched: a chain of
## N_ITER has MCSE proportional to N_ITER^(-1/2), and combining C such chains
## in quadrature, sqrt(sum_c mcse_c^2) / C, recovers the MCSE of a single
## chain of C * N_ITER. Four chains of 100,000 and one chain of 400,000 are
## therefore the same precision at the same cost. Several chains additionally
## admit R-hat, which compares between-chain with within-chain variance and is
## the standard convergence diagnostic; a single chain admits only
## within-chain substitutes.
##
## INIT_SCALE MUST BE POSITIVE whenever N_CHAINS > 1. R-hat is a ratio of
## between-chain to within-chain variance, so it can only detect a failure to
## converge if the chains begin at genuinely different points. Started from a
## common value, the between-chain variance is small by construction and R-hat
## is close to one whether or not the chains have converged. INIT_SCALE gives
## the standard deviation of the normal draw used for the initial values of
## beta_X and beta_XZ, and init_gamma alternates between chains so that the
## inclusion indicator is also dispersed.
N_CHAINS   <- 4
INIT_SCALE <- 1.0   # standard deviation of the dispersed start; must be > 0 when N_CHAINS > 1

stopifnot(
  "N_CHAINS must be a positive whole number" =
    is.numeric(N_CHAINS) && length(N_CHAINS) == 1L &&
    N_CHAINS >= 1 && N_CHAINS == round(N_CHAINS),
  "INIT_SCALE must be > 0 when N_CHAINS > 1, or R-hat cannot detect non-convergence" =
    N_CHAINS == 1L || INIT_SCALE > 0
)

## ---- Credible band for the effect curve ----------------------------------
## The receptor-modulated ligand effect, effect(z) = beta_X + beta_XZ * z, is
## reported with a pointwise credible band rather than as a single line
## through the posterior means. A single line carries no uncertainty, which
## matters when comparing sign-reversal thresholds between triplets.
##
## The band is built from PAIRED draws, for the same reason tau is: the curve
## depends on TWO parameters at once, so combining their marginal summaries
## would misstate the uncertainty wherever they are correlated -- and here
## they are strongly correlated.
##
## The pooled draws number roughly 400,000. A systematic subsample of
## BAND_DRAWS is far more than enough to place a 2.5%/97.5% quantile, and
## keeps the band a trivial calculation rather than a 250-million-element one.
BAND_DRAWS <- 2000L

## ---- Discovery rule ------------------------------------------------------
## "pip"  : PIP > pip_thresh (median probability model; Barbieri & Berger 2004)
## "bfdr" : largest set whose Bayesian false discovery rate does not exceed
##          bfdr_alpha, where  bFDR(S) = |S|^-1 * sum_{j in S} (1 - PIP_j)
##
## Both are always computed and reported; DISCOVERY_RULE only determines which
## one defines `discover` for the figures.
##
## NOTE ON THE DEFAULT. "pip" is kept as the primary rule because the median
## probability model is a standard, citable default, not because it is the
## more permissive of the two. Bayesian FDR is computed for every triplet and
## reported alongside, so the error rate incurred by the primary rule is
## stated explicitly rather than left implicit. This is the honest ordering:
## report a principled rule, then disclose its error rate, then show the
## FDR-controlled subsets as a nested hierarchy of confidence.
##
## Switching DISCOVERY_RULE to "bfdr" is legitimate but markedly more
## stringent at these PIP levels, and will retain far fewer triplets.
DISCOVERY_RULE <- "pip"
pip_thresh     <- 0.5     # used by the "pip" rule and always reported
bfdr_alpha     <- 0.10    # target Bayesian FDR for the "bfdr" rule

## ---- Instrument selection: linkage-disequilibrium clumping ----------------
## Instruments are the ten cis SNPs with the smallest marginal eQTL p-value,
## taken greedily in p-value order and skipping any SNP whose squared
## dosage correlation with an instrument already kept exceeds LD_R2_MAX.
##
## Without this step the top ten by p-value are often members of one
## linkage block, and a set can contain pairs of SNPs with |r| close to one,
## or effectively identical SNPs. Such sets report ten instruments while
## carrying the information of a few, which dilutes the per-instrument
## first-stage F and places the design outside the correlation range the
## simulations examine.
##
## The threshold r^2 = 0.8, i.e. |r| = 0.9, is the largest correlation the
## correlated-instrument simulation (S5) covers, so every retained set lies
## within the range for which the method's behaviour has been characterised.
## Clumping keeps the count at ten where enough independent SNPs exist, so
## redundant tags are replaced by further independent signal rather than
## simply removed. The lead SNP is always kept, so no gene loses all its
## instruments and the triplet universe is unchanged.
##
## Setting LD_R2_MAX to 1 disables clumping (the ten SNPs with the smallest
## p-values are taken). It can be overridden without editing this file by
## assigning .mrccc_ld_r2_max before sourcing, as with the other .mrccc_*
## settings.
##
## KEEP IN SYNC with LD_R2_MAX in instrument_overlap_sensitivity.R.
LD_R2_MAX <- if (exists(".mrccc_ld_r2_max")) as.numeric(.mrccc_ld_r2_max) else 0.8
stopifnot(is.numeric(LD_R2_MAX), LD_R2_MAX > 0, LD_R2_MAX <= 1)

# Root directory where the two Zenodo .rda files are saved.
#
# The default assumes ~/data/mrccc. To point at a different local copy
# without editing this file (and therefore without any risk of committing
# a machine-specific path), set an environment variable instead:
#
#   Sys.setenv(MRCCC_DATA_DIR = "/path/to/your/data")   # this session only
#
# or add a permanent line to ~/.Renviron:
#
#   MRCCC_DATA_DIR=/path/to/your/data
DATA_DIR <- Sys.getenv("MRCCC_DATA_DIR",
                       unset = file.path("~", "data", "mrccc"))
cat("DATA_DIR:", DATA_DIR, "\n")

############################################################
## LOAD DATA
############################################################
# Single-cell expression data: pseudo-bulk count matrices
# (genes x donors) for B cells, CD4+ T cells, CD8+ T cells,
# NK cells, and monocytes from the OneK1K cohort.
# Download from: https://doi.org/10.5281/zenodo.19675075
load(file.path(DATA_DIR, "B_T_NK_monocytes.rda"))

# Donor-level data from OneK1K.
# Contains:
#   donor    -- metadata (age, sex, ancestry PCs)
#   genotype -- character genotype matrix (SNP x donor)
#   vcf.ref  -- GRanges: SNP IDs and genomic coordinates
#   gene.ref -- GRanges: gene IDs, symbols, and coordinates
# Download from: https://doi.org/10.5281/zenodo.19675075
load(file.path(DATA_DIR, "donor.rda"))

############################################################
## SECTION 1: AGGREGATE CELL-LEVEL COUNTS TO DONOR LEVEL
##
## Each object (B_Cells, CD4_Cells, etc.) is a list of
## per-donor matrices (genes x cells). colSums aggregates
## across cells to give total read counts per gene per donor.
## Monocytes uses a helper because its list elements may
## have varying structure.
############################################################
RNA.count_BCells   <- sapply(B_Cells,   colSums)
RNA.count_CD4Cells <- sapply(CD4_Cells, colSums)
RNA.count_CD8Cells <- sapply(CD8_Cells, colSums)
RNA.count_NKCells  <- sapply(NK_Cells,  colSums)

# Monocytes list elements can be either a matrix or a named
# numeric vector; this helper handles both cases.
get_gene_counts <- function(x) {
  if (!is.null(dim(x)) && length(dim(x)) == 2) return(colSums(x))
  if (is.atomic(x) && is.null(dim(x)))           return(x)
  stop("Unexpected object in monocytes: class = ",
       paste(class(x), collapse = "+"))
}
RNA.count_MonocytesCells <- sapply(monocytes, get_gene_counts)

# The raw per-cell lists are the largest objects in the session (they hold one
# genes x cells matrix per donor) and nothing downstream uses them once the
# donor-level matrices exist. Releasing them here frees several gigabytes and
# is what allows the analysis to run alongside the genotype matrices without
# exhausting the vector memory limit.
rm(B_Cells, CD4_Cells, CD8_Cells, NK_Cells, monocytes)
invisible(gc())

# Select the count matrices for the active cell-type pair.
RNA.count_Cell1 <- get(paste0("RNA.count_", Cell1))
RNA.count_Cell2 <- get(paste0("RNA.count_", Cell2))

############################################################
## SECTION 2: QUALITY CONTROL AND LIBRARY-SIZE NORMALISATION
##
## Library size = total reads per donor, normalized to the
## median across donors. Donors with very low or very high
## library size (> 3 MADs above the median, or < 0.25)
## are excluded as outliers. After filtering, each count
## matrix is divided by its donor's library-size factor so
## that all donors are on a comparable scale. These normalised
## counts are used for every downstream quantity (ligand X, receptor
## Z, the eQTL screen and all six pathway scores).
############################################################
# Donors are aligned POSITIONALLY across the expression matrices, the
# genotype matrix and the donor table (the expression matrices carry no
# donor names). Equal donor counts are therefore the contract the whole
# pipeline relies on, and are checked before any subsetting.
stopifnot(
  "expression matrices and donor table differ in donor count" =
    ncol(RNA.count_Cell1) == nrow(donor) && ncol(RNA.count_Cell2) == nrow(donor),
  "genotype matrix and donor table differ in donor count" =
    ncol(genotype) == nrow(donor)
)

# Compute per-donor library-size factors
lib.size1 <- colSums(RNA.count_Cell1) / median(colSums(RNA.count_Cell1))
lib.size2 <- colSums(RNA.count_Cell2) / median(colSums(RNA.count_Cell2))

# Identify donors passing QC in both cell types
donor.filter1 <- (lib.size1 <= median(lib.size1) + 3 * mad(lib.size1)) &
  (lib.size1 > 0.25)
donor.filter2 <- (lib.size2 <= median(lib.size2) + 3 * mad(lib.size2)) &
  (lib.size2 > 0.25)
donor.filter  <- (donor.filter1 * donor.filter2 != 0)

# Apply filter to all donor-level objects
donor      <- donor[donor.filter, ]
genotype   <- genotype[, donor.filter]
RNA.count_Cell1_Filt <- RNA.count_Cell1[, donor.filter]
RNA.count_Cell2_Filt <- RNA.count_Cell2[, donor.filter]
rm(donor.filter)

# Recompute library sizes on filtered donors and normalise
lib.size1 <- colSums(RNA.count_Cell1_Filt) /
  median(colSums(RNA.count_Cell1_Filt))
RNA.count.adj_Cell1 <- sweep(RNA.count_Cell1_Filt, 2, lib.size1, "/")

lib.size2 <- colSums(RNA.count_Cell2_Filt) /
  median(colSums(RNA.count_Cell2_Filt))
RNA.count.adj_Cell2 <- sweep(RNA.count_Cell2_Filt, 2, lib.size2, "/")

rm(RNA.count_Cell1, RNA.count_Cell2,
   RNA.count_Cell1_Filt, RNA.count_Cell2_Filt)
invisible(gc())
cat("After QC:", ncol(RNA.count.adj_Cell1), "donors retained.\n")


# ---- Donor-level covariates V -------------------------------------------
# One covariate matrix, used everywhere a covariate adjustment is made: the
# eQTL screen, the first-stage F statistic and the model itself. Age is
# binned into 10-year intervals (1 = <30, ..., 7 = 80+) to regularise its
# effect given the moderate donor count.
#
# KEEP IN SYNC with build_covariates() in instrument_overlap_sensitivity.R.
build_covariates <- function(donor_df) {
  age <- donor_df$age
  cbind(
    MALE = as.integer(donor_df$sex == "male"),
    AGE  = ifelse(age < 30, 1L, ifelse(age < 40, 2L, ifelse(age < 50, 3L,
           ifelse(age < 60, 4L, ifelse(age < 70, 5L, ifelse(age < 80, 6L, 7L)))))),
    PC1  = donor_df$PC1, PC2 = donor_df$PC2, PC3 = donor_df$PC3
  )
}
V_donor <- build_covariates(donor)
if (anyNA(V_donor)) {
  stop("Missing donor covariates (age, sex or PC1-PC3) for ",
       sum(!stats::complete.cases(V_donor)), " retained donors.")
}

# ---- Column centring ------------------------------------------------------
# The model is written for centred variables and has no intercept in either
# first stage, so every column of G, H and V is centred before it reaches the
# sampler, exactly as X, Z and Y are. This makes the projections X* and Z*
# mean-zero for every draw, so that beta_X is the ligand effect at the
# average receptor level, and keeps the plug-in prior scale free of
# contributions from column means.
centre_cols <- function(M) {
  M <- as.matrix(M)
  sweep(M, 2, colMeans(M), "-")
}

############################################################
## SECTION 3: GENOTYPE PROCESSING
##
## The raw genotype matrix uses string encoding ("0/0", "0/1",
## "1/1", "./."). Convert to a numeric dosage matrix (0/1/2),
## with missing calls ("./.") left missing, then keep SNPs
## carried (at least one copy) by more than 5% of donors with
## a called genotype. Missing calls are then imputed by the
## SNP's mean dosage over the retained donors.
############################################################
# Promote gene reference to promoter windows (+/-200 kb)
# used for cis-eQTL instrument selection.
gene.promoter.ref <- promoters(gene.ref,
                               upstream   = 200000,
                               downstream = 200000)

# Convert string genotypes to numeric dosage; "./." stays NA
genotype.mat <- matrix(NA_real_, nrow = nrow(genotype), ncol = ncol(genotype))
genotype.mat[genotype == "0/0"] <- 0
genotype.mat[genotype == "0/1"] <- 1
genotype.mat[genotype == "1/1"] <- 2

# Carrier-frequency filter over called genotypes: keep SNPs with > 5% of
# donors carrying at least one copy. SNPs with no called genotype at all are
# dropped.
SNP.filter <- rowMeans(genotype.mat > 0, na.rm = TRUE) > 0.05
SNP.filter[is.na(SNP.filter)] <- FALSE
cat("SNPs retained after carrier-frequency filter:", sum(SNP.filter), "\n")
genotype     <- genotype[SNP.filter, ]
genotype.mat <- genotype.mat[SNP.filter, , drop = FALSE]
vcf.ref      <- vcf.ref[SNP.filter]

# Mean imputation of missing calls
.na_geno <- which(is.na(genotype.mat), arr.ind = TRUE)
if (nrow(.na_geno) > 0) {
  .snp_mean <- rowMeans(genotype.mat, na.rm = TRUE)
  genotype.mat[.na_geno] <- .snp_mean[.na_geno[, 1]]
}
cat("Missing genotype calls imputed by SNP mean:", nrow(.na_geno), "\n")
rm(.na_geno)

############################################################
## SECTION 4: FILTER LR PAIRS TO GENES PRESENT IN DATA
##
## lr_matrix (from build_lr_database.R) contains all
## CellPhoneDB Ligand-Receptor pairs with Ensembl IDs.
## We keep only pairs where both genes appear as rows of
## RNA.count_BCells, whose row names serve as the gene
## universe for every cell-type pair. key_cols (MSigDB
## Reactome) is filtered against the same gene list.
############################################################
valid_genes <- rownames(RNA.count_BCells)  # shared gene universe

lr_filtered <- lr_matrix %>%
  filter(ligand_ensembl   %in% valid_genes,
         receptor_ensembl %in% valid_genes)

key_cols <- key_cols %>%
  filter(ensembl_gene %in% valid_genes)

cat("LR pairs retained after gene filter:", nrow(lr_filtered), "\n")
cat("Pathways retained after gene filter:",
    n_distinct(lr_filtered$pathway_name), "\n")

## SMOKE-TEST OVERRIDE. Setting .mrccc_max_triplets truncates the triplet
## universe, which is the only way to make an end-to-end check quick: the cost
## per triplet is dominated by pathway scoring (roughly two minutes) rather
## than by sampling, so shortening the chain alone barely shortens the run.
## Truncation changes the reported result completely -- the discovery set, the
## Bayesian FDR and the selection-frequency table are all defined over the full
## universe -- so it is announced here and again in the Section 6a summary.
.mrccc_truncated <- FALSE
if (exists(".mrccc_max_triplets")) {
  .n_full <- nrow(lr_filtered)
  .keep_n <- min(as.integer(.mrccc_max_triplets), .n_full)
  lr_filtered <- lr_filtered[seq_len(.keep_n), , drop = FALSE]
  .mrccc_truncated <- TRUE
  cat("\n!!!! TRUNCATED TRIPLET UNIVERSE:", .keep_n, "of", .n_full,
      "triplets.\n",
      "     Discovery counts and bFDR are NOT comparable with a full run.",
      "Remove .mrccc_max_triplets to restore.\n\n")
}

############################################################
## SECTION 5: HELPER FUNCTIONS
############################################################
# Greedy LD clumping. `cand` holds genotype row indices already sorted by
# eQTL p-value, strongest first. Walks the list keeping a SNP only if its
# squared dosage correlation with every SNP already kept is at most r2_max,
# and stops once keep_max are kept. Returns the kept indices in p-value
# order.
#
# Because each decision depends only on SNPs earlier in the list, the first
# k kept for any k <= keep_max are exactly what a run with keep_max = k
# would keep; the instrument-count sensitivity relies on this.
#
# A SNP with no variation in the retained donors has an undefined
# correlation; it is skipped rather than kept, since it cannot serve as an
# instrument. r2_max = 1 keeps every SNP and reproduces the unclumped top-k.
#
# KEEP IN SYNC with clump_by_ld() in instrument_overlap_sensitivity.R.
clump_by_ld <- function(cand, geno, keep_max, r2_max) {
  if (r2_max >= 1) return(cand[seq_len(min(keep_max, length(cand)))])
  kept <- integer(0)
  for (s in cand) {
    if (length(kept) >= keep_max) break
    if (length(kept) == 0L) { kept <- s; next }
    r <- suppressWarnings(stats::cor(geno[s, ], t(geno[kept, , drop = FALSE])))
    if (any(is.na(r))) next
    if (max(r^2) <= r2_max) kept <- c(kept, s)
  }
  kept
}

# Build the cis-eQTL instrument matrix G for a given gene.
# Scans the gene's promoter window (+/-200 kb) for SNPs, ranks them by eQTL
# p-value, and returns up to 10 as columns of the instrument matrix, taken
# in p-value order with LD clumping at LD_R2_MAX (see USER SETTINGS).
# Returns NULL if fewer than min_snps qualifying SNPs exist.
get_SNP_matrix <- function(gene_id, RNA_mat_adj, min_snps = 5) {
  gene.ind <- which(gene.ref$gene_id == gene_id)
  if (length(gene.ind) == 0) return(NULL)
  
  snp.ind <- which(
    countOverlaps(vcf.ref, gene.promoter.ref[gene.ind]) > 0
  )
  if (length(snp.ind) < min_snps) return(NULL)
  
  # Compute marginal eQTL p-value for each candidate SNP, adjusting for the
  # same covariates V that enter the model
  snp.pval <- rep(NA_real_, length(snp.ind))
  for (j in seq_along(snp.ind)) {
    g_j <- genotype.mat[snp.ind[j], ]
    # A SNP with no variation among these donors has no eQTL p-value; lm()
    # would drop its term and row 2 would belong to a covariate.
    if (var(g_j) == 0) next
    lm.j <- lm(RNA_mat_adj[gene_id, ] ~ g_j + V_donor)
    cf <- summary(lm.j)$coefficients
    snp.pval[j] <- if ("g_j" %in% rownames(cf)) cf["g_j", 4] else NA_real_
  }
  
  ok <- which(!is.na(snp.pval))
  if (length(ok) < min_snps) return(NULL)
  
  # Up to 10 SNPs in eQTL p-value order, LD-clumped at LD_R2_MAX
  ranked    <- snp.ind[ok[order(snp.pval[ok])]]
  best_snps <- clump_by_ld(ranked, genotype.mat, keep_max = 10L,
                           r2_max = LD_R2_MAX)
  t(genotype.mat[best_snps, , drop = FALSE])  # n x k (donors x SNPs)
}

# ---- Orientation of the PC1 pathway score -----------------------------------
# A principal component is defined only up to sign. Orienting it against a
# reference score (the mean expression of the pathway genes) makes the sign of
# the reported effects reproducible and interpretable: higher PC1 means higher
# pathway expression.
#
# KEEP IN SYNC with the copies in pathway_representation_sensitivity.R and
# instrument_overlap_sensitivity.R.
orient_pc1 <- function(pc, ref) {
  r <- suppressWarnings(stats::cor(pc, ref, use = "complete.obs"))
  if (is.finite(r) && r < 0) -pc else pc
}

# ---- GSVA / ssGSEA scores, tolerant of both package APIs -------------------
# GSVA changed its interface at version 1.52 (Bioconductor 3.18). The former
# call
#     gsva(expr = ., gset.idx.list = ., method = "gsva"|"ssgsea", ...)
# was replaced by parameter objects
#     gsva(gsvaParam(exprData = ., geneSets = ., kcdf = .))
#     gsva(ssgseaParam(exprData = ., geneSets = .))
# and the old signature now fails with
#     "unable to find an inherited method for function 'gsva'".
#
# Both forms are supported here so that the code runs unchanged on older and
# newer installations. Note that `abs.ranking` applied only to the
# "gsva" method in the legacy interface and was silently ignored for ssGSEA;
# the new interface names it `absRanking` and it is left at its default.
#
# KEEP IN SYNC with the copy in instrument_overlap_sensitivity.R, which
# defines it only when this one is not already present. The other refit
# scripts inherit this definition.
gsva_scores <- function(expr_full, geneSets,
                        method = c("gsva", "ssgsea"),
                        kcdf = "Gaussian") {
  method <- match.arg(method)
  has_new_api <- exists("gsvaParam", where = asNamespace("GSVA"),
                        inherits = FALSE)
  if (has_new_api) {
    prm <- if (method == "gsva") {
      GSVA::gsvaParam(exprData = expr_full, geneSets = geneSets, kcdf = kcdf)
    } else {
      GSVA::ssgseaParam(exprData = expr_full, geneSets = geneSets)
    }
    as.matrix(GSVA::gsva(prm, verbose = FALSE))
  } else {
    as.matrix(GSVA::gsva(expr = expr_full, gset.idx.list = geneSets,
                         method = method, kcdf = kcdf, verbose = FALSE))
  }
}

# ---- First-stage F-statistic ----------------------------------------------
# Partial F comparing the covariate-only model with the model that also
# contains all instruments. This is the conventional diagnostic for
# instrument strength.
#
# Two cautions on interpretation:
#   * the familiar "F > 10" rule of thumb derives from two-stage least squares
#     with a single exposure and has no established analogue for the
#     interaction term used here, whose design column X* o Z* is a PRODUCT of
#     two first-stage projections and therefore inherits weakness from both;
#   * instruments are selected by p-value in these same data, so F is
#     optimistic (winner's curse).
# F is therefore reported as a diagnostic, not as a certificate.
#
# The covariates are the model's own V (see build_covariates()), so that the
# F statistic describes the same first stage the sampler fits.
first_stage_F <- function(Gmat, y, Vmat) {
  if (is.null(Gmat) || ncol(Gmat) == 0) return(NA_real_)
  Gm <- as.matrix(Gmat)
  colnames(Gm) <- paste0("snp", seq_len(ncol(Gm)))
  Vm <- as.matrix(Vmat)
  colnames(Vm) <- paste0("cov", seq_len(ncol(Vm)))
  dat <- data.frame(y = as.numeric(y), Vm, Gm, check.names = FALSE)
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  if (nrow(dat) <= ncol(dat) + 1L) return(NA_real_)
  m0 <- lm(y ~ ., data = dat[, c("y", colnames(Vm))])
  m1 <- lm(y ~ ., data = dat)
  a  <- tryCatch(anova(m0, m1), error = function(e) NULL)
  if (!is.null(a) && nrow(a) >= 2 && "F" %in% colnames(a)) a[2, "F"] else NA_real_
}

# ---- Monte Carlo error of an MCMC average ---------------------------------
# Non-overlapping batch-means estimator. Returns the Monte Carlo standard
# error and the implied effective sample size.
#
# For the PIP this is the quantity that matters: the PIP is a chain average of
# a binary indicator, and the gamma chain is autocorrelated, so the naive
# binomial standard error badly understates the true uncertainty.
mcse_batch <- function(x, n_batch = 30L) {
  n <- length(x)
  b <- floor(n / n_batch)
  if (n < 4L || b < 2L) return(c(mcse = NA_real_, ess = NA_real_))
  k  <- floor(n / b)
  bm <- vapply(seq_len(k), function(i) mean(x[((i - 1) * b + 1):(i * b)]),
               numeric(1))
  sigma2 <- b * stats::var(bm)          # asymptotic variance of the chain mean
  vx     <- stats::var(x)
  if (!is.finite(sigma2) || sigma2 <= 0) return(c(mcse = NA_real_, ess = NA_real_))
  c(mcse = sqrt(sigma2 / n), ess = n * vx / sigma2)
}

# ---- Gelman-Rubin R-hat ----------------------------------------------------
# R-hat for a list of equal-length chains (Gelman & Rubin, 1992). Values near
# 1 indicate that the between-chain and within-chain variances agree; a common
# working rule is R-hat < 1.01.
#
# This is character-for-character the function used in
# convergence_diagnostics.R. The two must not be allowed to drift apart: the
# diagnostics script exists to certify the numbers this script reports, which
# it can only do if both compute R-hat the same way.
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

# ---- Monte Carlo error of an average pooled over independent chains --------
# Given the per-chain batch-means standard errors, the pooled mean is the
# average of the chain means, so its variance is the average of the per-chain
# variances divided by the number of chains:
#
#   mcse_pooled = sqrt( sum_c mcse_c^2 ) / C
#
# Effective sample sizes add, because the chains are independent.
#
# Applying mcse_batch() to the concatenated draws instead would place batches
# across the joins between chains, which started from different values, and
# would inflate the estimated variance.
mcse_pooled <- function(mc_list) {
  ok <- vapply(mc_list, function(z) is.finite(z["mcse"]), logical(1))
  if (!any(ok)) return(c(mcse = NA_real_, ess = NA_real_))
  se  <- vapply(mc_list[ok], function(z) unname(z["mcse"]), numeric(1))
  es  <- vapply(mc_list[ok], function(z) unname(z["ess"]),  numeric(1))
  # Divide by the TOTAL number of chains, not the number with a finite MCSE.
  # mcse_batch() returns NA for a chain that is constant over all retained
  # draws (PIP at exactly 0 or 1), whose true MCSE is zero; treating it as
  # absent instead of as zero would overstate the pooled MCSE by up to C-fold.
  c(mcse = sqrt(sum(se^2)) / length(mc_list), ess = sum(es, na.rm = TRUE))
}

# ---- Bayesian false discovery rate ----------------------------------------
# For a declared set S, the posterior expected number of false discoveries is
# sum_{j in S} (1 - PIP_j), so the Bayesian FDR is that quantity divided by
# |S|. Ranking triplets by PIP and accumulating gives, for each possible
# cut-off, the FDR incurred by declaring everything above it.
#
# Because PIPs are posterior probabilities, this controls the error rate
# without a separate multiplicity correction.
bayes_fdr <- function(pip) {
  o   <- order(pip, decreasing = TRUE)
  f   <- cumsum(1 - pip[o]) / seq_along(pip)
  out <- numeric(length(pip))
  out[o] <- f
  out
}

# Largest set whose Bayesian FDR does not exceed alpha (a nested, monotone
# rule: if no set qualifies, nothing is declared).
bfdr_threshold_set <- function(pip, alpha) {
  o <- order(pip, decreasing = TRUE)
  f <- cumsum(1 - pip[o]) / seq_along(pip)
  k <- suppressWarnings(max(which(f <= alpha)))
  keep <- rep(FALSE, length(pip))
  if (is.finite(k) && k >= 1) keep[o[seq_len(k)]] <- TRUE
  keep
}

############################################################
## SECTION 6: MAIN LOOP OVER LR-PATHWAY TRIPLETS
##
## For each (ligand, receptor, pathway) triplet in lr_filtered:
##   1. Construct Y  -- pathway activity score for Cell2
##   2. Construct X  -- ligand expression in Cell1
##   3. Construct Z  -- receptor expression in Cell2
##   4. Construct G  -- sender cis-eQTL instruments
##   5. Construct H  -- receiver cis-eQTL instruments
##   6. Construct V  -- shared donor covariates
##   7. Run MR-CCC Gibbs sampler
##   8. Standardize effect sizes and store output row
############################################################
Output_list <- vector("list", nrow(lr_filtered))
out_idx     <- 0

for (run in seq_len(nrow(lr_filtered))) {
  
  # Progress. The sampler itself is silent (verbose = FALSE below), because a
  # line every 1000 iterations across N_CHAINS chains and every triplet would
  # run to tens of thousands of lines. One line per triplet, with the elapsed
  # time, is enough to see where an hour-long run has got to and to notice a
  # triplet that is taking far longer than its neighbours.
  if (run == 1L) .t0_run <- Sys.time()
  cat(sprintf("[%3d/%3d] %-12s %-12s %-28s  %5.1f min elapsed\n",
              run, nrow(lr_filtered),
              lr_filtered$ligand_symbol[run],
              lr_filtered$receptor_symbol[run],
              substr(lr_filtered$pathway_name[run], 1L, 28L),
              as.numeric(difftime(Sys.time(), .t0_run, units = "mins"))))
  utils::flush.console()

  pathway_name       <- lr_filtered$pathway_name[run]
  ligand_gene_name   <- lr_filtered$ligand_symbol[run]
  ligand_gene_id     <- lr_filtered$ligand_ensembl[run]
  receptor_gene_name <- lr_filtered$receptor_symbol[run]
  receptor_gene_id   <- lr_filtered$receptor_ensembl[run]
  
  # ---- Step 1: Build pathway activity Y for Cell2 ----
  # Extract keyword from pathway name to find corresponding
  # genes in MSigDB Reactome (key_cols).
  keyword          <- sub(".*by\\s+", "", pathway_name)
  # An empty keyword would turn the regex below into "", which matches EVERY
  # gene set and would silently build Y from the whole transcriptome. Skip
  # such rows here, exactly as convergence_diagnostics.R and
  # pathway_representation_sensitivity.R already do, so that every triplet
  # the primary emits can be reproduced by those scripts.
  if (!nzchar(keyword)) next
  hits             <- key_cols %>%
    dplyr::filter(str_detect(gs_name, regex(keyword, ignore_case = TRUE)))
  pathway_gene_ids <- unique(hits$ensembl_gene)
  pathway_gene_ids <- pathway_gene_ids[
    pathway_gene_ids %in% rownames(RNA.count.adj_Cell2)
  ]
  if (length(pathway_gene_ids) == 0) next
  
  # ---- Step 2: X and Z (ligand and receptor expression) ----
  if (!ligand_gene_id   %in% rownames(RNA.count.adj_Cell1)) next
  if (!receptor_gene_id %in% rownames(RNA.count.adj_Cell2)) next
  
  X <- matrix(RNA.count.adj_Cell1[ligand_gene_id,   ], ncol = 1)
  Z <- matrix(RNA.count.adj_Cell2[receptor_gene_id, ], ncol = 1)
  
  # ---- Step 3: Compute six pathway activity scores for Y ----
  # Six scoring methods are computed and the one most
  # correlated with X is selected as the outcome variable.
  # This avoids committing to a single method without
  # biological justification.
  expr      <- RNA.count.adj_Cell2[pathway_gene_ids, , drop = FALSE]
  expr_full <- RNA.count.adj_Cell2
  
  # (a) PC1 of the pathway gene expression matrix
  pc_fit <- prcomp(t(expr), center = TRUE, scale. = TRUE)
  Y_PC1  <- pc_fit$x[, 1]
  
  # (b) AUCell: area under the recovery curve
  geneRanks  <- AUCell_buildRankings(expr_full, plotStats = FALSE)
  path_genes <- intersect(rownames(expr_full), pathway_gene_ids)
  geneSets   <- list(pathway = path_genes)
  nGenes     <- nrow(expr_full)
  auc        <- AUCell_calcAUC(geneSets, geneRanks,
                               aucMaxRank = ceiling(0.05 * nGenes))
  Y_AUCell   <- as.numeric(getAUC(auc)["pathway", ])
  
  # (c) Simple mean of pathway gene expression
  Y_rowMean <- colMeans(expr)

  # The sign of a principal component is arbitrary (prcomp may return PC1 or
  # -PC1), and flipping Y flips the signs of beta_X and beta_XZ. PC1 is
  # therefore oriented so that higher scores mean higher average expression
  # of the pathway genes.
  Y_PC1 <- orient_pc1(Y_PC1, Y_rowMean)

  # (d) UCell: Mann-Whitney U-based enrichment score
  ucell_scores <- ScoreSignatures_UCell(
    expr_full, features = list(pathway = path_genes)
  )
  Y_UCell <- as.numeric(ucell_scores[, "pathway_UCell"])
  
  # (e) ssGSEA: single-sample GSEA enrichment score
  Y_ssGSEA <- as.numeric(
    gsva_scores(expr_full, geneSets, method = "ssgsea")["pathway", ])

  # (f) GSVA: Gaussian kernel GSVA score
  Y_GSVA <- as.numeric(
    gsva_scores(expr_full, geneSets, method = "gsva", kcdf = "Gaussian")["pathway", ])
  
  # Select Y with highest absolute correlation with X
  Y_list    <- list(PC1 = Y_PC1, AUCell = Y_AUCell, Mean = Y_rowMean,
                    UCell = Y_UCell, ssGSEA = Y_ssGSEA, GSVA = Y_GSVA)
  cors      <- sapply(Y_list, function(y) cor(y, X, use = "complete.obs"))
  best_name <- names(which.max(abs(cors)))
  Y         <- matrix(Y_list[[best_name]], ncol = 1)

  # Record WHICH representation won and how decisively. Recording best_name
  # makes the adaptive rule auditable: the output states which of the six
  # representations was used for every triplet.
  #
  # Two uses depend on this. Characterising how the "best of six" choice
  # behaves, and what it does to calibration, requires the selection
  # frequencies. And any script that reproduces a triplet outside this loop
  # must pick the same representation, which can only be checked against a
  # recorded value.
  #
  # best_cor_gap is the margin over the runner-up: a gap near zero means the
  # choice was nearly arbitrary and the result should be insensitive to it.
  .abs_cors    <- sort(abs(cors), decreasing = TRUE)
  best_cor     <- unname(cors[best_name])
  best_cor_gap <- if (length(.abs_cors) >= 2L) {
    unname(.abs_cors[1] - .abs_cors[2])
  } else NA_real_
  
  # ---- Step 4: Standardization scales (pre-centering) ----
  # SDs are unchanged by centering and are used to convert
  # posterior means to interpretable standardized effect sizes.
  sd_X <- sd(as.numeric(X))
  sd_Z <- sd(as.numeric(Z))
  sd_Y <- sd(as.numeric(Y))
  
  # ---- Centering: enforce E[X] = E[Z] = E[Y] = 0 ----
  # Required by the centering assumption in Proposition 1. G, H and V are
  # centred below in the same way.
  X_c <- matrix(as.numeric(X) - mean(X), ncol = 1)
  Z_c <- matrix(as.numeric(Z) - mean(Z), ncol = 1)
  Y_c <- matrix(as.numeric(Y) - mean(Y), ncol = 1)

  # ---- Step 5: Instrument matrices G and H ----
  # get_SNP_matrix() returns raw dosages, which are what the clumping and the
  # F statistic use; the sampler receives the centred columns.
  G_raw <- get_SNP_matrix(ligand_gene_id,   RNA.count.adj_Cell1, min_snps = 1)
  if (is.null(G_raw)) next
  H_raw <- get_SNP_matrix(receptor_gene_id, RNA.count.adj_Cell2, min_snps = 1)
  if (is.null(H_raw)) next
  G <- centre_cols(G_raw)
  H <- centre_cols(H_raw)

  # ---- Step 6: Covariate matrix V ----
  # The donor-level covariates built once in Section 2, centred.
  V <- centre_cols(V_donor)

  # ---- Step 7: Run MR-CCC Gibbs sampler ------------------------------------
  # g-prior scales follow the paper: min(n, 100).
  # Default hyperparameters: a_sigma=3, b_sigma=2,
  # a_rho=3, b_rho=1, nu1=1e-4.
  #
  # The chains' STARTING VALUES are dispersed, which is what makes the
  # Gelman-Rubin statistic informative: init_gamma alternates 1, 0, 1, 0 so
  # the inclusion indicator begins in the slab for odd chains and in the spike
  # for even ones, and the two causal coefficients are drawn from
  # N(0, INIT_SCALE^2). Their random streams differ because the chains run
  # consecutively through R's RNG. No per-chain seed is set, so the single
  # set.seed(SEED) at the top of the script keeps the whole analysis
  # reproducible.
  #
  # The loop remains correct at N_CHAINS = 1, where it reduces to one chain
  # from init_gamma = 1, but R-hat is then undefined and reported as NA. The
  # assertion beside INIT_SCALE above prevents the combination that would be
  # misleading: several chains from a common start.
  n       <- nrow(X_c)
  g_scale <- min(n, 100.0)

  fits <- lapply(seq_len(N_CHAINS), function(cc) {
    mr_ccc_gibbs(
      X_c, Z_c, Y_c, G, H, V,
      n_iter  = N_ITER, burn_in = BURN_IN, thin = THIN,
      a_sigma = 3.0,   b_sigma = 2.0,
      a_rho   = 3.0,   b_rho   = 1.0,
      nu1     = 1e-4,
      gG      = g_scale, gH = g_scale,
      gV      = g_scale, gZ = g_scale,
      gBeta   = g_scale,
      ridge   = 1e-8,
      init_gamma = if (cc %% 2L == 0L) 0L else 1L,
      init_scale = INIT_SCALE,
      verbose    = FALSE
    )
  })

  # ---- Step 8: Standardize effect sizes ----
  # Convert from raw expression units to SD units:
  #   beta_X^(s)  = beta_X  * sd(X) / sd(Y)
  #   beta_XZ^(s) = beta_XZ * sd(X) * sd(Z) / sd(Y)
  #   beta_Z^(s)  = beta_Z  * sd(Z) / sd(Y)
  #
  # The rescaling is applied to the retained DRAWS, so that credible
  # intervals are reported on the interpretable standardized scale and the
  # posterior means below are simple averages of those rescaled draws.
  sx_y  <- sd_X / sd_Y
  sxz_y <- sd_X * sd_Z / sd_Y
  sz_y  <- sd_Z / sd_Y

  # Draws kept SEPARATELY BY CHAIN (needed for R-hat) and POOLED (needed for
  # the posterior summaries). All chains have the same length, so the pooled
  # average equals the average of the chain averages.
  bX_ch  <- lapply(fits, function(f) as.numeric(f$Beta_X_draws)  * sx_y)
  bXZ_ch <- lapply(fits, function(f) as.numeric(f$Beta_XZ_draws) * sxz_y)
  bZ_ch  <- lapply(fits, function(f) as.numeric(f$Beta_Z_draws)  * sz_y)
  g_ch   <- lapply(fits, function(f) as.numeric(f$gamma_draws))
  # The joint log-likelihood is on its own scale and is never standardised.
  ll_ch  <- lapply(fits, function(f) as.numeric(f$loglik_draws))

  bX_dr  <- unlist(bX_ch,  use.names = FALSE)
  bXZ_dr <- unlist(bXZ_ch, use.names = FALSE)
  g_dr   <- unlist(g_ch,   use.names = FALSE)

  Beta_X_std  <- mean(bX_dr)
  Beta_XZ_std <- mean(bXZ_dr)
  Beta_Z_std  <- mean(unlist(bZ_ch, use.names = FALSE))

  # PIP pooled over all chains.
  PIP <- mean(g_dr)

  # ---- Step 8e: convergence diagnostics ------------------------------------
  # The Gelman-Rubin statistic for each monitored quantity, and the largest of
  # them as the headline, since one unconverged quantity is enough to make the
  # triplet's summaries unreliable. R-hat requires at least two chains and is
  # NA when N_CHAINS == 1.
  #
  # The joint log-likelihood is monitored alongside the three parameters
  # because it is a single scalar that every parameter in the model feeds
  # into, including the first-stage blocks and the three residual variances,
  # none of which is otherwise diagnosed. It is an addition to the
  # parameter-level statistics and not a substitute for them: the
  # log-likelihood is dominated by the variances, so it can settle while the
  # inclusion indicator is still moving.
  rhat_X  <- gelman_rubin(bX_ch)
  rhat_XZ <- gelman_rubin(bXZ_ch)
  rhat_g  <- gelman_rubin(g_ch)
  rhat_ll <- gelman_rubin(ll_ch)
  rhat_max <- suppressWarnings(
    max(c(rhat_X, rhat_XZ, rhat_g, rhat_ll), na.rm = TRUE))
  if (!is.finite(rhat_max)) rhat_max <- NA_real_

  # ---- Step 8a: Credible intervals -----------------------------------------
  ci_X  <- stats::quantile(bX_dr,  c(0.025, 0.5, 0.975), names = FALSE)
  ci_XZ <- stats::quantile(bXZ_dr, c(0.025, 0.5, 0.975), names = FALSE)

  # ---- Step 8b: Sign-reversal threshold tau --------------------------------
  # The receptor-modulated ligand effect is  effect(z) = bX + bXZ * z, which
  # changes sign at  tau = -bX / bXZ  (in SD units of receptor expression).
  #
  # tau is a RATIO of two parameters, so its posterior can only be obtained
  # from PAIRED draws -- posterior means are not sufficient. Draws for which
  # bXZ is numerically negligible are dropped, since tau diverges there.
  ok_tau  <- is.finite(bXZ_dr) & abs(bXZ_dr) > 1e-8
  tau_dr  <- if (any(ok_tau)) -bX_dr[ok_tau] / bXZ_dr[ok_tau] else NA_real_
  tau_q   <- if (any(ok_tau)) {
    stats::quantile(tau_dr, c(0.025, 0.5, 0.975), names = FALSE)
  } else rep(NA_real_, 3)

  # Probability that the sign reversal falls inside the observed receptor
  # range: this, not tau itself, determines whether a reversal is actually
  # realised in the donor population.
  z_obs_rng   <- range(as.numeric(Z_c) / sd_Z, na.rm = TRUE)
  p_tau_inrng <- if (any(ok_tau)) {
    mean(tau_dr >= z_obs_rng[1] & tau_dr <= z_obs_rng[2])
  } else NA_real_

  # ---- Step 8c: Monte Carlo error of the PIP -------------------------------
  # The PIP is a chain average of a binary indicator. Because the gamma chain
  # is autocorrelated, the naive binomial standard error understates the true
  # Monte Carlo uncertainty; a batch-means estimator is used instead.
  #
  # The estimator is applied WITHIN each chain and the results combined in
  # quadrature, rather than applied to the concatenated draws, so that no
  # batch straddles the join between two chains that began at different
  # starting values.
  mc <- mcse_pooled(lapply(g_ch, mcse_batch))

  # ---- Step 8d: Instrument strength ----------------------------------------
  F_L <- first_stage_F(G_raw, RNA.count.adj_Cell1[ligand_gene_id,   ], V_donor)
  F_R <- first_stage_F(H_raw, RNA.count.adj_Cell2[receptor_gene_id, ], V_donor)

  # Z_vec: receptor expression on the SD scale (for plotting)
  # Effect_vec: total ligand effect as a function of receptor level
  #   effect(z) = beta_X^(s) + beta_XZ^(s) * z,  z = Z_c / sd(Z)
  z_c        <- as.numeric(Z_c)
  Z_vec      <- z_c / sd_Z
  Effect_vec <- Beta_X_std + Beta_XZ_std * Z_vec

  # ---- Step 8f: pointwise 95% credible band for the effect curve -----------
  # effect(z) is evaluated on every draw, then summarised by quantiles at each
  # z. Using PAIRED draws is essential: beta_X and beta_XZ are strongly
  # correlated in the posterior, so a band built by combining their separate
  # marginal intervals would be far too wide.
  #
  # outer() gives a (draws x z) matrix; adding bX recycles down the columns,
  # which is correct here because the matrix has one row per draw.
  .bi <- if (length(bX_dr) > BAND_DRAWS) {
    round(seq(1, length(bX_dr), length.out = BAND_DRAWS))
  } else seq_along(bX_dr)
  .eff_mat <- outer(bXZ_dr[.bi], Z_vec) + bX_dr[.bi]
  .eff_q   <- apply(.eff_mat, 2, stats::quantile,
                    probs = c(0.025, 0.975), names = FALSE)
  Effect_lo <- .eff_q[1, ]
  Effect_hi <- .eff_q[2, ]
  rm(.eff_mat)

  # ---- Step 9: Store output row ----
  out_idx <- out_idx + 1
  Output_list[[out_idx]] <- tibble::tibble(
    ligand_col_name   = ligand_gene_name,
    receptor_col_name = receptor_gene_name,
    pathway_name      = pathway_name,
    # Adaptive pathway representation
    representation    = best_name,
    representation_cor     = best_cor,
    representation_cor_gap = best_cor_gap,
    Beta_X_mean       = Beta_X_std,
    Beta_XZ_mean      = Beta_XZ_std,
    Beta_Z_mean       = Beta_Z_std,
    gamma_mean        = PIP,             # PIP = P(gamma=1 | data), pooled
    # Credible intervals (standardized scale)
    Beta_X_lo         = ci_X[1],  Beta_X_med  = ci_X[2],  Beta_X_hi  = ci_X[3],
    Beta_XZ_lo        = ci_XZ[1], Beta_XZ_med = ci_XZ[2], Beta_XZ_hi = ci_XZ[3],
    # Sign-reversal threshold
    tau_lo            = tau_q[1], tau_med = tau_q[2], tau_hi = tau_q[3],
    tau_in_range_prob = p_tau_inrng,
    z_min             = z_obs_rng[1],
    z_max             = z_obs_rng[2],
    # MCMC reliability
    PIP_mcse          = unname(mc["mcse"]),
    PIP_ess           = unname(mc["ess"]),
    n_chains          = N_CHAINS,
    n_iter            = N_ITER,            # per chain; records the protocol
    ld_r2_max         = LD_R2_MAX,         # instrument clumping threshold
    n_keep            = length(g_dr),      # pooled across chains
    rhat_Beta_X       = rhat_X,            # NA unless N_CHAINS >= 2
    rhat_Beta_XZ      = rhat_XZ,
    rhat_gamma        = rhat_g,
    rhat_loglik       = rhat_ll,
    rhat_max          = rhat_max,          # max over the four above
    # Instrument strength
    n_iv_L            = ncol(G),
    n_iv_R            = ncol(H),
    F_L               = F_L,
    F_R               = F_R,
    Z_vec             = list(Z_vec),
    Effect_vec        = list(Effect_vec),
    Effect_lo         = list(Effect_lo),
    Effect_hi         = list(Effect_hi)
  )
}

# Combine all output rows into a single data frame
Output_MR_CCC <- dplyr::bind_rows(Output_list[seq_len(out_idx)])

cat("\nDone.", out_idx, "triplets analysed.\n")

############################################################
## SECTION 6a: ERROR-RATE CONTROL AND DIAGNOSTIC SUMMARY
##
## Two discovery rules are computed for every triplet:
##
##   PIP > pip_thresh   the median probability model
##                      (Barbieri & Berger, 2004), the standard default
##                      in Bayesian variable selection;
##
##   Bayesian FDR       the largest set whose posterior expected false
##                      discovery proportion does not exceed bfdr_alpha.
##
## Both are reported so that the sensitivity of the findings to the choice
## of rule is visible rather than implicit. DISCOVERY_RULE selects which one
## drives the figures.
############################################################
Output_MR_CCC <- Output_MR_CCC %>%
  mutate(
    bFDR         = bayes_fdr(gamma_mean),
    disc_pip     = gamma_mean > pip_thresh,
    disc_bfdr    = bfdr_threshold_set(gamma_mean, bfdr_alpha)
  )

## A truncated or shortened run is flagged again here, beside the numbers it
## invalidates, so that console output copied out of context cannot be mistaken
## for a reportable result.
if (isTRUE(.mrccc_truncated) || N_ITER != 100000L) {
  cat("\n!!!! THIS RUN IS NOT REPORTABLE:",
      if (isTRUE(.mrccc_truncated)) "truncated triplet universe" else "",
      if (isTRUE(.mrccc_truncated) && N_ITER != 100000L) "and" else "",
      if (N_ITER != 100000L) paste0("N_ITER = ", N_ITER) else "", "\n")
}

cat("\n--- Discovery rules ---\n")
cat("  PIP >", pip_thresh, ":", sum(Output_MR_CCC$disc_pip), "triplets\n")
if (sum(Output_MR_CCC$disc_pip) > 0) {
  cat("      Bayesian FDR of that set:",
      round(mean(1 - Output_MR_CCC$gamma_mean[Output_MR_CCC$disc_pip]), 3), "\n")
}
cat("  Bayesian FDR <=", bfdr_alpha, ":", sum(Output_MR_CCC$disc_bfdr),
    "triplets\n")

# Nested view: how many discoveries at a range of FDR levels
for (a in c(0.01, 0.05, 0.10, 0.20)) {
  cat("      FDR <=", format(a, nsmall = 2), ":",
      sum(bfdr_threshold_set(Output_MR_CCC$gamma_mean, a)), "\n")
}

cat("\n--- MCMC reliability ---\n")
cat("  chains per triplet        :", N_CHAINS,
    if (N_CHAINS > 1) "(dispersed starts)\n" else "\n")
cat("  retained draws per triplet:", unique(Output_MR_CCC$n_keep)[1], "\n")

n_trip <- nrow(Output_MR_CCC)

# The Gelman-Rubin statistic is the reported convergence diagnostic. It
# compares the variance between chains with the variance within them, so it is
# defined only for N_CHAINS >= 2 and only informative when the chains start
# from dispersed values, which INIT_SCALE > 0 and the alternating init_gamma
# above ensure. R-hat is reported for beta_X, beta_XZ, gamma and the joint
# log-likelihood, with the largest of the four as the headline.
#
# 1.01 is used as the threshold rather than the older 1.1, following the
# recommendation that 1.1 is too permissive for the precision that modern
# chain lengths afford.
if (N_CHAINS > 1 && any(is.finite(Output_MR_CCC$rhat_max))) {
  cat("  Gelman-Rubin R-hat : max over all triplets and parameters",
      round(max(Output_MR_CCC$rhat_max, na.rm = TRUE), 4), "\n")
  n_bad_rhat <- sum(Output_MR_CCC$rhat_max > 1.01, na.rm = TRUE)
  cat("  triplets with R-hat > 1.01:", n_bad_rhat, "of", n_trip,
      if (n_bad_rhat == 0) "(all converged)\n" else
        "(RERUN THESE WITH A LONGER CHAIN)\n")

  # The INCLUSION INDICATOR deserves its own line, because it is the only
  # parameter whose convergence bears on a reported decision: discovery is
  # declared by comparing its chain average with a threshold. A flag on
  # beta_X or beta_XZ affects an interval that is reported with its full
  # width; a flag on gamma would affect a claim. These are reported
  # separately so the two are not conflated.
  # max() of an empty vector returns -Inf with a warning, so the
  # discovery line is guarded: a pair with no discoveries is an ordinary
  # outcome across the twenty ordered pairs, not an error.
  rhat_g_disc <- Output_MR_CCC$rhat_gamma[Output_MR_CCC$disc_pip]
  rhat_g_disc <- rhat_g_disc[is.finite(rhat_g_disc)]
  cat("  gamma (inclusion indicator) only:\n")
  cat("      max R-hat over ALL triplets  :",
      round(max(Output_MR_CCC$rhat_gamma, na.rm = TRUE), 4),
      sprintf("(%d of %d above 1.01)\n",
              sum(Output_MR_CCC$rhat_gamma > 1.01, na.rm = TRUE), n_trip))
  if (length(rhat_g_disc) > 0L) {
    cat("      max R-hat over DISCOVERIES   :",
        round(max(rhat_g_disc), 4),
        sprintf("(%d of %d above 1.01)\n",
                sum(rhat_g_disc > 1.01), length(rhat_g_disc)))
  } else {
    cat("      max R-hat over DISCOVERIES   : none declared\n")
  }
} else {
  cat("  Gelman-Rubin R-hat : not available (requires N_CHAINS >= 2)\n")
}

# Mixing quality underlies every Monte Carlo summary reported above: the
# effective sample size, and therefore the MCSE, is the number of retained
# draws divided by the integrated autocorrelation time. Reporting the implied
# IACT makes that link visible rather than leaving the precision unexplained.
iact <- Output_MR_CCC$n_keep / Output_MR_CCC$PIP_ess
cat("  implied integrated autocorrelation time for gamma: median",
    round(median(iact, na.rm = TRUE)), " max", round(max(iact, na.rm = TRUE)),
    "\n")

cat("  PIP Monte Carlo SE : median",
    round(median(Output_MR_CCC$PIP_mcse, na.rm = TRUE), 4),
    " max", round(max(Output_MR_CCC$PIP_mcse, na.rm = TRUE), 4), "\n")
cat("  PIP effective size : median",
    round(median(Output_MR_CCC$PIP_ess, na.rm = TRUE), 1),
    " min", round(min(Output_MR_CCC$PIP_ess, na.rm = TRUE), 1), "\n")

# A triplet whose PIP sits within ~2 MCSE of the threshold cannot be
# reliably classified: its discovery status is within Monte Carlo noise.
borderline <- with(Output_MR_CCC,
                   abs(gamma_mean - pip_thresh) < 2 * PIP_mcse)
cat("  triplets within 2 MCSE of the PIP threshold:",
    sum(borderline, na.rm = TRUE),
    "(their discovery status is not resolved by this chain length)\n")

# Which of the six representations the adaptive rule actually selected.
# The selection frequencies characterise how the adaptive rule behaves and
# are reported alongside the results; a small median margin over the
# runner-up is evidence that the choice is not doing much work.
cat("\n--- Pathway representation (adaptive rule) ---\n")
print(table(Output_MR_CCC$representation))
cat("  median margin over runner-up (|cor| gap):",
    round(median(Output_MR_CCC$representation_cor_gap, na.rm = TRUE), 4), "\n")

cat("\n--- Instrument strength ---\n")
cat("  F (ligand)  : median", round(median(Output_MR_CCC$F_L, na.rm = TRUE), 2),
    " range", paste(round(range(Output_MR_CCC$F_L, na.rm = TRUE), 2),
                    collapse = " - "), "\n")
cat("  F (receptor): median", round(median(Output_MR_CCC$F_R, na.rm = TRUE), 2),
    " range", paste(round(range(Output_MR_CCC$F_R, na.rm = TRUE), 2),
                    collapse = " - "), "\n")
cat("  triplets with both F > 10:",
    sum(Output_MR_CCC$F_L > 10 & Output_MR_CCC$F_R > 10, na.rm = TRUE),
    "of", nrow(Output_MR_CCC), "\n")

############################################################
## SECTION 6b: SAVE AND RELOAD RESULTS
##
## Serialise Output_MR_CCC to an RDS file so that plots
## can be regenerated without re-running the Gibbs sampler.
## The filename encodes the sender-receiver pair.
############################################################
dir.create("Results", showWarnings = FALSE)

## Output names. A run made under either smoke-test override writes to
## names prefixed "smoke_", so that a check can never overwrite the result
## file of a reportable run, and so that regenerate_figures.R and run_pair.R,
## which look for the canonical names, never pick up a check by mistake.
## A run at any clumping threshold other than the canonical 0.8 is a check
## too, for the same reason, and is prefixed in the same way.
.out_prefix <- if (isTRUE(.mrccc_truncated) || N_ITER != 100000L ||
                   !isTRUE(all.equal(LD_R2_MAX, 0.8))) "smoke_" else ""
.rds_path   <- paste0("Results/", .out_prefix, Cell1, "_", Cell2, "_MR_CCC.rds")

saveRDS(Output_MR_CCC, file = .rds_path)
Output_MR_CCC <- readRDS(.rds_path)
cat("Results written:", .rds_path, "\n")

############################################################
## SECTION 7: PLOTS
##
## Three publication-quality figures:
##   p1 -- Bubble plot of |beta_X| and |beta_XZ| by pathway
##   p2 -- Ranked lollipop of all triplets by PIP
##   p3 -- Receptor-modulated ligand effect curves
##          (only PIP > pip_thresh pairs displayed)
##
## Each figure is built by a function in figure_functions.R (sourced
## through figure_style.R), which is the single definition of that figure;
## regenerate_figures.R calls the same functions on the saved RDS. Theme,
## palettes and notation come from figure_style.R.
############################################################

# ---- Prepare shared plotting variables ----
# pair_label is the HTML-formatted "sender -> receiver" title fragment.
# add_plot_columns() adds the ligand-receptor label and the `discover` flag,
# which follows DISCOVERY_RULE; both underlying rules remain available as
# disc_pip and disc_bfdr.
pair_label    <- pair_label_for(Cell1, Cell2)
Output_MR_CCC <- add_plot_columns(Output_MR_CCC, DISCOVERY_RULE)

n_disc  <- sum(Output_MR_CCC$discover)
n_total <- nrow(Output_MR_CCC)

# Plot 1 -- Bubble plot of |beta_X| and |beta_XZ|
p1 <- plot_bubble(Output_MR_CCC, pair_label, pip_thresh)
print(p1)

# Plot 2 -- Ranked lollipop of all triplets by PIP
p2 <- plot_pip_ranking(Output_MR_CCC, pair_label, pip_thresh)
print(p2)

# Plot 3 -- Receptor-modulated ligand effect curves (discoveries only)
p3 <- plot_effect_curves(Output_MR_CCC, pair_label, pip_thresh)
print(p3)

############################################################
## SECTION 8: SAVE PLOTS
##
## Filenames are auto-built from Cell1 and Cell2 so that
## running the script for different pairs produces separate
## files without overwriting earlier results. save_fig()
## writes vector PDF through cairo_pdf into Plots/.
############################################################
save_fig(p1 + labs(title = NULL, subtitle = NULL),
         paste0(.out_prefix, "Supp_", Cell1, "_", Cell2, "_bubble"),      12, 16)
save_fig(p2, paste0(.out_prefix, "Supp_", Cell1, "_", Cell2, "_pip_ranking"), 14, 13)
save_fig(p3, paste0(.out_prefix, "Supp_", Cell1, "_", Cell2, "_curves"),      14, 9)

cat("Plots saved to Plots/ directory.\n")
