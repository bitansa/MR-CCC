############################################################
## instrument_overlap_sensitivity.R
##
## INSTRUMENT OVERLAP DIAGNOSTIC AND SENSITIVITY ANALYSIS
##
## Assesses whether the ligand and receptor instrument sets
## are separated. Because instruments
## are selected independently for each gene from a +/-200 kb
## cis window, a single SNP may act as an eQTL for BOTH the
## ligand and the receptor. Such a SNP violates the exclusion
## restriction: it influences the receiver pathway through two
## routes rather than one.
##
## Note that simple window overlap is NOT the right criterion.
## A SNP can lie in only one gene's window and still be
## strongly associated with the other gene (through linkage
## disequilibrium or longer-range regulation). The diagnostic
## implemented here therefore tests, for every selected
## instrument, its association with the OTHER gene.
##
## The script runs in two tiers:
##
##   TIER 1 (cheap; no MCMC)
##     For every ordered cell-type pair and every triplet,
##     record: whether the two cis windows intersect, whether
##     the selected instrument sets share SNPs, and how many
##     instruments are cross-associated with the other gene.
##     Also records first-stage F-statistics, which are
##     required separately for the instrument-strength assessment.
##
##   TIER 2 (expensive; re-runs the Gibbs sampler)
##     For one focal pair, re-runs MR-CCC after removing
##     cross-associated SNPs from BOTH instrument sets, and
##     compares PIPs and discovery decisions against the
##     primary analysis.
##
## IMPORTANT: cross-associated SNPs are dropped from both sets,
## never reassigned to the more strongly associated gene.
## Reassignment would leave the second pathway unmodelled and
## would make the exclusion-restriction violation worse rather
## than better.
##
## PREREQUISITES (run before sourcing this file):
##   Rcpp::sourceCpp("Github Codes/mr_ccc_gibbs.cpp")
##   source("Github Codes/build_lr_database.R")   # produces lr_matrix, key_cols
##
## DATA:
##   Same two Zenodo files used by real_data_analysis.R:
##     B_T_NK_monocytes.rda, donor.rda
##   https://doi.org/10.5281/zenodo.19675075
##
## OUTPUT (written to Results/):
##   overlap_detail_all_pairs.csv     -- one row per triplet, all pairs
##   overlap_summary_all_pairs.csv    -- one row per ordered pair
##   instrument_count_sensitivity.csv -- one row per k in K_GRID
##   selection_variant_summary.csv    -- one row per selection variant
##   selection_variant_detail.csv     -- one row per triplet and variant
##   leave_one_out_F_profile_<Cell1>_<Cell2>.csv
##                                    -- first-stage F under single drops
##   sensitivity_<Cell1>_<Cell2>.csv  -- primary vs disjoint run
##   sessionInfo_instrument_overlap_sensitivity.txt
## The first three cover all twenty pairs; the rest the focal pair only.
## Each name is prefixed "smoke_" when LD_R2_MAX is not the canonical 0.8.
##
## Run in a FRESH R session: the script reloads all data itself.
############################################################
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)
library(Matrix)

## GenomicRanges is a BIOCONDUCTOR package and is required by both tiers.
## AUCell, UCell and GSVA are also Bioconductor packages but are needed
## ONLY by Tier 2 (pathway activity scores), so they are loaded lazily
## further down. This means the Tier 1 diagnostic can be run with just
## GenomicRanges installed.
##
## Bioconductor packages are NOT on CRAN -- install.packages() will fail
## with "package is not available for this version of R". Use:
##
##   install.packages("BiocManager")
##   BiocManager::install(c("GenomicRanges", "GenomeInfoDb",
##                          "AUCell", "UCell", "GSVA"),
##                        dependencies = TRUE)
##
## GenomeInfoDb is a runtime dependency of GenomicRanges (countOverlaps()
## fails without it) and is easy to miss, because GenomicRanges itself
## will load successfully while the dependency is absent.
##
## If that fails, check that R and Bioconductor versions are compatible:
##   R.version.string ; BiocManager::version() ; BiocManager::valid()
library(GenomicRanges)

# Fail early with an actionable message rather than deep inside the sweep.
for (.pkg in c("GenomeInfoDb", "IRanges", "S4Vectors")) {
  if (!requireNamespace(.pkg, quietly = TRUE)) {
    stop("Required Bioconductor package '", .pkg, "' is not installed.\n",
         "  BiocManager::install(\"", .pkg, "\", dependencies = TRUE)")
  }
}

############################################################
## USER SETTINGS
############################################################

# Root directory containing the two Zenodo .rda files
# (B_T_NK_monocytes.rda and donor.rda), downloadable from
# https://doi.org/10.5281/zenodo.19675075
#
# The default assumes ~/data/mrccc. To run against a different
# local copy WITHOUT editing (and therefore without any risk of committing
# a machine-specific path), set the environment variable instead:
#
#   Sys.setenv(MRCCC_DATA_DIR = "/path/to/your/data")   # this session only
#
# or add a permanent line to ~/.Renviron:
#
#   MRCCC_DATA_DIR=/path/to/your/data
DATA_DIR <- Sys.getenv("MRCCC_DATA_DIR",
                       unset = file.path("~", "data", "mrccc"))
cat("DATA_DIR:", DATA_DIR, "\n")

# Cell types to include in the Tier 1 sweep.
ALL_CELLS <- c("BCells", "CD4Cells", "CD8Cells",
               "NKCells", "MonocytesCells")

# Focal pair (main-text axis). Used by Tier 2 and by the selection-variant
# sweep in Section 4c, and MUST be included in TIER1_PAIRS so that the
# baseline self-check in 4c has something to compare against.
FOCAL_CELL1 <- "NKCells"
FOCAL_CELL2 <- "MonocytesCells"

# Ordered pairs for the Tier 1 sweep. NULL means all 20 (5 x 4).
# This is the ONLY assignment of TIER1_PAIRS in the script; it must not
# be assigned again further down, because a later assignment would
# silently override this one.
TIER1_PAIRS <- NULL

# Run the expensive Tier 2 re-run? Set FALSE to do diagnostics only.
# Tier 1 is cheap (no MCMC) and can be run on its own; Tier 2 re-runs
# the Gibbs sampler and should be enabled only when the diagnostic has
# identified triplets worth re-running. Tier 2 also compares against the
# primary results in Results/<Cell1>_<Cell2>_MR_CCC.rds, so it must run after
# those exist. Assigning .mrccc_run_tier2 <- FALSE before sourcing runs the
# diagnostics alone without editing this file.
RUN_TIER2 <- if (exists(".mrccc_run_tier2")) isTRUE(.mrccc_run_tier2) else TRUE

pip_thresh <- 0.5      # discovery threshold

# A selected instrument is called "cross-associated" when its
# eQTL p-value for the OTHER gene falls below this threshold.
# Both a nominal and a Bonferroni-corrected rule are recorded so
# that conclusions can be shown not to depend on the choice.
cross_p_nominal <- 0.05

# ---- Instrument-count sensitivity ----------------------------------
# The primary analysis keeps up to 10 SNPs in order of marginal cis-eQTL
# p-value, LD-clumped at LD_R2_MAX. K_GRID re-evaluates the selection diagnostics at
# one count below and one above that, so the choice of 10 is bracketed
# from both sides rather than defended from one.
#
# This costs almost nothing. get_snp_selection() returns the candidate
# SNPs SORTED by p-value, so if it is asked for the top max(K_GRID) the
# top-k set for any smaller k is the first k entries of the same vector:
# no re-selection, no re-fitting of the eQTL scan, and no MCMC. Only the
# first-stage F has to be recomputed, which is a pair of linear models.
#
# PRIMARY_K is the count the primary analysis uses. Every "primary"
# column in the Tier 1 diagnostic is computed from the first PRIMARY_K
# entries, so that those columns are unchanged by widening K_GRID.
#
# The trade-off being measured is that more instruments raise explained
# variance but also add weak-instrument dilution and more opportunities
# for cross-association, so F need not increase monotonically in k.
PRIMARY_K <- 10L
K_GRID    <- c(5L, 10L, 20L)

# ---- LD clumping of the selected instruments ------------------------
# Instruments are taken in eQTL p-value order, skipping any SNP whose
# squared dosage correlation with one already kept exceeds LD_R2_MAX.
# The rationale and the choice of 0.8 (|r| = 0.9, the largest correlation
# covered by simulation S5) are set out beside LD_R2_MAX in
# real_data_analysis.R. Clumping is greedy in p-value order, so the first k
# kept are the same for every k and the K_GRID prefixes remain valid.
# LD_R2_MAX = 1 disables clumping.
#
# KEEP IN SYNC with LD_R2_MAX in real_data_analysis.R.
LD_R2_MAX <- if (exists(".mrccc_ld_r2_max")) as.numeric(.mrccc_ld_r2_max) else 0.8
stopifnot(is.numeric(LD_R2_MAX), LD_R2_MAX > 0, LD_R2_MAX <= 1)
stopifnot(PRIMARY_K %in% K_GRID)

# ---- Output names ----------------------------------------------------
# A run at a clumping threshold other than the canonical 0.8 is a check, and
# writes to names prefixed "smoke_", by the rule used in
# real_data_analysis.R, so that it never overwrites a reported result. A
# "smoke_" prefix already set in this session is honoured as well.
.out_prefix <- if ((exists(".out_prefix") && identical(.out_prefix, "smoke_")) ||
                   !isTRUE(all.equal(LD_R2_MAX, 0.8))) "smoke_" else ""

# ---- Selection-variant sensitivity ---------------------------------
# The primary analysis selects instruments from a +/-200 kb window with
# NO significance threshold. Two alternative strategies are evaluated.
#
#   window_bp    half-width of the cis window, in base pairs. The primary
#                value is 200000. A wider window admits more candidate
#                SNPs, which can raise first-stage strength, but also
#                raises the chance that one SNP instruments BOTH genes --
#                so the two effects must be reported together, never
#                separately.
#
#   lead_p_max   if finite, a gene is retained only when its STRONGEST
#                cis-eQTL reaches this p-value. NA (the primary setting)
#                imposes no threshold at all, which is one reason
#                instrument strength is low: many genes have no lead
#                cis-eQTL at p < 0.05. Imposing the filter raises mean strength but shrinks the
#                triplet universe; the question is whether the discovery
#                set survives.
#
# WINDOW_BP sets the window used by the PRIMARY diagnostic in Section 4
# and by Tier 2; it must stay at 200000 to match real_data_analysis.R.
# The alternative settings live in SELECTION_VARIANTS below and are applied
# only within Section 4c, so nothing else in this script is affected.
WINDOW_BP <- 200000

# Variants swept by Section 4c: the cis window at one width below and one
# above the primary, so that +/-200 kb is bracketed from both sides. The
# middle row MUST reproduce the primary analysis, so that every comparison
# is against a reproduced baseline rather than against remembered numbers.
#
# A lead-eQTL p-value filter is supported through the lead_p_max argument
# of get_snp_selection() (default NA, i.e. no filter) but is not swept
# here; sensitivity to instrument strength is examined directly in
# simulation scenario S4.
SELECTION_VARIANTS <- list(
  list(label = "narrow_100kb",  window_bp = 100000, lead_p_max = NA_real_),
  list(label = "primary_200kb", window_bp = 200000, lead_p_max = NA_real_),
  list(label = "wide_500kb",    window_bp = 500000, lead_p_max = NA_real_)
)

# Fixed seed (the value used by the other scripts). Each script draws its own
# stream; no stream is shared between scripts. donor.rda, loaded below,
# contains a saved .Random.seed, which load() restores in the global
# environment, so the draws made after loading (the Tier 2 sampler and AUCell
# rankings) follow that saved state and are reproducible from it.
set.seed(20260921)

############################################################
## LOAD DATA
############################################################
load(file.path(DATA_DIR, "B_T_NK_monocytes.rda"))
load(file.path(DATA_DIR, "donor.rda"))

# Keep pristine copies: the QC filter below is pair-specific and
# must be re-applied from scratch for every ordered pair.
donor_raw    <- donor
genotype_raw <- genotype

## The genotype dosage conversion is performed AFTER Section 0, so
## that input validation fails fast rather than after several minutes
## of pointless work on a very large matrix.

############################################################
## SECTION 0: INPUT VALIDATION
##
## The pipeline assumes that `donor` (rows), `genotype`
## (columns) and the per-cell-type expression matrices
## (columns) describe THE SAME DONORS IN THE SAME ORDER: a
## single logical vector, donor.filter, is used to subset all
## three. If a donor.rda from a different project or cohort
## subset is used, every downstream step will still execute
## and will silently return misaligned -- and therefore
## meaningless -- results.
##
## These checks fail loudly instead.
############################################################
.required <- c("donor", "genotype", "vcf.ref", "gene.ref",
               "B_Cells", "CD4_Cells", "CD8_Cells", "NK_Cells", "monocytes")
.missing  <- .required[!vapply(.required, exists, logical(1))]
if (length(.missing) > 0) {
  stop("Objects missing after loading the .rda files: ",
       paste(.missing, collapse = ", "))
}

# Prerequisites from OTHER scripts, checked here so that a missing one
# fails immediately rather than after the multi-minute genotype conversion.
.prereq <- c("lr_matrix", "key_cols")
if (RUN_TIER2) .prereq <- c(.prereq, "mr_ccc_gibbs")
.missing <- .prereq[!vapply(.prereq, exists, logical(1))]
if (length(.missing) > 0) {
  stop("Missing: ", paste(.missing, collapse = ", "),
       ".\n  Run  source(\"Github Codes/build_lr_database.R\")",
       if (RUN_TIER2) "  and  Rcpp::sourceCpp(\"Github Codes/mr_ccc_gibbs.cpp\")",
       "  first.")
}

# Donor metadata must carry the covariates used as instruments' controls.
.need_cols <- c("age", "sex", "PC1", "PC2", "PC3")
.miss_cols <- setdiff(.need_cols, colnames(donor))
if (length(.miss_cols) > 0) {
  stop("`donor` is missing required column(s): ",
       paste(.miss_cols, collapse = ", "),
       "\n  Present columns: ", paste(colnames(donor), collapse = ", "))
}

if (!"gene_id" %in% colnames(GenomicRanges::mcols(gene.ref))) {
  stop("`gene.ref` has no `gene_id` metadata column; found: ",
       paste(colnames(GenomicRanges::mcols(gene.ref)), collapse = ", "))
}

cat("\n--- Input dimensions ---\n")
cat("  donor:    ", nrow(donor),      "rows x", ncol(donor), "cols\n")
cat("  genotype: ", nrow(genotype),   "SNPs x", ncol(genotype), "donors\n")
cat("  vcf.ref:  ", length(vcf.ref),  "ranges\n")
cat("  gene.ref: ", length(gene.ref), "ranges\n")

# 1. genotype columns must correspond one-to-one with donor rows
if (ncol(genotype) != nrow(donor)) {
  stop("Donor mismatch: genotype has ", ncol(genotype),
       " columns but `donor` has ", nrow(donor), " rows. ",
       "These must be the same donors in the same order.")
}

# 2. vcf.ref must describe every row of the genotype matrix
if (length(vcf.ref) != nrow(genotype)) {
  stop("SNP mismatch: `vcf.ref` has ", length(vcf.ref),
       " ranges but `genotype` has ", nrow(genotype), " rows.")
}

# 3. expression data must describe the same donors, for every cell type.
# Each object is a list with one element per donor, so its length is the
# donor count.
.n_expr <- vapply(list(BCells = B_Cells, CD4Cells = CD4_Cells,
                       CD8Cells = CD8_Cells, NKCells = NK_Cells,
                       MonocytesCells = monocytes),
                  length, integer(1))
if (any(.n_expr != nrow(donor))) {
  stop("Donor mismatch: the expression data have ",
       paste(names(.n_expr), .n_expr, sep = " = ", collapse = ", "),
       " donors but `donor` has ", nrow(donor), " rows.\n",
       "  This usually means donor.rda and B_T_NK_monocytes.rda come ",
       "from different cohort subsets. Results would be silently wrong.")
}

# 4. Identifier check (informational only).
#
# The pipeline aligns donors POSITIONALLY: a single logical vector,
# donor.filter, subsets donor rows, genotype columns and both expression
# matrices together. real_data_analysis.R never consults an identifier --
# it uses only donor$age, donor$sex and donor$PC1-PC3. Positional
# alignment is therefore the contract, and equal donor counts (checks 1-3
# above) are what actually has to hold.
#
# Identifiers are nevertheless compared where they exist, because a
# mismatch is a useful hint that two files came from different cohort
# subsets. In donor.rda the genotype column names are PLINK-style "FID_IID"
# labels, `donor` carries its identifier in an `id` COLUMN (its row names
# are R's default 1..n) in the same "FID_IID" form, and the expression
# matrices have no donor names. The labels can therefore be compared
# directly. When they are not identical, the IID (the text after the first
# underscore, which assumes that the FID contains no underscore) is
# compared position by position, and a mismatch there stops the run.
.gn <- colnames(genotype)

# Treat R's default row names ("1","2",...) as absent rather than as IDs.
.dn <- rownames(donor)
if (!is.null(.dn) && identical(.dn, as.character(seq_len(nrow(donor))))) {
  .dn <- NULL
}
# Prefer an explicit id column when present.
.did <- if ("id" %in% colnames(donor)) as.character(donor$id) else .dn

if (!is.null(.gn) && !is.null(.did)) {
  if (identical(.gn, .did)) {
    cat("  donor identifiers match exactly between genotype and donor\n")
  } else if (setequal(.gn, .did)) {
    warning("genotype columns and donor ids contain the same donors but in ",
            "DIFFERENT ORDER. Because the pipeline aligns positionally, ",
            "reorder before proceeding:\n",
            "  genotype <- genotype[, match(as.character(donor$id), colnames(genotype))]")
  } else if (all(grepl("^[^_]+_.+$", .gn))) {
    # Genotype columns in FID_IID form: compare their IIDs with the donor
    # identifiers, themselves reduced to an IID when they have the same form.
    .iid_of <- function(x) sub("^[^_]+_", "", x)
    .g_iid  <- .iid_of(.gn)
    .d_iid  <- if (all(grepl("^[^_]+_.+$", .did))) .iid_of(.did) else .did
    if (!identical(.g_iid, .d_iid)) {
      .bad <- which(.g_iid != .d_iid)
      stop("Donor mismatch: the IIDs parsed from the genotype column names ",
           "differ from the donor identifiers at ", length(.bad),
           " position(s), first at position ", .bad[1], " (\"", .gn[.bad[1]],
           "\" vs \"", .did[.bad[1]], "\"). The pipeline aligns donors ",
           "positionally, so the files must list the same donors in the ",
           "same order.")
    }
    cat("  genotype IIDs (from FID_IID column names) match the donor ",
        "identifiers position by position\n", sep = "")
  } else {
    cat("  NOTE: genotype column names and donor$id use different labelling\n",
        "        (e.g. \"", utils::head(.gn, 1), "\" vs \"",
        utils::head(.did, 1), "\"). Alignment is positional and donor\n",
        "        counts agree, matching real_data_analysis.R.\n", sep = "")
  }
} else {
  cat("  NOTE: donors are identified positionally (no usable names on\n",
      "        genotype/donor/expression). Counts agree, which is the\n",
      "        same guarantee the primary pipeline relies on.\n", sep = "")
}
cat("--- Input validation passed ---\n\n")

# ---- Convert genotypes to numeric dosage ONCE ----------------------
# The string-to-dosage conversion does not depend on which donors
# survive QC, only the column subset does. Doing it once here rather
# than inside the per-pair function avoids repeating an expensive
# conversion of a very large matrix twenty times over.
cat("Converting genotypes to dosage (once)...\n")
genotype.mat.full <- matrix(NA_real_,
                            nrow = nrow(genotype_raw),
                            ncol = ncol(genotype_raw))
# Missing calls ("./.") stay NA here; they are imputed per pair, after the
# pair's donor QC, by the SNP's mean dosage -- as in real_data_analysis.R.
genotype.mat.full[genotype_raw == "0/0"] <- 0
genotype.mat.full[genotype_raw == "0/1"] <- 1
genotype.mat.full[genotype_raw == "1/1"] <- 2
cat("  dosage matrix:", nrow(genotype.mat.full), "SNPs x",
    ncol(genotype.mat.full), "donors\n")

# Any genotype code other than "./.", "0/0", "0/1", "1/1" (for example a
# phased code such as "0|1") would also leave an NA and be imputed as if it
# were a missing call. Such codes are counted separately so that the
# situation is visible if it ever arises.
.n_bad_geno <- sum(is.na(genotype.mat.full) & genotype_raw != "./.")
if (.n_bad_geno > 0) {
  warning("Genotype dosage contains ", .n_bad_geno,
          " unrecognised genotype codes. Observed codes: ",
          paste(utils::head(unique(as.vector(genotype_raw)), 10),
                collapse = ", "),
          ". They are treated as missing calls.")
} else {
  cat("  all genotype codes recognised;", sum(is.na(genotype.mat.full)),
      "missing calls\n")
}

############################################################
## SECTION 1: DONOR-LEVEL AGGREGATION (all five cell types)
############################################################
RNA.count_BCells   <- sapply(B_Cells,   colSums)
RNA.count_CD4Cells <- sapply(CD4_Cells, colSums)
RNA.count_CD8Cells <- sapply(CD8_Cells, colSums)
RNA.count_NKCells  <- sapply(NK_Cells,  colSums)

get_gene_counts <- function(x) {
  if (!is.null(dim(x)) && length(dim(x)) == 2) return(colSums(x))
  if (is.atomic(x) && is.null(dim(x)))           return(x)
  stop("Unexpected object in monocytes: class = ",
       paste(class(x), collapse = "+"))
}
RNA.count_MonocytesCells <- sapply(monocytes, get_gene_counts)

valid_genes <- rownames(RNA.count_BCells)

## ---- Shared preprocessing helpers ---------------------------------------
## Copies of the definitions in real_data_analysis.R, which this script does
## not source. KEEP IN SYNC with them.
##
## build_covariates(): donor-level covariates V (sex, 10-year age bins,
##   PC1-PC3), used in the eQTL screen, the F statistic and the model.
## centre_cols(): column centring; the model has no first-stage intercept, so
##   G, H and V are centred before they reach the sampler, as X, Z, Y are.
## orient_pc1(): fixes the arbitrary sign of the PC1 pathway score so that
##   higher scores mean higher mean expression of the pathway genes.
build_covariates <- function(donor_df) {
  age <- donor_df$age
  cbind(
    MALE = as.integer(donor_df$sex == "male"),
    AGE  = ifelse(age < 30, 1L, ifelse(age < 40, 2L, ifelse(age < 50, 3L,
           ifelse(age < 60, 4L, ifelse(age < 70, 5L, ifelse(age < 80, 6L, 7L)))))),
    PC1  = donor_df$PC1, PC2 = donor_df$PC2, PC3 = donor_df$PC3
  )
}
centre_cols <- function(M) {
  M <- as.matrix(M)
  sweep(M, 2, colMeans(M), "-")
}
orient_pc1 <- function(pc, ref) {
  r <- suppressWarnings(stats::cor(pc, ref, use = "complete.obs"))
  if (is.finite(r) && r < 0) -pc else pc
}

############################################################
## SECTION 2: PER-PAIR PREPROCESSING
##
## Reproduces exactly the QC, normalisation and genotype
## processing used in real_data_analysis.R, but wrapped in a
## function so that it can be applied to any ordered pair.
## Returns everything the downstream steps need.
############################################################
prepare_pair <- function(Cell1, Cell2) {

  RNA.count_Cell1 <- get(paste0("RNA.count_", Cell1))
  RNA.count_Cell2 <- get(paste0("RNA.count_", Cell2))

  # ---- Library-size QC (identical rule to the primary pipeline) ----
  lib.size1 <- colSums(RNA.count_Cell1) / median(colSums(RNA.count_Cell1))
  lib.size2 <- colSums(RNA.count_Cell2) / median(colSums(RNA.count_Cell2))

  donor.filter1 <- (lib.size1 <= median(lib.size1) + 3 * mad(lib.size1)) &
    (lib.size1 > 0.25)
  donor.filter2 <- (lib.size2 <= median(lib.size2) + 3 * mad(lib.size2)) &
    (lib.size2 > 0.25)
  donor.filter  <- donor.filter1 & donor.filter2

  donor_p <- donor_raw[donor.filter, ]
  C1f <- RNA.count_Cell1[, donor.filter]
  C2f <- RNA.count_Cell2[, donor.filter]

  # ---- Library-size normalisation ----
  ls1 <- colSums(C1f) / median(colSums(C1f))
  ls2 <- colSums(C2f) / median(colSums(C2f))
  # Normalised counts, as in the primary pipeline
  adj1 <- sweep(C1f, 2, ls1, "/")
  adj2 <- sweep(C2f, 2, ls2, "/")

  # ---- Subset the pre-computed dosage matrix, then carrier filter ----
  # The filter must be recomputed here because carrier rates depend on
  # which donors survived this pair's QC. It is evaluated over called
  # genotypes; missing calls are then imputed by the SNP's mean dosage.
  gmat <- genotype.mat.full[, donor.filter, drop = FALSE]
  keep <- rowMeans(gmat > 0, na.rm = TRUE) > 0.05
  keep[is.na(keep)] <- FALSE
  gmat <- gmat[keep, , drop = FALSE]
  vcf  <- vcf.ref[keep]
  na_g <- which(is.na(gmat), arr.ind = TRUE)
  if (nrow(na_g) > 0) gmat[na_g] <- rowMeans(gmat, na.rm = TRUE)[na_g[, 1]]

  V_p <- build_covariates(donor_p)
  if (anyNA(V_p)) stop("Missing donor covariates in pair ", Cell1, " -> ", Cell2)

  # SNP identifiers of the rows of gmat, from the row names of the character
  # genotype matrix where present (NULL otherwise). Used as labels only.
  snp_id <- if (is.null(rownames(genotype_raw))) NULL else
    rownames(genotype_raw)[keep]

  list(donor = donor_p, V = V_p, gmat = gmat, vcf = vcf, snp_id = snp_id,
       adj1 = adj1, adj2 = adj2, n = ncol(adj1))
}

# Promoter windows used for cis instrument selection.
#
# The flank is recorded as an attribute on the returned object. This is not
# decoration: get_snp_selection() caches its results, and the cache key is
# built from that attribute, so a selection computed under one window can
# never be silently reused under another.
promoter_ref_for <- function(bp) {
  pr <- promoters(gene.ref, upstream = bp, downstream = bp)
  attr(pr, "flank_bp") <- bp
  pr
}

# Primary window (+/-200 kb), matching real_data_analysis.R.
gene.promoter.ref <- promoter_ref_for(WINDOW_BP)

############################################################
## Pathway gene lookup
##
## Reproduces exactly the rule used in real_data_analysis.R.
## CellPhoneDB stores a classification string such as
## "Signaling by Interleukin"; the trailing keyword is extracted
## and matched (case-insensitively) against MSigDB Reactome
## gene-set names. Note this is a REGEX MATCH on gs_name, not an
## equality test -- one CellPhoneDB class maps to many Reactome
## sets.
##
## A triplet whose pathway yields no genes present in the
## receiver matrix is skipped by the primary pipeline, so the
## same filter is applied here to keep the triplet universe
## identical to the primary analysis.
############################################################
pathway_gene_ids_for <- function(pathway_name, RNA_mat_receiver, key_tbl) {
  keyword <- sub(".*by\\s+", "", pathway_name)
  # An empty keyword would match every gene set; the primary skips such rows.
  if (!nzchar(keyword)) return(character(0))
  hits <- key_tbl %>%
    dplyr::filter(str_detect(gs_name, regex(keyword, ignore_case = TRUE)))
  ids <- unique(hits$ensembl_gene)
  ids[ids %in% rownames(RNA_mat_receiver)]
}

############################################################
## SECTION 3: INSTRUMENT SELECTION HELPERS
##
## get_snp_selection() mirrors get_SNP_matrix() in
## real_data_analysis.R, but additionally returns the SELECTED
## SNP INDICES and their eQTL p-values. The indices are what
## make it possible to compare the ligand and receptor
## instrument sets against one another.
##
## Selection is cached per (ordered pair, gene): many triplets
## within a pair share the same ligand or receptor, and the
## marginal eQTL scan is by far the most expensive part of the
## diagnostic.
##
## The cache key MUST include the ordered pair, not just the
## cell type. The donor QC filter depends on BOTH cell types,
## so the surviving donor set -- and therefore the eQTL scan --
## differs between, say, NK -> monocytes and NK -> B cells.
## Caching on cell type alone would silently reuse a selection
## computed on the wrong donor set.
############################################################
.sel_cache <- new.env(parent = emptyenv())

# Greedy LD clumping in p-value order; see the fuller comment on the copy in
# real_data_analysis.R. `cand` holds row indices into `geno`, strongest
# first; returns the kept indices in the same order.
#
# KEEP IN SYNC with clump_by_ld() in real_data_analysis.R.
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

get_snp_selection <- function(gene_id, RNA_mat_adj, P,
                              cache_tag, min_snps = 1L, keep_max = 10L,
                              promoter_ref = gene.promoter.ref,
                              lead_p_max   = NA_real_) {

  # The cache key carries the ordered pair (via cache_tag), the window
  # half-width and the lead-p threshold. The last two are read from the
  # arguments rather than supplied by the caller, so a selection variant
  # cannot accidentally inherit a cached result from a different variant
  # because someone forgot to change a tag.
  fb <- attr(promoter_ref, "flank_bp")
  if (is.null(fb)) {
    stop("promoter_ref lacks the flank_bp attribute; build it with ",
         "promoter_ref_for(). Without it two windows would share one cache ",
         "entry and the variant comparison would silently report identical ",
         "selections.")
  }
  # format(scientific = FALSE) so that 200000 and 200000L key identically
  # (as.character(2e5) is "2e+05", as.character(200000L) is "200000").
  #
  # keep_max is part of the key too. Section 4 asks for the top 20 so that
  # the k-sensitivity can subset; Tier 2 asks for the top 10. Without
  # keep_max in the key, whichever ran first would be served to the other.
  # The LD threshold is in the key for the same reason: a selection clumped
  # at one threshold must never be served to a request made at another.
  ck <- paste0(cache_tag, "::", gene_id,
               "::w", format(fb, scientific = FALSE, trim = TRUE),
               "::p", ifelse(is.na(lead_p_max), "none", lead_p_max),
               "::k", keep_max,
               "::r2", format(LD_R2_MAX, scientific = FALSE, trim = TRUE))
  # exists()/get()/assign() rather than [[ ]], so that a cached
  # "no valid selection" result (stored as NULL) is remembered
  # rather than silently recomputed on every triplet.
  if (exists(ck, envir = .sel_cache, inherits = FALSE)) {
    return(get(ck, envir = .sel_cache, inherits = FALSE))
  }

  gene.ind <- which(gene.ref$gene_id == gene_id)
  if (length(gene.ind) == 0) {
    assign(ck, NULL, envir = .sel_cache); return(NULL)
  }

  snp.ind <- which(countOverlaps(P$vcf, promoter_ref[gene.ind]) > 0)
  if (length(snp.ind) < min_snps) {
    assign(ck, NULL, envir = .sel_cache); return(NULL)
  }

  # Marginal cis-eQTL p-value for each candidate SNP
  snp.pval <- rep(NA_real_, length(snp.ind))
  for (j in seq_along(snp.ind)) {
    g_j <- P$gmat[snp.ind[j], ]
    # A SNP with no variation among these donors has no eQTL p-value; lm()
    # would drop its term and row 2 would belong to a covariate.
    if (var(g_j) == 0) next
    lm.j <- lm(RNA_mat_adj[gene_id, ] ~ g_j + P$V)
    cf <- summary(lm.j)$coefficients
    snp.pval[j] <- if ("g_j" %in% rownames(cf)) cf["g_j", 4] else NA_real_
  }

  ok <- which(!is.na(snp.pval))
  if (length(ok) < min_snps) {
    assign(ck, NULL, envir = .sel_cache); return(NULL)
  }

  # Optional significance filter on the LEAD SNP. Applied to the gene, not
  # to individual instruments: if the strongest cis-eQTL for this gene does
  # not reach lead_p_max, the gene contributes no instruments at all and
  # every triplet involving it drops out. This is an alternative selection
  # strategy; the primary analysis passes lead_p_max = NA and therefore
  # skips this block entirely.
  if (is.finite(lead_p_max) && min(snp.pval[ok]) > lead_p_max) {
    assign(ck, NULL, envir = .sel_cache); return(NULL)
  }

  # Up to keep_max SNPs in eQTL p-value order, LD-clumped at LD_R2_MAX.
  # Identical to get_SNP_matrix() in real_data_analysis.R, which is what
  # the primary analysis fits.
  ranked    <- snp.ind[ok[order(snp.pval[ok])]]
  best_snps <- clump_by_ld(ranked, P$gmat, keep_max = keep_max,
                           r2_max = LD_R2_MAX)
  pos       <- match(best_snps, snp.ind)   # back to positions in snp.pval

  out <- list(
    idx        = best_snps,                 # global row indices into P$gmat
    pval       = snp.pval[pos],             # eQTL p-values, own gene
    window_idx = snp.ind,                   # all candidates in the window
    gene_ind   = gene.ind
  )
  assign(ck, out, envir = .sel_cache)
  out
}

# Association p-value of a set of SNPs with an ARBITRARY gene.
# Used to test a ligand instrument against the receptor, and
# vice versa: this is the cross-association screen.
cross_assoc_pvals <- function(snp_idx, gene_id, RNA_mat_adj, P) {
  vapply(snp_idx, function(s) {
    g_s <- P$gmat[s, ]
    if (var(g_s) == 0) return(NA_real_)
    lm.s <- lm(RNA_mat_adj[gene_id, ] ~ g_s + P$V)
    cf <- summary(lm.s)$coefficients
    if ("g_s" %in% rownames(cf)) cf["g_s", 4] else NA_real_
  }, numeric(1))
}

# First-stage partial F-statistic for an instrument set:
# compares the covariate-only model with the model that also
# contains all instruments. Reported for the instrument-strength
# assessment, and used here to separate "bias removed" from
# "power lost" when instruments are dropped.
first_stage_F <- function(snp_idx, gene_id, RNA_mat_adj, P) {
  if (length(snp_idx) == 0) return(NA_real_)
  Gm <- t(P$gmat[snp_idx, , drop = FALSE])
  colnames(Gm) <- paste0("snp", seq_len(ncol(Gm)))

  # Both models must be fitted to the SAME rows. Genotype dosages can
  # contain NA (any code other than 0/0, 0/1, 1/1, ./.), and lm() drops
  # incomplete rows per model -- so fitting m0 and m1 separately on the
  # raw data would leave them with different sample sizes and make
  # anova() fail. Restrict to complete cases up front.
  # Covariates are the model's own V (see build_covariates()).
  Vm <- P$V
  colnames(Vm) <- paste0("cov", seq_len(ncol(Vm)))
  dat <- data.frame(y = RNA_mat_adj[gene_id, ], Vm, Gm, check.names = FALSE)
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  if (nrow(dat) <= ncol(dat) + 1L) return(NA_real_)

  cov_names <- colnames(Vm)
  m0 <- lm(y ~ ., data = dat[, c("y", cov_names)])
  m1 <- lm(y ~ ., data = dat)
  a  <- tryCatch(anova(m0, m1), error = function(e) NULL)
  if (!is.null(a) && nrow(a) >= 2 && "F" %in% colnames(a)) {
    a[2, "F"]
  } else NA_real_
}

## Largest absolute pairwise dosage correlation between two instrument
## index sets. Called with the same set twice it returns the strongest
## within-set correlation, excluding the diagonal; called with the two
## different sets it returns the strongest cross-set correlation, where
## the diagonal carries no special meaning and is kept. Returns NA when
## there is no pair to correlate, or when a SNP is monomorphic in the
## retained donors and its correlation is undefined.
max_abs_cor <- function(idx_a, idx_b, gmat) {
  if (length(idx_a) == 0L || length(idx_b) == 0L) return(NA_real_)
  same <- identical(idx_a, idx_b)
  if (same && length(idx_a) < 2L) return(NA_real_)
  R <- suppressWarnings(stats::cor(t(gmat[idx_a, , drop = FALSE]),
                                   t(gmat[idx_b, , drop = FALSE])))
  if (same) R <- R[upper.tri(R)]
  if (all(is.na(R))) return(NA_real_)
  max(abs(R), na.rm = TRUE)
}

############################################################
## SECTION 4: TIER 1 -- OVERLAP DIAGNOSTIC ACROSS ALL PAIRS
##
## No MCMC is run here. For each ordered pair and each triplet
## the following are recorded:
##   window_overlap   -- do the two +/-200 kb windows intersect?
##   n_shared_sel     -- SNPs present in BOTH selected sets
##   n_cross_L, n_cross_R -- instruments cross-associated with
##                           the other gene (nominal and
##                           Bonferroni rules)
##   F_L, F_R         -- first-stage F-statistics
############################################################
diagnose_pair <- function(Cell1, Cell2,
                          promoter_ref = gene.promoter.ref,
                          lead_p_max   = NA_real_) {

  t0 <- Sys.time()
  cat("\n[Tier 1] ", Cell1, "->", Cell2, " ... ")
  P <- prepare_pair(Cell1, Cell2)

  lr_p <- lr_matrix %>%
    filter(ligand_ensembl   %in% valid_genes,
           receptor_ensembl %in% valid_genes) %>%
    distinct(ligand_symbol, ligand_ensembl,
             receptor_symbol, receptor_ensembl, pathway_name)

  rows <- vector("list", nrow(lr_p)); k <- 0

  key_cols_p <- key_cols %>% filter(ensembl_gene %in% valid_genes)

  for (i in seq_len(nrow(lr_p))) {
    lig <- lr_p$ligand_ensembl[i]
    rec <- lr_p$receptor_ensembl[i]
    if (!lig %in% rownames(P$adj1)) next
    if (!rec %in% rownames(P$adj2)) next

    # Apply the same pathway filter as the primary pipeline so that
    # the triplet universe here matches that of the primary analysis.
    if (length(pathway_gene_ids_for(lr_p$pathway_name[i],
                                    P$adj2, key_cols_p)) == 0) next

    # Ask for the top max(K_GRID) so the k-sensitivity can subset. The
    # PRIMARY analysis uses the first PRIMARY_K of these, and every
    # "primary" column below is computed from that prefix only.
    selL <- get_snp_selection(lig, P$adj1, P,
                              cache_tag    = paste0(Cell1, "_", Cell2, "_L"),
                              keep_max     = max(K_GRID),
                              promoter_ref = promoter_ref,
                              lead_p_max   = lead_p_max)
    selR <- get_snp_selection(rec, P$adj2, P,
                              cache_tag    = paste0(Cell1, "_", Cell2, "_R"),
                              keep_max     = max(K_GRID),
                              promoter_ref = promoter_ref,
                              lead_p_max   = lead_p_max)
    if (is.null(selL) || is.null(selR)) next

    # Primary instrument sets: the first PRIMARY_K by p-value.
    primL <- selL$idx[seq_len(min(PRIMARY_K, length(selL$idx)))]
    primR <- selR$idx[seq_len(min(PRIMARY_K, length(selR$idx)))]

    # Two DIFFERENT questions, recorded in two separate columns:
    #
    #   window_overlap      do the two GENOMIC windows intersect? Computed
    #                       from the GRanges directly. This is what the
    #                       column name says and is the quantity summarised
    #                       as n_window_overlap in the per-pair table.
    #
    #   n_shared_candidates how many post-filter candidate SNPs lie in BOTH
    #                       windows. Can be zero even when the windows
    #                       overlap, if the shared region holds no SNP that
    #                       survives the carrier-rate filter.
    #
    # The two must not be conflated. The distinction matters because they
    # answer different questions: "do the windows overlap" has the genomic
    # answer, and "could one SNP instrument both genes" has the candidate
    # answer.
    # Strand is ignored: windows of adjacent genes on opposite strands
    # overlap genomically even though their ranges carry different strands.
    win_ov <- as.logical(countOverlaps(promoter_ref[selL$gene_ind],
                                       promoter_ref[selR$gene_ind],
                                       ignore.strand = TRUE) > 0)
    n_shared_cand <- length(intersect(selL$window_idx, selR$window_idx))

    # SNPs appearing in both PRIMARY instrument sets
    shared_sel <- intersect(primL, primR)

    # Cross-association screen (the criterion that actually matters),
    # computed once on the full top-max(K_GRID) sets. Because selL$idx is
    # ordered by p-value and cross_assoc_pvals() preserves order, the
    # p-values for any prefix are the same prefix of these vectors.
    pL_on_R_all <- cross_assoc_pvals(selL$idx, rec, P$adj2, P)  # ligand IVs vs receptor
    pR_on_L_all <- cross_assoc_pvals(selR$idx, lig, P$adj1, P)  # receptor IVs vs ligand

    # Primary-set slices of the above.
    pL_on_R <- pL_on_R_all[seq_along(primL)]
    pR_on_L <- pR_on_L_all[seq_along(primR)]

    bonfL <- cross_p_nominal / max(1L, length(primL))
    bonfR <- cross_p_nominal / max(1L, length(primR))

    # ---- Instrument-count sensitivity ------------------------------
    # The top-kk instruments are the first kk entries of the p-value-
    # sorted selection, and the matching cross-association p-values are
    # the same prefix of *_all, so nothing is recomputed except F.
    #
    # The loop variable is kk, NOT k: k is the row counter for this
    # pair and must not be overwritten.
    kv <- list()
    for (kk in K_GRID) {
      iL <- selL$idx[seq_len(min(kk, length(selL$idx)))]
      iR <- selR$idx[seq_len(min(kk, length(selR$idx)))]
      kv[[paste0("F_L_k",      kk)]] <- first_stage_F(iL, lig, P$adj1, P)
      kv[[paste0("F_R_k",      kk)]] <- first_stage_F(iR, rec, P$adj2, P)
      kv[[paste0("n_iv_L_k",   kk)]] <- length(iL)
      kv[[paste0("n_iv_R_k",   kk)]] <- length(iR)
      kv[[paste0("n_shared_k", kk)]] <- length(intersect(iL, iR))
      kv[[paste0("n_cross_k",  kk)]] <-
        sum(pL_on_R_all[seq_along(iL)] < cross_p_nominal, na.rm = TRUE) +
        sum(pR_on_L_all[seq_along(iR)] < cross_p_nominal, na.rm = TRUE)
    }

    k <- k + 1
    rows[[k]] <- dplyr::bind_cols(
      tibble(
        Cell1 = Cell1, Cell2 = Cell2,
        ligand = lr_p$ligand_symbol[i], receptor = lr_p$receptor_symbol[i],
        pathway = lr_p$pathway_name[i],
        # All "primary" columns use the first PRIMARY_K instruments, which
        # is exactly what real_data_analysis.R fits.
        n_iv_L = length(primL), n_iv_R = length(primR),
        window_overlap      = win_ov,
        n_shared_candidates = n_shared_cand,
        n_shared_selected   = length(shared_sel),
        n_cross_L_nom  = sum(pL_on_R < cross_p_nominal, na.rm = TRUE),
        n_cross_R_nom  = sum(pR_on_L < cross_p_nominal, na.rm = TRUE),
        n_cross_L_bonf = sum(pL_on_R < bonfL, na.rm = TRUE),
        n_cross_R_bonf = sum(pR_on_L < bonfR, na.rm = TRUE),
        F_L = first_stage_F(primL, lig, P$adj1, P),
        F_R = first_stage_F(primR, rec, P$adj2, P),
        # Pairwise dosage correlation among the selected instruments.
        # Reported because clumping bounds but does not remove the
        # correlation, so the degree actually present has to be stated
        # rather than assumed: scenario S5 sweeps AR(1) correlation up to 0.9, and
        # these columns say where the real instrument sets sit on that
        # axis. max_abs_r_LR is the cross-set quantity and is the one that
        # bears on identification -- the two exposures can only be
        # separated to the extent that their instruments are not shared
        # information, and unlike window_overlap it is a continuous
        # measure rather than a yes/no flag.
        max_abs_r_L  = max_abs_cor(primL, primL, P$gmat),
        max_abs_r_R  = max_abs_cor(primR, primR, P$gmat),
        max_abs_r_LR = max_abs_cor(primL, primR, P$gmat)
      ),
      tibble::as_tibble(kv)
    )
  }

  cat(k, "triplets, ",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
  bind_rows(rows[seq_len(k)])
}

# Sweep every ordered pair (sender != receiver): 5 x 4 = 20 pairs,
# unless TIER1_PAIRS restricts the sweep to a subset.
if (is.null(TIER1_PAIRS)) {
  pairs_grid <- expand.grid(Cell1 = ALL_CELLS, Cell2 = ALL_CELLS,
                            stringsAsFactors = FALSE) %>%
    filter(Cell1 != Cell2)
} else {
  pairs_grid <- data.frame(
    Cell1 = vapply(TIER1_PAIRS, `[`, character(1), 1),
    Cell2 = vapply(TIER1_PAIRS, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )
}
cat("\nTier 1 sweep over", nrow(pairs_grid), "ordered pair(s).\n")

t_start <- Sys.time()
detail_all <- purrr::map2_dfr(pairs_grid$Cell1, pairs_grid$Cell2,
                              diagnose_pair)
cat("\nTier 1 elapsed:",
    round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1),
    "minutes\n")

dir.create("Results", showWarnings = FALSE)

# Per-pair summary
overlap_summary <- detail_all %>%
  group_by(Cell1, Cell2) %>%
  summarise(
    n_triplets            = n(),
    n_window_overlap      = sum(window_overlap),
    n_shared_cand_gt0     = sum(n_shared_candidates > 0),
    n_shared_selected_gt0 = sum(n_shared_selected > 0),
    n_cross_any_nom       = sum((n_cross_L_nom + n_cross_R_nom) > 0),
    n_cross_any_bonf      = sum((n_cross_L_bonf + n_cross_R_bonf) > 0),
    median_F_L            = median(F_L, na.rm = TRUE),
    median_F_R            = median(F_R, na.rm = TRUE),
    # Instrument correlation: the median over triplets of the strongest
    # correlation within each set, and the maximum over triplets of the
    # strongest CROSS-set correlation. The cross-set column is summarised
    # by its maximum rather than its median because the claim it supports
    # is a worst-case one -- that no triplet has instrument sets sharing
    # substantial information.
    median_max_r_L        = median(max_abs_r_L,  na.rm = TRUE),
    median_max_r_R        = median(max_abs_r_R,  na.rm = TRUE),
    max_max_r_LR          = max(max_abs_r_LR,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_cross_nom  = round(100 * n_cross_any_nom  / n_triplets, 1),
    pct_cross_bonf = round(100 * n_cross_any_bonf / n_triplets, 1)
  )

write.csv(detail_all,
          paste0("Results/", .out_prefix, "overlap_detail_all_pairs.csv"),
          row.names = FALSE)
write.csv(overlap_summary,
          paste0("Results/", .out_prefix, "overlap_summary_all_pairs.csv"),
          row.names = FALSE)

cat("\n================ TIER 1 SUMMARY ================\n")
print(as.data.frame(overlap_summary))
cat("\nWritten to Results/", .out_prefix, "overlap_summary_all_pairs.csv\n",
    sep = "")

############################################################
## SECTION 4b: INSTRUMENT-COUNT SENSITIVITY (k in K_GRID)
##
## Collapses the per-k columns produced above into one row per
## instrument count, so that the effect of k on instrument
## strength and on cross-association can be read directly.
##
## This is the evidence behind the choice of 10 instruments.
## Two quantities move in opposite directions as k grows:
## explained variance in the first stage (which helps) and the
## number of instruments that are also associated with the
## other gene (which hurts). The table shows both.
############################################################
k_sensitivity <- purrr::map_dfr(K_GRID, function(kk) {
  fl <- detail_all[[paste0("F_L_k",      kk)]]
  fr <- detail_all[[paste0("F_R_k",      kk)]]
  sh <- detail_all[[paste0("n_shared_k", kk)]]
  cr <- detail_all[[paste0("n_cross_k",  kk)]]
  tibble(
    k               = kk,
    n_triplets      = nrow(detail_all),
    median_F_L      = round(median(fl, na.rm = TRUE), 3),
    median_F_R      = round(median(fr, na.rm = TRUE), 3),
    max_F_L         = round(max(fl, na.rm = TRUE), 3),
    max_F_R         = round(max(fr, na.rm = TRUE), 3),
    n_both_F_gt10   = sum(fl > 10 & fr > 10, na.rm = TRUE),
    n_shared_gt0    = sum(sh > 0, na.rm = TRUE),
    n_cross_gt0     = sum(cr > 0, na.rm = TRUE),
    pct_cross       = round(100 * sum(cr > 0, na.rm = TRUE) / nrow(detail_all), 1)
  )
})

write.csv(k_sensitivity,
          paste0("Results/", .out_prefix, "instrument_count_sensitivity.csv"),
          row.names = FALSE)

cat("\n========= INSTRUMENT-COUNT SENSITIVITY =========\n")
print(as.data.frame(k_sensitivity))
cat("\nWritten to Results/", .out_prefix, "instrument_count_sensitivity.csv\n",
    sep = "")

############################################################
## SECTION 4c: SELECTION-VARIANT SENSITIVITY
##
## The alternative selection settings are run together
## because they are the same operation with different
## parameters: re-select the instruments under an alternative
## rule and re-measure both instrument STRENGTH and instrument
## OVERLAP.
##
## Reporting the two together is the point. A wider window
## buys strength by admitting more candidate SNPs, but the
## same widening makes it more likely that one SNP sits in
## both genes' windows -- which is precisely the exclusion-
## restriction problem diagnosed in Section 4. A table showing
## only the improvement in F would be misleading.
##
## No MCMC runs here. This section establishes WHICH variants
## change the instrument sets enough to be worth re-running;
## Tier 2 does the re-running.
##
## The second variant, "primary_200kb", reproduces the
## primary analysis, so the baseline in the comparison is
## recomputed rather than quoted. Its numbers must match
## Section 4's.
############################################################
cat("\nSelection-variant sweep:", length(SELECTION_VARIANTS),
    "variants on", FOCAL_CELL1, "->", FOCAL_CELL2, "\n")

# Clear the selection cache so the baseline variant is a genuine
# re-derivation. Without this, "primary_200kb" would build the same cache
# key as Section 4 and simply read the stored selection back, and the
# self-check below could never detect a regression in the selection step
# itself -- it would be checking a copy against itself. Cost: one extra
# eQTL scan for the focal pair. Tier 2 re-hits the cache afterwards.
rm(list = ls(envir = .sel_cache), envir = .sel_cache)

variant_detail <- list()
variant_summary <- purrr::map_dfr(SELECTION_VARIANTS, function(v) {
  cat("  [", v$label, "] window +/-", v$window_bp / 1000, "kb, lead p ",
      ifelse(is.na(v$lead_p_max), "unrestricted",
             paste0("< ", v$lead_p_max)), "\n", sep = "")

  d <- diagnose_pair(FOCAL_CELL1, FOCAL_CELL2,
                     promoter_ref = promoter_ref_for(v$window_bp),
                     lead_p_max   = v$lead_p_max)

  variant_detail[[v$label]] <<- dplyr::mutate(d, variant = v$label)

  if (nrow(d) == 0) {
    return(tibble(variant = v$label, window_kb = v$window_bp / 1000,
                  lead_p = v$lead_p_max, n_triplets = 0L))
  }
  tibble(
    variant        = v$label,
    window_kb      = v$window_bp / 1000,
    lead_p         = v$lead_p_max,
    n_triplets     = nrow(d),
    median_F_L     = round(median(d$F_L, na.rm = TRUE), 3),
    median_F_R     = round(median(d$F_R, na.rm = TRUE), 3),
    max_F          = round(max(c(d$F_L, d$F_R), na.rm = TRUE), 3),
    n_both_F_gt10  = sum(d$F_L > 10 & d$F_R > 10, na.rm = TRUE),
    # The cost side of the trade-off, reported alongside the benefit.
    n_window_overlap = sum(d$window_overlap, na.rm = TRUE),
    n_shared_gt0     = sum(d$n_shared_selected > 0, na.rm = TRUE),
    n_cross_gt0      = sum((d$n_cross_L_nom + d$n_cross_R_nom) > 0, na.rm = TRUE),
    pct_cross        = round(100 * sum((d$n_cross_L_nom + d$n_cross_R_nom) > 0,
                                       na.rm = TRUE) / nrow(d), 1)
  )
})

variant_detail_all <- dplyr::bind_rows(variant_detail)

write.csv(variant_summary,
          paste0("Results/", .out_prefix, "selection_variant_summary.csv"),
          row.names = FALSE)
write.csv(variant_detail_all,
          paste0("Results/", .out_prefix, "selection_variant_detail.csv"),
          row.names = FALSE)

cat("\n========= SELECTION-VARIANT SENSITIVITY =========\n")
print(as.data.frame(variant_summary))

# Which triplets survive every variant? This is the quantity that
# determines whether the findings depend on the selection rule: if the
# discovery set is a subset of the triplets common to all variants, they
# do not.
if (length(variant_detail) > 1L) {
  key_of <- function(d) paste(d$ligand, d$receptor, d$pathway, sep = "|")
  common <- Reduce(intersect, lapply(variant_detail, key_of))
  cat("\nTriplets present under ALL", length(variant_detail), "variants:",
      length(common), "\n")
  cat("Baseline triplet count (primary rule)   :",
      nrow(variant_detail[[1]]), "\n")
}

# The baseline variant must reproduce Section 4 exactly. If it does not,
# something in the variant machinery has altered the primary analysis,
# which would invalidate every comparison below it.
.base <- variant_detail[["primary_200kb"]]
.ref  <- detail_all[detail_all$Cell1 == FOCAL_CELL1 &
                      detail_all$Cell2 == FOCAL_CELL2, ]
if (!is.null(.base) && nrow(.ref) > 0) {
  if (nrow(.base) != nrow(.ref) ||
      !isTRUE(all.equal(sort(.base$F_L), sort(.ref$F_L)))) {
    warning("Baseline variant does NOT reproduce the Section 4 diagnostic ",
            "for the focal pair (", nrow(.base), " vs ", nrow(.ref),
            " triplets). Do not trust the variant comparison until this ",
            "is resolved.")
  } else {
    cat("Baseline variant reproduces Section 4 exactly.\n")
  }
} else {
  cat("Baseline self-check SKIPPED: focal pair not present in the Tier 1 ",
      "sweep. Add it to TIER1_PAIRS.\n", sep = "")
}

cat("\nWritten to Results/", .out_prefix, "selection_variant_summary.csv\n",
    sep = "")

############################################################
## SECTION 4e: LEAVE-ONE-INSTRUMENT-OUT FIRST-STAGE F PROFILE
##
## leave_one_instrument_out.R reports how each discovery's PIP moves when a
## single instrument is dropped. A discovery that falls below the threshold
## when one instrument is removed admits two very different readings, and the
## PIP alone cannot separate them:
##
##   (a) the first-stage signal is CONCENTRATED in that instrument, so
##       removing it leaves the exposure barely instrumented and the
##       posterior reverts toward the prior. A strength problem.
##
##   (b) the instrument was contributing something other than exposure
##       variation -- a direct path to the pathway, for instance -- so
##       removing it changes the answer while the first stage is untouched.
##       A validity problem: the exclusion restriction, not instrument
##       strength.
##
## The discriminating quantity is the first-stage partial F recomputed on
## the reduced instrument set. If F collapses alongside the PIP, reading (a)
## is demonstrated rather than asserted; if F is unchanged while the PIP
## moves, reading (b) is the live possibility and must be reported as such.
##
## No MCMC is involved: this is two linear models per dropped instrument, so
## the whole profile costs seconds. It is computed for every triplet of the
## focal pair rather than only the discoveries, so that the discoveries can
## be read against the distribution of the rest.
############################################################
loo_F_profile <- function(Cell1 = FOCAL_CELL1, Cell2 = FOCAL_CELL2) {

  cat("\n[4e] Leave-one-out first-stage F profile:", Cell1, "->", Cell2, "\n")
  P <- prepare_pair(Cell1, Cell2)

  key_cols_l <- key_cols %>% filter(ensembl_gene %in% valid_genes)
  lr_l <- lr_matrix %>%
    filter(ligand_ensembl   %in% valid_genes,
           receptor_ensembl %in% valid_genes) %>%
    distinct(ligand_symbol, ligand_ensembl,
             receptor_symbol, receptor_ensembl, pathway_name)

  rows <- list(); k <- 0

  for (i in seq_len(nrow(lr_l))) {
    lig <- lr_l$ligand_ensembl[i]
    rec <- lr_l$receptor_ensembl[i]
    if (!lig %in% rownames(P$adj1)) next
    if (!rec %in% rownames(P$adj2)) next
    if (length(pathway_gene_ids_for(lr_l$pathway_name[i],
                                    P$adj2, key_cols_l)) == 0) next

    selL <- get_snp_selection(lig, P$adj1, P,
                              cache_tag = paste0(Cell1, "_", Cell2, "_L"))
    selR <- get_snp_selection(rec, P$adj2, P,
                              cache_tag = paste0(Cell1, "_", Cell2, "_R"))
    if (is.null(selL) || is.null(selR)) next

    primL <- selL$idx[seq_len(min(PRIMARY_K, length(selL$idx)))]
    primR <- selR$idx[seq_len(min(PRIMARY_K, length(selR$idx)))]

    F_L_full <- first_stage_F(primL, lig, P$adj1, P)
    F_R_full <- first_stage_F(primR, rec, P$adj2, P)

    # One row per dropped instrument. `side` names the block the instrument
    # belongs to; F_same is that block's F after the drop and is the column
    # that matters, while F_other is recorded unchanged so that a reader can
    # confirm the drop touched only one first stage. `dropped` is the SNP
    # identifier where the genotype matrix carries one, and otherwise the
    # side prefix and the instrument's rank (L1, L2, ...; R1, R2, ...).
    snp_lab <- function(idx, s, prefix) {
      id <- if (is.null(P$snp_id)) NA_character_ else P$snp_id[idx[s]]
      if (is.na(id) || !nzchar(id)) paste0(prefix, s) else id
    }
    for (s in seq_along(primL)) {
      k <- k + 1
      rows[[k]] <- tibble(
        ligand = lr_l$ligand_symbol[i], receptor = lr_l$receptor_symbol[i],
        pathway = lr_l$pathway_name[i], side = "ligand",
        dropped = snp_lab(primL, s, "L"),
        F_same_full = F_L_full,
        F_same_drop = first_stage_F(primL[-s], lig, P$adj1, P),
        F_other     = F_R_full,
        n_iv_same   = length(primL))
    }
    for (s in seq_along(primR)) {
      k <- k + 1
      rows[[k]] <- tibble(
        ligand = lr_l$ligand_symbol[i], receptor = lr_l$receptor_symbol[i],
        pathway = lr_l$pathway_name[i], side = "receptor",
        dropped = snp_lab(primR, s, "R"),
        F_same_full = F_R_full,
        F_same_drop = first_stage_F(primR[-s], rec, P$adj2, P),
        F_other     = F_L_full,
        n_iv_same   = length(primR))
    }
  }

  out <- bind_rows(rows) %>%
    mutate(F_ratio = F_same_drop / F_same_full)

  write.csv(out, paste0("Results/", .out_prefix, "leave_one_out_F_profile_",
                        Cell1, "_", Cell2, ".csv"), row.names = FALSE)
  cat("   ", nrow(out), "single-instrument drops profiled;",
      "written to Results/", .out_prefix, "leave_one_out_F_profile_",
      Cell1, "_", Cell2,
      ".csv\n", sep = "")

  # The headline number: how far the first stage falls in the worst case for
  # each triplet. Joined against the PIP profile by hand when writing up,
  # since that file is produced by a different script.
  print(as.data.frame(
    out %>% group_by(ligand, receptor) %>%
      summarise(min_F_ratio = round(min(F_ratio, na.rm = TRUE), 3),
                worst_drop  = dropped[which.min(F_ratio)],
                .groups = "drop") %>%
      arrange(min_F_ratio) %>% head(12)
  ))
  invisible(out)
}

loo_F <- loo_F_profile()

############################################################
## SECTION 5: TIER 2 -- SENSITIVITY RE-RUN (focal pair)
##
## Repeats the primary analysis for the focal pair with all
## cross-associated SNPs removed from BOTH instrument sets,
## then compares PIPs and discovery decisions against the
## primary results stored in Results/<Cell1>_<Cell2>_MR_CCC.rds.
##
## WHY THIS IS A STRESS TEST AND NOT A CORRECTION
##
## The removal criterion here is deliberately more aggressive than the
## instrument conditions require, so the result is a lower bound on
## robustness rather than a better estimate. Two facts set the reading.
##
##   1. Whether a discovery has SNPs in both cis-windows is recorded in
##      Section 4 (n_shared_candidates). Where it is zero, every
##      cross-association removed below is LONG-RANGE -- a ligand-window
##      SNP associated with the receptor gene on another part of the
##      genome.
##
##   2. A long-range association of that kind is expected under chance at
##      this sample size: both instrument sets are screened against the
##      partner gene at a nominal 0.05, and across the screen the share of
##      triplets with at least one hit rises with the number of instruments
##      k (column pct_cross of Results/instrument_count_sensitivity.csv). It is
##      ALSO the signature of a real ligand-to-receptor effect: the
##      outcome model conditions on receptor expression, so a ligand
##      instrument that acts on the receptor THROUGH the ligand is
##      mediation along the path under study, not a violated instrument
##      condition. The three instrument conditions -- association with
##      the exposure, independence of confounders, and no effect on the
##      outcome except through the exposures -- are untouched by it.
##
## Removing these SNPs can therefore delete causal signal as readily as
## confounding, which is why the primary analysis retains them and this
## section is reported as a stress test. A discovery that survives it has
## been shown not to depend on cross-associated instruments at all; one
## that does not survive is fragile rather than refuted, and the
## leave-one-out profile is the diagnostic that distinguishes the two.
##
## Everything other than G and H (X, Z, Y, V, priors, sampler
## settings) is identical to real_data_analysis.R, so any
## difference is attributable to the instrument change alone.
############################################################
if (RUN_TIER2) {

  # Pathway-scoring packages are needed only here. Loading them at this
  # point rather than at the top of the file lets Tier 1 run on a machine
  # where only GenomicRanges is installed.
  for (pkg in c("AUCell", "UCell", "GSVA")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Tier 2 requires the Bioconductor package '", pkg,
           "'. Install with:\n",
           '  install.packages("BiocManager")\n',
           '  BiocManager::install(c("AUCell", "UCell", "GSVA"))')
    }
  }
  library(AUCell); library(UCell); library(GSVA)

  # GSVA changed its interface at GSVA 1.50 (Bioconductor 3.18): the old
  # gsva(expr=, gset.idx.list=, method=) signature was replaced by parameter
  # objects. This wrapper supports both so the script runs on either
  # installation. KEEP IN SYNC with the copy in real_data_analysis.R.
  if (!exists("gsva_scores")) {
    gsva_scores <- function(expr_full, geneSets,
                            method = c("gsva", "ssgsea"),
                            kcdf = "Gaussian") {
      method <- match.arg(method)
      if (exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
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
  }

  Cell1 <- FOCAL_CELL1
  Cell2 <- FOCAL_CELL2
  cat("\n[Tier 2] Sensitivity re-run:", Cell1, "->", Cell2, "\n")

  P <- prepare_pair(Cell1, Cell2)

  key_cols_f <- key_cols %>% filter(ensembl_gene %in% valid_genes)
  lr_f <- lr_matrix %>%
    filter(ligand_ensembl   %in% valid_genes,
           receptor_ensembl %in% valid_genes)

  out_list <- vector("list", nrow(lr_f)); oi <- 0

  for (i in seq_len(nrow(lr_f))) {

    lig <- lr_f$ligand_ensembl[i]
    rec <- lr_f$receptor_ensembl[i]
    pth <- lr_f$pathway_name[i]

    if (!lig %in% rownames(P$adj1)) next
    if (!rec %in% rownames(P$adj2)) next

    pw_genes <- pathway_gene_ids_for(pth, P$adj2, key_cols_f)
    if (length(pw_genes) == 0) next

    # ---- X and Z ----
    X <- matrix(P$adj1[lig, ], ncol = 1)
    Z <- matrix(P$adj2[rec, ], ncol = 1)

    # ---- Y: six representations, best absolute correlation with X ----
    # (identical rule to the primary analysis)
    expr      <- P$adj2[pw_genes, , drop = FALSE]
    expr_full <- P$adj2

    # Oriented against the mean score, as in the primary pipeline.
    Y_PC1 <- orient_pc1(prcomp(t(expr), center = TRUE, scale. = TRUE)$x[, 1],
                        colMeans(expr))

    geneRanks <- AUCell_buildRankings(expr_full, plotStats = FALSE)
    geneSets  <- list(pathway = intersect(rownames(expr_full), pw_genes))
    auc       <- AUCell_calcAUC(geneSets, geneRanks,
                                aucMaxRank = ceiling(0.05 * nrow(expr_full)))
    Y_AUCell  <- as.numeric(getAUC(auc)["pathway", ])

    Y_rowMean <- colMeans(expr)
    Y_UCell   <- as.numeric(ScoreSignatures_UCell(
      expr_full, features = geneSets)[, "pathway_UCell"])
    Y_ssGSEA  <- as.numeric(
      gsva_scores(expr_full, geneSets, method = "ssgsea")["pathway", ])
    Y_GSVA    <- as.numeric(
      gsva_scores(expr_full, geneSets, method = "gsva",
                  kcdf = "Gaussian")["pathway", ])

    Y_list <- list(PC1 = Y_PC1, AUCell = Y_AUCell, Mean = Y_rowMean,
                   UCell = Y_UCell, ssGSEA = Y_ssGSEA, GSVA = Y_GSVA)
    cors   <- sapply(Y_list, function(y) cor(y, X, use = "complete.obs"))
    Y      <- matrix(Y_list[[names(which.max(abs(cors)))]], ncol = 1)
    # Recorded so that it can be compared with the primary analysis below.
    rep_t2 <- names(which.max(abs(cors)))

    sd_X <- sd(as.numeric(X)); sd_Z <- sd(as.numeric(Z)); sd_Y <- sd(as.numeric(Y))
    X_c <- matrix(as.numeric(X) - mean(X), ncol = 1)
    Z_c <- matrix(as.numeric(Z) - mean(Z), ncol = 1)
    Y_c <- matrix(as.numeric(Y) - mean(Y), ncol = 1)

    # ---- Instruments, with cross-associated SNPs REMOVED ----
    selL <- get_snp_selection(lig, P$adj1, P,
                              cache_tag = paste0(Cell1, "_", Cell2, "_L"))
    selR <- get_snp_selection(rec, P$adj2, P,
                              cache_tag = paste0(Cell1, "_", Cell2, "_R"))
    if (is.null(selL) || is.null(selR)) next

    pL_on_R <- cross_assoc_pvals(selL$idx, rec, P$adj2, P)
    pR_on_L <- cross_assoc_pvals(selR$idx, lig, P$adj1, P)

    keepL <- selL$idx[!(pL_on_R < cross_p_nominal) | is.na(pL_on_R)]
    keepR <- selR$idx[!(pR_on_L < cross_p_nominal) | is.na(pR_on_L)]

    # Also drop any SNP shared between the two selected sets.
    shared <- intersect(keepL, keepR)
    keepL  <- setdiff(keepL, shared)
    keepR  <- setdiff(keepR, shared)

    n_dropL <- length(selL$idx) - length(keepL)
    n_dropR <- length(selR$idx) - length(keepR)

    # A triplet with no surviving instrument on either side cannot
    # be estimated; it is recorded rather than silently skipped.
    if (length(keepL) == 0 || length(keepR) == 0) {
      oi <- oi + 1
      out_list[[oi]] <- tibble(
        ligand = lr_f$ligand_symbol[i], receptor = lr_f$receptor_symbol[i],
        pathway = pth, n_iv_L_new = length(keepL), n_iv_R_new = length(keepR),
        n_dropped_L = n_dropL, n_dropped_R = n_dropR,
        F_L_new = NA_real_, F_R_new = NA_real_,
        Beta_X_new = NA_real_, Beta_XZ_new = NA_real_, PIP_new = NA_real_,
        status = "no_instruments_left",
        representation = rep_t2
      )
      next
    }

    # Centred instruments and covariates, as in the primary pipeline.
    G <- centre_cols(t(P$gmat[keepL, , drop = FALSE]))
    H <- centre_cols(t(P$gmat[keepR, , drop = FALSE]))
    V <- centre_cols(P$V)

    # ---- Gibbs sampler: identical settings to the primary run ----
    #
    # CANONICAL MCMC SETTINGS -- must match real_data_analysis.R exactly:
    # FOUR chains of 100000 iterations, 2000 burn-in, no thinning, from
    # overdispersed starting values.
    #
    # The sampler configuration is part of that contract, not a detail. This
    # run is a sensitivity analysis whose ONLY intended difference from the
    # primary analysis is the instrument set. Any difference in chain length
    # or start would let a PIP gap between the two be Monte Carlo noise
    # rather than the effect of dropping cross-associated SNPs -- which is
    # exactly the quantity this analysis is designed to measure.
    #
    # KEEP IN SYNC with real_data_analysis.R (N_ITER, BURN_IN, THIN,
    # N_CHAINS, INIT_SCALE). Defined locally because this script does not
    # source it (only build_lr_database.R is required).
    N_ITER_T2     <- 100000L
    BURN_IN_T2    <- 2000L
    N_CHAINS_T2   <- 4L
    INIT_SCALE_T2 <- 1.0

    n       <- nrow(X_c)
    g_scale <- min(n, 100.0)
    fits <- lapply(seq_len(N_CHAINS_T2), function(cc) {
      mr_ccc_gibbs(
        X_c, Z_c, Y_c, G, H, V,
        n_iter = N_ITER_T2, burn_in = BURN_IN_T2, thin = 1,
        a_sigma = 3.0, b_sigma = 2.0, a_rho = 3.0, b_rho = 1.0,
        nu1 = 1e-4, gG = g_scale, gH = g_scale, gV = g_scale,
        gZ = g_scale, gBeta = g_scale, ridge = 1e-8,
        init_gamma = if (cc %% 2L == 0L) 0L else 1L,
        init_scale = INIT_SCALE_T2
      )
    })
    bX_dr  <- unlist(lapply(fits, function(f) as.numeric(f$Beta_X_draws)),
                     use.names = FALSE)
    bXZ_dr <- unlist(lapply(fits, function(f) as.numeric(f$Beta_XZ_draws)),
                     use.names = FALSE)
    g_dr   <- unlist(lapply(fits, function(f) as.numeric(f$gamma_draws)),
                     use.names = FALSE)

    oi <- oi + 1
    out_list[[oi]] <- tibble(
      ligand = lr_f$ligand_symbol[i], receptor = lr_f$receptor_symbol[i],
      pathway = pth,
      n_iv_L_new = length(keepL), n_iv_R_new = length(keepR),
      n_dropped_L = n_dropL, n_dropped_R = n_dropR,
      F_L_new = first_stage_F(keepL, lig, P$adj1, P),
      F_R_new = first_stage_F(keepR, rec, P$adj2, P),
      Beta_X_new  = mean(bX_dr)  * sd_X / sd_Y,
      Beta_XZ_new = mean(bXZ_dr) * sd_X * sd_Z / sd_Y,
      PIP_new     = mean(g_dr),
      n_chains    = N_CHAINS_T2,
      status      = "ok",
      representation = rep_t2
    )
  }

  sens <- bind_rows(out_list[seq_len(oi)])

  # ---- Compare against the primary analysis ----
  primary_file <- paste0("Results/", Cell1, "_", Cell2, "_MR_CCC.rds")
  if (file.exists(primary_file)) {
    primary_raw <- readRDS(primary_file)
    primary <- primary_raw %>%
      transmute(ligand = ligand_col_name, receptor = receptor_col_name,
                pathway = pathway_name,
                Beta_X_old = Beta_X_mean, Beta_XZ_old = Beta_XZ_mean,
                PIP_old = gamma_mean)

    # The join is keyed on gene symbols and pathway. A key that occurs more
    # than once on either side would duplicate rows, so it is reported
    # instead of joined.
    .key_dups <- function(d) {
      key <- c("ligand", "receptor", "pathway")
      if (!all(key %in% names(d))) return(0L)
      sum(duplicated(as.data.frame(d)[, key, drop = FALSE]))
    }
    if (.key_dups(primary) > 0L || .key_dups(sens) > 0L) {
      stop("The (ligand, receptor, pathway) key is not unique (",
           .key_dups(primary), " repeated key(s) in the primary results, ",
           .key_dups(sens), " in the Tier 2 results); the comparison would ",
           "duplicate rows.", call. = FALSE)
    }

    # The adaptive representation is re-evaluated here. A choice that
    # differs from the primary analysis means a different outcome variable,
    # and is reported; the comparison proceeds.
    if ("representation" %in% names(primary_raw) &&
        "representation" %in% names(sens)) {
      .rep_cmp <- primary_raw %>%
        transmute(ligand = ligand_col_name, receptor = receptor_col_name,
                  pathway = pathway_name,
                  rep_old = as.character(representation)) %>%
        inner_join(sens %>% dplyr::select(ligand, receptor, pathway,
                                          rep_new = representation),
                   by = c("ligand", "receptor", "pathway")) %>%
        filter(rep_old != rep_new)
      if (nrow(.rep_cmp) > 0L) {
        warning(nrow(.rep_cmp), " triplet(s) use a different pathway ",
                "representation in Tier 2 than in the primary analysis, e.g. ",
                .rep_cmp$ligand[1], "-", .rep_cmp$receptor[1], ": '",
                .rep_cmp$rep_new[1], "' vs '", .rep_cmp$rep_old[1], "'.",
                call. = FALSE)
      }
    }

    sens <- primary %>%
      left_join(sens, by = c("ligand", "receptor", "pathway")) %>%
      mutate(
        disc_old    = PIP_old > pip_thresh,
        disc_new    = PIP_new > pip_thresh,
        PIP_change  = PIP_new - PIP_old,
        decision_changed = disc_old != disc_new
      ) %>%
      arrange(desc(PIP_old))

    cat("\n============== TIER 2 SUMMARY ==============\n")
    cat("Discoveries (primary):        ", sum(sens$disc_old, na.rm = TRUE), "\n")
    cat("Discoveries (disjoint IVs):   ", sum(sens$disc_new, na.rm = TRUE), "\n")
    cat("Decisions changed:            ",
        sum(sens$decision_changed, na.rm = TRUE), "\n")
    cat("Triplets with >=1 SNP dropped:",
        sum((sens$n_dropped_L + sens$n_dropped_R) > 0, na.rm = TRUE), "\n")
    print(as.data.frame(
      sens %>% filter(disc_old) %>%
        dplyr::select(ligand, receptor, PIP_old, PIP_new,
                      n_dropped_L, n_dropped_R,
                      F_L_new, F_R_new, decision_changed)
    ))
  } else {
    warning("Primary results not found at ", primary_file,
            " -- reporting sensitivity run only.")
  }

  write.csv(sens,
            paste0("Results/", .out_prefix, "sensitivity_", Cell1, "_", Cell2,
                   ".csv"),
            row.names = FALSE)
  cat("\nWritten to Results/", .out_prefix, "sensitivity_", Cell1, "_", Cell2,
      ".csv\n", sep = "")
}

############################################################
## INTERPRETATION NOTE
##
## When comparing the two runs, read the PIP change alongside
## the F-statistics. Removing instruments necessarily reduces
## first-stage strength, so an attenuated effect may reflect
## lost power rather than removed bias. A discovery that
## survives with a comparable F is strong evidence that the
## primary finding was not driven by shared instruments; a
## discovery that weakens while F collapses is inconclusive
## rather than refuted.
############################################################

## Package versions used for this run.
writeLines(utils::capture.output(sessionInfo()),
           paste0("Results/", .out_prefix,
                  "sessionInfo_instrument_overlap_sensitivity.txt"))
