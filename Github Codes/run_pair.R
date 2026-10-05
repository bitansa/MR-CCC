############################################################
## run_pair.R
##
## Run the MR-CCC primary analysis for ONE ordered cell-type
## pair, in a FRESH R process.
##
## Usage (from the project root, not from Github Codes/):
##
##   Rscript "Github Codes/run_pair.R" NKCells MonocytesCells
##
## and to sweep every ordered pair, see run_all_pairs.sh.
##
## WHY A SEPARATE PROCESS PER PAIR.
## The analysis holds a genes x donors expression matrix for
## each of the five cell types plus a SNPs x donors genotype
## matrix, which together run to several gigabytes. Sourcing
## the analysis repeatedly inside one long-lived session lets
## that memory accumulate across pairs, and a sweep of twenty
## pairs takes roughly a day -- long enough that an
## out-of-memory failure late in the sweep would be expensive.
## A fresh process per pair guarantees the memory is returned
## to the operating system between pairs, at a cost of a
## minute or two of reloading per pair against about an
## hour of computation.
##
## Each pair writes its own
##   Results/<Cell1>_<Cell2>_MR_CCC.rds
## so a sweep can be resumed after an interruption: pairs
## whose output already exists are skipped (see SKIP_EXISTING).
############################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript \"Github Codes/run_pair.R\" <Cell1> <Cell2>")
}

.mrccc_cell1 <- args[1]
.mrccc_cell2 <- args[2]

## Set to FALSE to force a re-run of pairs that already have output.
SKIP_EXISTING <- TRUE

VALID <- c("BCells", "CD4Cells", "CD8Cells", "NKCells", "MonocytesCells")
if (!.mrccc_cell1 %in% VALID || !.mrccc_cell2 %in% VALID) {
  stop("Cell types must be two of: ", paste(VALID, collapse = ", "))
}
if (identical(.mrccc_cell1, .mrccc_cell2)) {
  stop("Sender and receiver must differ.")
}

out_rds <- file.path("Results",
                     paste0(.mrccc_cell1, "_", .mrccc_cell2, "_MR_CCC.rds"))
if (SKIP_EXISTING && file.exists(out_rds)) {
  cat("SKIP", .mrccc_cell1, "->", .mrccc_cell2,
      "-- output already exists at", out_rds, "\n")
  quit(save = "no", status = 0)
}

t0 <- Sys.time()
cat("START", .mrccc_cell1, "->", .mrccc_cell2, "at",
    format(t0, "%Y-%m-%d %H:%M:%S"), "\n")

## Prerequisites, in the same order as an interactive run.
Rcpp::sourceCpp("Github Codes/mr_ccc_gibbs.cpp")
source("Github Codes/build_lr_database.R")
source("Github Codes/real_data_analysis.R")

el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat("DONE ", .mrccc_cell1, "->", .mrccc_cell2,
    sprintf(" in %.1f minutes\n", el))
