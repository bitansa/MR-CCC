############################################################
## build_triplet_inputs.R
##
## Reassemble the model inputs (X, Z, Y, G, H, V) for ONE
## ligand-receptor-pathway triplet, exactly as the main loop
## of real_data_analysis.R does.
##
## WHY THIS IS A SEPARATE FILE.
## Several scripts need to refit an individual triplet after
## the primary analysis has run -- convergence diagnostics,
## posterior sign probabilities, leave-one-instrument-out.
## If each carried its own copy of this logic, the copies
## could diverge: a copy offering fewer than the six pathway
## representations would, for any triplet where AUCell or
## UCell wins the correlation contest, silently fit a
## different outcome variable than the primary analysis. A
## single definition removes that class of inconsistency;
## the representation actually chosen in a refit is also
## checked against the primary output (see below).
##
## KEEP IN SYNC with Steps 1-6 of real_data_analysis.R. If
## the primary pipeline changes how a triplet is assembled,
## this file must change with it.
##
## PREREQUISITES (all inherited from real_data_analysis.R):
##   lr_filtered, key_cols, RNA.count.adj_Cell1,
##   RNA.count.adj_Cell2, V_donor, get_SNP_matrix, gsva_scores,
##   centre_cols, orient_pc1; optionally Output_MR_CCC, against
##   which the re-evaluated pathway representation is checked.
############################################################

## ---- Comparison with the primary representation ----------------------------
## The adaptive rule is re-evaluated whenever a triplet is rebuilt. AUCell
## breaks ties in its gene rankings at random, so the choice can in principle
## differ from the one recorded in the primary output, in which case the
## refit uses a different outcome variable. This helper reports such a
## difference with a warning; it does not stop the refit. The primary output
## is looked up in Output_MR_CCC when that object and its `representation`
## column exist, and the check is skipped otherwise.
check_primary_representation <- function(lig_sym, rec_sym, pth, chosen) {
  if (!exists("Output_MR_CCC") ||
      !"representation" %in% names(Output_MR_CCC)) return(invisible(NA))
  hit <- which(Output_MR_CCC$ligand_col_name   == lig_sym &
               Output_MR_CCC$receptor_col_name == rec_sym &
               Output_MR_CCC$pathway_name      == pth)
  if (length(hit) != 1L) return(invisible(NA))
  primary <- Output_MR_CCC$representation[hit]
  same <- identical(as.character(primary), as.character(chosen))
  if (!same) {
    warning("Pathway representation for ", lig_sym, "-", rec_sym, " (", pth,
            ") is '", chosen, "' in this refit but '", primary,
            "' in the primary analysis.", call. = FALSE)
  }
  invisible(same)
}

build_triplet_inputs <- function(lig_sym, rec_sym, pth, verbose = TRUE) {

  # Triplets are identified by gene SYMBOL. A symbol triple that matches more
  # than one row of lr_filtered (for example one symbol mapped to two Ensembl
  # identifiers) cannot be resolved here and is reported rather than resolved
  # by taking the first match.
  runs <- which(lr_filtered$ligand_symbol   == lig_sym &
                lr_filtered$receptor_symbol == rec_sym &
                lr_filtered$pathway_name    == pth)
  if (length(runs) > 1L) {
    stop("The triplet ", lig_sym, "-", rec_sym, " (", pth, ") matches ",
         length(runs), " rows of lr_filtered (Ensembl IDs: ",
         paste(lr_filtered$ligand_ensembl[runs], lr_filtered$receptor_ensembl[runs],
               sep = "/", collapse = ", "),
         "); it cannot be identified by gene symbol alone.", call. = FALSE)
  }
  run <- runs[1]
  if (is.na(run)) return(NULL)

  lig <- lr_filtered$ligand_ensembl[run]
  rec <- lr_filtered$receptor_ensembl[run]

  # ---- Step 1: pathway gene set --------------------------------------------
  # An empty keyword would make the regex match EVERY gene set and build Y
  # from the whole transcriptome, so such rows are skipped here exactly as
  # the primary pipeline skips them.
  keyword <- sub(".*by\\s+", "", pth)
  if (!nzchar(keyword)) return(NULL)
  pw <- unique(key_cols$ensembl_gene[
    stringr::str_detect(key_cols$gs_name,
                        stringr::regex(keyword, ignore_case = TRUE))])
  pw <- pw[pw %in% rownames(RNA.count.adj_Cell2)]
  if (length(pw) == 0) return(NULL)

  # ---- Step 2: X and Z ------------------------------------------------------
  if (!lig %in% rownames(RNA.count.adj_Cell1)) return(NULL)
  if (!rec %in% rownames(RNA.count.adj_Cell2)) return(NULL)
  X <- matrix(RNA.count.adj_Cell1[lig, ], ncol = 1)
  Z <- matrix(RNA.count.adj_Cell2[rec, ], ncol = 1)

  # ---- Step 3: all SIX pathway representations ------------------------------
  # The adaptive rule selects whichever is most correlated with X. All six
  # must be offered here so that the refit uses the same outcome variable as
  # the primary analysis.
  expr       <- RNA.count.adj_Cell2[pw, , drop = FALSE]
  expr_full  <- RNA.count.adj_Cell2
  path_genes <- intersect(rownames(expr_full), pw)
  geneSets   <- list(pathway = path_genes)

  geneRanks <- AUCell::AUCell_buildRankings(expr_full, plotStats = FALSE)
  auc       <- AUCell::AUCell_calcAUC(
                 geneSets, geneRanks,
                 aucMaxRank = ceiling(0.05 * nrow(expr_full)))

  ucell_scores <- UCell::ScoreSignatures_UCell(
                    expr_full, features = list(pathway = path_genes))

  # PC1 is oriented against the mean score, as in the primary pipeline, so
  # that its sign (and hence the sign of the effects) is reproducible.
  Y_list <- list(
    PC1    = orient_pc1(prcomp(t(expr), center = TRUE, scale. = TRUE)$x[, 1],
                        colMeans(expr)),
    AUCell = as.numeric(AUCell::getAUC(auc)["pathway", ]),
    Mean   = colMeans(expr),
    UCell  = as.numeric(ucell_scores[, "pathway_UCell"]),
    ssGSEA = as.numeric(gsva_scores(expr_full, geneSets, "ssgsea")["pathway", ]),
    GSVA   = as.numeric(gsva_scores(expr_full, geneSets, "gsva")["pathway", ])
  )
  cors      <- vapply(Y_list, function(y) cor(y, X, use = "complete.obs"),
                      numeric(1))
  best_name <- names(which.max(abs(cors)))
  if (verbose) {
    cat("    [", lig_sym, "-", rec_sym, "] representation: ", best_name,
        "\n", sep = "")
  }
  check_primary_representation(lig_sym, rec_sym, pth, best_name)
  Y <- matrix(Y_list[[best_name]], ncol = 1)

  # ---- Steps 4-6: scales, instruments, covariates ---------------------------
  sd_X <- sd(as.numeric(X)); sd_Z <- sd(as.numeric(Z)); sd_Y <- sd(as.numeric(Y))
  if (!is.finite(sd_Y) || sd_Y == 0) return(NULL)

  # Raw dosages from get_SNP_matrix(); the model receives centred columns,
  # as do the covariates (V_donor, built once in real_data_analysis.R).
  G <- get_SNP_matrix(lig, RNA.count.adj_Cell1, min_snps = 1)
  H <- get_SNP_matrix(rec, RNA.count.adj_Cell2, min_snps = 1)
  if (is.null(G) || is.null(H)) return(NULL)
  G <- centre_cols(G)
  H <- centre_cols(H)
  V <- centre_cols(V_donor)

  list(X = matrix(as.numeric(X) - mean(X), ncol = 1),
       Z = matrix(as.numeric(Z) - mean(Z), ncol = 1),
       Y = matrix(as.numeric(Y) - mean(Y), ncol = 1),
       G = G, H = H, V = V,
       sd_X = sd_X, sd_Z = sd_Z, sd_Y = sd_Y,
       representation = best_name)
}
