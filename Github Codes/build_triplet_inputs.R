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
## single definition removes that whole class of
## inconsistency.
##
## KEEP IN SYNC with Steps 1-6 of real_data_analysis.R. If
## the primary pipeline changes how a triplet is assembled,
## this file must change with it.
##
## PREREQUISITES (all inherited from real_data_analysis.R):
##   lr_filtered, key_cols, RNA.count.adj_Cell1,
##   RNA.count.adj_Cell2, V_donor, get_SNP_matrix, gsva_scores,
##   centre_cols, orient_pc1
############################################################

build_triplet_inputs <- function(lig_sym, rec_sym, pth, verbose = TRUE) {

  run <- which(lr_filtered$ligand_symbol   == lig_sym &
               lr_filtered$receptor_symbol == rec_sym &
               lr_filtered$pathway_name    == pth)[1]
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
