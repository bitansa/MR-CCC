############################################################
## real_data_analysis.R
##
## MR-CCC real data analysis: one ordered cell-type pair
## (Cell1 -> Cell2), iterating over all ligand-receptor-
## pathway triplets in lr_filtered.
##
## PREREQUISITES (run before sourcing this file):
##   Rcpp::sourceCpp("mr_ccc_gibbs.cpp")
##   source("build_lr_database.R")   # produces lr_matrix, key_cols
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
##                       PIP (gamma_mean), Z_vec, Effect_vec
##   Three ggplot objects: p1 (bubble), p2 (lollipop), p3 (curves)
##   PDF files saved to Plots/ subdirectory.
##   RDS file saved to Results/ subdirectory for re-use.
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

############################################################
## USER SETTINGS
## Set Cell1 (sender) and Cell2 (receiver) before running.
##
## Valid values:
##   "BCells", "CD4Cells", "CD8Cells", "NKCells", "MonocytesCells"
############################################################
Cell1 <- "NKCells"         # sender cell type
Cell2 <- "MonocytesCells"  # receiver cell type
pip_thresh <- 0.5          # PIP threshold for communication discovery

# Root directory where the two Zenodo .rda files are saved.
# Update this path before running.
DATA_DIR <- file.path("~", "data", "mrccc")

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
## that all donors are on a comparable scale.
############################################################
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
cat("After QC:", ncol(RNA.count.adj_Cell1), "donors retained.\n")

############################################################
## SECTION 3: GENOTYPE PROCESSING
##
## The raw genotype matrix uses string encoding ("0/0", "0/1",
## "1/1", "./."). Convert to a numeric dosage matrix (0/1/2),
## then filter out rare SNPs (minor-allele frequency < 5%,
## i.e. fewer than 5% of donors carry at least one copy).
############################################################
# Promote gene reference to promoter windows (+/-200 kb)
# used for cis-eQTL instrument selection.
gene.promoter.ref <- promoters(gene.ref,
                               upstream   = 200000,
                               downstream = 200000)

# Convert string genotypes to numeric dosage
genotype.mat <- matrix(nrow = nrow(genotype), ncol = ncol(genotype))
genotype.mat[genotype == "./."] <- 0
genotype.mat[genotype == "0/0"] <- 0
genotype.mat[genotype == "0/1"] <- 1
genotype.mat[genotype == "1/1"] <- 2

# Minor-allele frequency filter: keep SNPs with > 5% carrier rate
SNP.filter   <- apply(genotype.mat, 1, function(x) mean(x > 0)) > 0.05
cat("SNPs retained after MAF filter:", sum(SNP.filter), "\n")
genotype     <- genotype[SNP.filter, ]
genotype.mat <- genotype.mat[SNP.filter, ]
vcf.ref      <- vcf.ref[SNP.filter]

############################################################
## SECTION 4: FILTER LR PAIRS TO GENES PRESENT IN DATA
##
## lr_matrix (from build_lr_database.R) contains all
## CellPhoneDB Ligand-Receptor pairs with Ensembl IDs.
## We keep only pairs where both genes appear as rows in
## the count matrix for the sender (Cell1) cell type.
## key_cols (MSigDB Reactome) is filtered similarly.
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

############################################################
## SECTION 5: HELPER FUNCTIONS
############################################################
# Build the cis-eQTL instrument matrix G for a given gene.
# Scans the gene's promoter window (+/-200 kb) for SNPs,
# ranks them by eQTL p-value, and returns the top SNPs
# (up to 10) as columns of the instrument matrix.
# Returns NULL if fewer than min_snps qualifying SNPs exist.
get_SNP_matrix <- function(gene_id, RNA_mat_adj, min_snps = 5) {
  gene.ind <- which(gene.ref$gene_id == gene_id)
  if (length(gene.ind) == 0) return(NULL)
  
  snp.ind <- which(
    countOverlaps(vcf.ref, gene.promoter.ref[gene.ind]) > 0
  )
  if (length(snp.ind) < min_snps) return(NULL)
  
  # Compute marginal eQTL p-value for each candidate SNP
  snp.pval <- rep(NA_real_, length(snp.ind))
  for (j in seq_along(snp.ind)) {
    lm.j <- lm(
      RNA_mat_adj[gene_id, ] ~ genotype.mat[snp.ind[j], ] +
        donor$age + donor$sex + donor$PC1 + donor$PC2 + donor$PC3
    )
    snp.pval[j] <- summary(lm.j)$coefficients[2, 4]
  }
  
  ok <- which(!is.na(snp.pval))
  if (length(ok) < min_snps) return(NULL)
  
  # Keep the top SNPs by eQTL p-value (up to 10)
  keep_k    <- min(10L, length(ok))
  best_snps <- snp.ind[ok[order(snp.pval[ok])[1:keep_k]]]
  t(genotype.mat[best_snps, , drop = FALSE])  # n x k (donors x SNPs)
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
  
  if (run %% 10 == 0) {
    cat("Processed", run, "of", nrow(lr_filtered), "triplets\n")
  }
  
  pathway_name       <- lr_filtered$pathway_name[run]
  ligand_gene_name   <- lr_filtered$ligand_symbol[run]
  ligand_gene_id     <- lr_filtered$ligand_ensembl[run]
  receptor_gene_name <- lr_filtered$receptor_symbol[run]
  receptor_gene_id   <- lr_filtered$receptor_ensembl[run]
  
  # ---- Step 1: Build pathway activity Y for Cell2 ----
  # Extract keyword from pathway name to find corresponding
  # genes in MSigDB Reactome (key_cols).
  keyword          <- sub(".*by\\s+", "", pathway_name)
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
  
  # (d) UCell: Mann-Whitney U-based enrichment score
  ucell_scores <- ScoreSignatures_UCell(
    expr_full, features = list(pathway = path_genes)
  )
  Y_UCell <- as.numeric(ucell_scores[, "pathway_UCell"])
  
  # (e) ssGSEA: single-sample GSEA enrichment score
  ssgsea_res <- gsva(expr = expr_full, gset.idx.list = geneSets,
                     method = "ssgsea", kcdf = "Gaussian",
                     abs.ranking = TRUE, verbose = FALSE)
  Y_ssGSEA   <- as.numeric(ssgsea_res["pathway", ])
  
  # (f) GSVA: Gaussian kernel GSVA score
  gsva_res <- gsva(expr = expr_full, gset.idx.list = geneSets,
                   method = "gsva", kcdf = "Gaussian", verbose = FALSE)
  Y_GSVA   <- as.numeric(gsva_res["pathway", ])
  
  # Select Y with highest absolute correlation with X
  Y_list    <- list(PC1 = Y_PC1, AUCell = Y_AUCell, Mean = Y_rowMean,
                    UCell = Y_UCell, ssGSEA = Y_ssGSEA, GSVA = Y_GSVA)
  cors      <- sapply(Y_list, function(y) cor(y, X, use = "complete.obs"))
  best_name <- names(which.max(abs(cors)))
  Y         <- matrix(Y_list[[best_name]], ncol = 1)
  
  # ---- Step 4: Standardization scales (pre-centering) ----
  # SDs are unchanged by centering and are used to convert
  # posterior means to interpretable standardized effect sizes.
  sd_X <- sd(as.numeric(X))
  sd_Z <- sd(as.numeric(Z))
  sd_Y <- sd(as.numeric(Y))
  
  # ---- Centering: enforce E[X] = E[Z] = E[Y] = 0 ----
  # Required by the centering assumption in Proposition 1.
  # Beta_X and Beta_XZ are location-invariant; only mu changes.
  X_c <- matrix(as.numeric(X) - mean(X), ncol = 1)
  Z_c <- matrix(as.numeric(Z) - mean(Z), ncol = 1)
  Y_c <- matrix(as.numeric(Y) - mean(Y), ncol = 1)
  
  # ---- Step 5: Instrument matrices G and H ----
  G <- get_SNP_matrix(ligand_gene_id,   RNA.count.adj_Cell1, min_snps = 1)
  if (is.null(G)) next
  H <- get_SNP_matrix(receptor_gene_id, RNA.count.adj_Cell2, min_snps = 1)
  if (is.null(H)) next
  
  # ---- Step 6: Covariate matrix V ----
  # Age is binned into 10-year intervals (1 = <30, 7 = 80+)
  # to regularize its effect given the moderate donor count.
  V <- donor %>%
    mutate(
      MALE = as.integer(sex == "male"),
      AGE  = case_when(
        age < 30 ~ 1L, age < 40 ~ 2L, age < 50 ~ 3L,
        age < 60 ~ 4L, age < 70 ~ 5L, age < 80 ~ 6L,
        TRUE     ~ 7L
      )
    ) %>%
    dplyr::select(MALE, AGE, PC1, PC2, PC3) %>%
    as.matrix()
  
  # ---- Step 7: Run MR-CCC Gibbs sampler ----
  # g-prior scales follow the paper: min(n, 100).
  # Default hyperparameters: a_sigma=3, b_sigma=2,
  # a_rho=3, b_rho=1, nu1=1e-4.
  n       <- nrow(X_c)
  g_scale <- min(n, 100.0)
  
  res <- mr_ccc_gibbs(
    X_c, Z_c, Y_c, G, H, V,
    n_iter  = 20000, burn_in = 2000, thin = 10,
    a_sigma = 3.0,   b_sigma = 2.0,
    a_rho   = 3.0,   b_rho   = 1.0,
    nu1     = 1e-4,
    gG      = g_scale, gH = g_scale,
    gV      = g_scale, gZ = g_scale,
    gBeta   = g_scale,
    ridge   = 1e-8
  )
  
  # ---- Step 8: Standardize effect sizes ----
  # Convert from raw expression units to SD units:
  #   beta_X^(s)  = beta_X  * sd(X) / sd(Y)
  #   beta_XZ^(s) = beta_XZ * sd(X) * sd(Z) / sd(Y)
  #   beta_Z^(s)  = beta_Z  * sd(Z) / sd(Y)
  Beta_X_std  <- res$Beta_X_mean  * sd_X / sd_Y
  Beta_XZ_std <- res$Beta_XZ_mean * sd_X * sd_Z / sd_Y
  Beta_Z_std  <- res$Beta_Z_mean  * sd_Z / sd_Y
  
  # Z_vec: receptor expression on the SD scale (for plotting)
  # Effect_vec: total ligand effect as a function of receptor level
  #   effect(z) = beta_X^(s) + beta_XZ^(s) * z,  z = Z_c / sd(Z)
  z_c        <- as.numeric(Z_c)
  Z_vec      <- z_c / sd_Z
  Effect_vec <- Beta_X_std + Beta_XZ_std * Z_vec
  
  # ---- Step 9: Store output row ----
  out_idx <- out_idx + 1
  Output_list[[out_idx]] <- tibble::tibble(
    ligand_col_name   = ligand_gene_name,
    receptor_col_name = receptor_gene_name,
    pathway_name      = pathway_name,
    Beta_X_mean       = Beta_X_std,
    Beta_XZ_mean      = Beta_XZ_std,
    Beta_Z_mean       = Beta_Z_std,
    gamma_mean        = res$gamma_mean,  # PIP = P(gamma=1 | data)
    Z_vec             = list(Z_vec),
    Effect_vec        = list(Effect_vec)
  )
}

# Combine all output rows into a single data frame
Output_MR_CCC <- dplyr::bind_rows(Output_list[seq_len(out_idx)])

cat("\nDone.", out_idx, "triplets analysed.\n")
cat("Triplets with PIP >=", pip_thresh, ":",
    sum(Output_MR_CCC$gamma_mean >= pip_thresh), "\n")

############################################################
## SECTION 6b: SAVE AND RELOAD RESULTS
##
## Serialise Output_MR_CCC to an RDS file so that plots
## can be regenerated without re-running the Gibbs sampler.
## The filename encodes the sender-receiver pair.
############################################################
dir.create("Results", showWarnings = FALSE)
saveRDS(Output_MR_CCC,
        file = paste0("Results/", Cell1, "_", Cell2, "_MR_CCC.rds"))
Output_MR_CCC <- readRDS(
  paste0("Results/", Cell1, "_", Cell2, "_MR_CCC.rds")
)

############################################################
## SECTION 7: PLOTS
##
## Three publication-quality figures:
##   p1 -- Bubble plot of |beta_X| and |beta_XZ| by pathway
##   p2 -- Ranked lollipop of all triplets by PIP
##   p3 -- Receptor-modulated ligand effect curves
##          (only PIP > pip_thresh pairs displayed)
############################################################

# ---- Helper: HTML-formatted cell type display names ----
fmt_cell <- function(x) {
  dplyr::case_when(
    x == "BCells"         ~ "B cells",
    x == "CD4Cells"       ~ "CD4<sup>+</sup> T cells",
    x == "CD8Cells"       ~ "CD8<sup>+</sup> T cells",
    x == "NKCells"        ~ "NK cells",
    x == "MonocytesCells" ~ "Monocytes",
    TRUE                  ~ x
  )
}
sender     <- fmt_cell(Cell1)
receiver   <- fmt_cell(Cell2)
pair_label <- paste0(sender, " \u2192 ", receiver)

# ---- Prepare shared plotting variables ----
Output_MR_CCC <- Output_MR_CCC %>%
  mutate(
    LR       = paste(ligand_col_name, receptor_col_name, sep = "\u2013"),
    discover = gamma_mean >= pip_thresh
  )

n_disc  <- sum(Output_MR_CCC$discover)
n_total <- nrow(Output_MR_CCC)

# Pathway colour palette (up to 11 distinct pathways)
pathway_levels <- sort(unique(Output_MR_CCC$pathway_name))
pathway_pal <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8"
)
pathway_colors <- setNames(
  pathway_pal[seq_along(pathway_levels)],
  pathway_levels
)

# ---- Shared publication theme ----
theme_pub <- theme_bw(base_size = 15) +
  theme(
    plot.title        = element_markdown(size = 17, face = "bold", hjust = 0),
    plot.subtitle     = element_text(size = 12, color = "grey40", hjust = 0),
    axis.title        = element_text(size = 14, face = "bold"),
    axis.text         = element_text(size = 11, colour = "black"),
    axis.text.x       = element_text(angle = 40, hjust = 1),
    axis.text.y       = element_text(size = 10),
    strip.text        = element_text(size = 12, face = "bold"),
    strip.background  = element_rect(fill = "grey93", color = NA),
    legend.title      = element_text(size = 12, face = "bold"),
    legend.text       = element_text(size = 10),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4)
  )

# ============================================================
# Plot 1 -- Bubble plot of |beta_X| and |beta_XZ|
#
# Each bubble: one (ligand-receptor, pathway) combination.
# Size   = absolute standardized posterior mean effect.
# Fill   = PIP (blue gradient: light = low, dark = high).
# Border = black for PIP > pip_thresh (discovered triplets).
# ============================================================
lr_order <- Output_MR_CCC %>%
  arrange(gamma_mean) %>%
  pull(LR) %>%
  unique()

df_long <- Output_MR_CCC %>%
  mutate(
    Pathway   = fct_reorder(pathway_name, gamma_mean, .fun = max),
    LR_factor = factor(LR, levels = lr_order),
    pip_group = if_else(gamma_mean >= pip_thresh, "high", "low")
  ) %>%
  pivot_longer(
    cols      = c(Beta_X_mean, Beta_XZ_mean),
    names_to  = "effect",
    values_to = "beta"
  ) %>%
  mutate(
    effect   = recode(effect,
                      Beta_X_mean  = "abs(hat(beta)[X]^{(s)})",
                      Beta_XZ_mean = "abs(hat(beta)[XZ]^{(s)})"),
    size_val = abs(beta)
  )

p1 <- ggplot(df_long, aes(x = Pathway, y = LR_factor)) +
  # Low-PIP: faint, grey border
  geom_point(
    data  = df_long %>% filter(pip_group == "low"),
    aes(size = size_val, fill = gamma_mean),
    shape = 21, color = "grey75", stroke = 0.25, alpha = 0.45
  ) +
  # High-PIP: vivid, black border
  geom_point(
    data  = df_long %>% filter(pip_group == "high"),
    aes(size = size_val, fill = gamma_mean),
    shape = 21, color = "black", stroke = 0.7, alpha = 0.95
  ) +
  facet_wrap(~effect, nrow = 1, labeller = label_parsed) +
  scale_fill_gradient(
    low = "#deebf7", high = "#08306b", limits = c(0, 1),
    name = expression(paste("PIP  (", hat(gamma), ")"))
  ) +
  scale_size_continuous(
    range = c(1.5, 9),
    name  = "Standardized\neffect magnitude"
  ) +
  labs(
    title    = paste0("Communication effects: ", pair_label, "  (MR-CCC)"),
    subtitle = paste0(
      "Point size = |standardized posterior mean|;  color = PIP;  ",
      "black border = PIP > ", pip_thresh
    ),
    x = "Pathway",
    y = "Ligand\u2013Receptor"
  ) +
  theme_pub +
  theme(panel.grid.major.y = element_blank())

print(p1)

# ============================================================
# Plot 2 -- Ranked lollipop of all triplets by PIP
#
# All (ligand-receptor-pathway) triplets ranked by PIP.
# Segments and points coloured by pathway.
# Black ring and PIP label on discovered triplets.
# Shaded region = discovery zone (PIP > pip_thresh).
# ============================================================
df_lollipop <- Output_MR_CCC %>%
  mutate(LR_ranked = fct_reorder(LR, gamma_mean))

p2 <- ggplot(df_lollipop, aes(x = gamma_mean, y = LR_ranked)) +
  # Shaded discovery zone
  annotate("rect",
           xmin = pip_thresh, xmax = 1.02, ymin = -Inf, ymax = Inf,
           fill = "#08306b", alpha = 0.04) +
  # Lollipop segments coloured by pathway
  geom_segment(
    aes(x = 0, xend = gamma_mean, yend = LR_ranked,
        color = pathway_name),
    linewidth = 0.65, alpha = 0.75
  ) +
  # Points coloured by pathway
  geom_point(aes(color = pathway_name), size = 3, alpha = 0.9) +
  # Black ring for discovered triplets
  geom_point(
    data  = df_lollipop %>% filter(discover),
    color = "black", size = 4.8, shape = 1, stroke = 1.1
  ) +
  # PIP value label for discovered triplets
  geom_text(
    data  = df_lollipop %>% filter(discover),
    aes(label = sprintf("%.2f", gamma_mean)),
    hjust = -0.30, size = 3.4, fontface = "bold", color = "grey15"
  ) +
  # Threshold line and annotation
  geom_vline(xintercept = pip_thresh,
             linetype = "dashed", color = "grey30", linewidth = 0.8) +
  annotate("text",
           x = pip_thresh - 0.02, y = 1.6,
           label = paste0("PIP = ", pip_thresh),
           hjust = 1, size = 3.5, color = "grey30", fontface = "italic") +
  scale_x_continuous(
    limits = c(0, 1.12),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = c("0", "0.25", "0.50", "0.75", "1.00")
  ) +
  scale_color_manual(values = pathway_colors, name = "Pathway") +
  labs(
    title    = paste0("PIP ranking across all triplets: ",
                      pair_label, "  (MR-CCC)"),
    subtitle = paste0(
      n_total, " ligand\u2013receptor\u2013pathway triplets ranked by PIP;  ",
      n_disc,  " with PIP > ", pip_thresh,
      " (black ring);  shaded region = discovery zone"
    ),
    x = "Posterior inclusion probability (PIP)",
    y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.y        = element_text(size = 8.5),
    axis.text.x        = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.4),
    legend.position    = "right"
  )

print(p2)

# ============================================================
# Plot 3 -- Receptor-modulated ligand effect curves
#           (discovery pairs only, PIP > pip_thresh)
#
# Each curve shows how the total standardized ligand effect
#   effect(z) = beta_X^(s) + beta_XZ^(s) * z
# varies with standardized receptor expression z = Z / SD(Z).
#
# Only pairs that exceed the PIP threshold are plotted.
# Only pathway panels that contain at least one such discovery
# are rendered (pathways with no discoveries are suppressed).
# Every shown curve is labelled at its rightmost data point.
# ============================================================

# Identify pathways that contain at least one discovered pair.
# Panels for pathways with zero discoveries are not drawn.
sig_pathways <- Output_MR_CCC %>%
  filter(gamma_mean >= pip_thresh) %>%
  pull(pathway_name) %>%
  unique()

# Build long-format data for plotting:
#   - restrict to discovered pairs only (PIP > pip_thresh)
#   - restrict to panels that have at least one discovery
#   - unnest Z_vec / Effect_vec into individual (Z, effect) rows
#   - sort by Z within each triplet for correct line drawing
df_lines <- Output_MR_CCC %>%
  filter(pathway_name %in% sig_pathways,
         gamma_mean   >= pip_thresh) %>%
  mutate(
    triplet_id = row_number(),
    paired     = purrr::map2(Z_vec, Effect_vec,
                             ~tibble::tibble(Z = .x, effect = .y))
  ) %>%
  select(triplet_id, pathway_name, LR, gamma_mean, paired) %>%
  tidyr::unnest(paired) %>%
  group_by(triplet_id) %>%
  arrange(Z, .by_group = TRUE) %>%
  ungroup()

# Label data: one label per curve, placed at the rightmost
# observed receptor expression value for that triplet.
# Because all displayed curves are discoveries, every curve
# receives a label.
df_labels <- df_lines %>%
  group_by(triplet_id, pathway_name, LR, gamma_mean) %>%
  slice_max(order_by = Z, n = 1, with_ties = FALSE) %>%
  ungroup()

p3 <- ggplot() +
  
  # Single bold layer: all displayed curves are high-PIP discoveries.
  # Color encodes PIP within the discovery range (0.5--1.0).
  geom_line(
    data      = df_lines,
    aes(x = Z, y = effect, group = triplet_id, color = gamma_mean),
    linewidth = 1.15, alpha = 0.95
  ) +
  
  # Label every curve at its rightmost point.
  geom_text_repel(
    data               = df_labels,
    aes(x = Z, y = effect, label = LR),
    color              = "grey10",
    size               = 3.3,
    fontface           = "bold",
    min.segment.length = 0,
    segment.alpha      = 0.45,
    segment.size       = 0.35,
    max.overlaps       = Inf,
    show.legend        = FALSE
  ) +
  
  # One panel per pathway; only panels with discoveries are shown.
  facet_wrap(~pathway_name, scales = "free") +
  
  scale_color_gradient(
    low    = "#deebf7",
    high   = "#08306b",
    limits = c(0, 1),
    name   = expression(paste("PIP  (", hat(gamma), ")"))
  ) +
  
  labs(
    title    = paste0("Receptor-modulated ligand effects: ",
                      pair_label, "  (MR-CCC)"),
    subtitle = bquote(
      paste("Each curve: ",
            hat(beta)[X]^(s), " + ", hat(beta)[XZ]^(s),
            "  \u00b7  (Z / SD(Z));  ",
            .(n_disc), " pairs with PIP > ", .(pip_thresh), " shown")
    ),
    x = "Standardized receptor expression  (Z / SD(Z))",
    y = "Standardized ligand effect  (SD units of Y per SD of X)"
  ) +
  theme_pub

print(p3)

############################################################
## SECTION 8: SAVE PLOTS
##
## Filenames are auto-built from Cell1 and Cell2 so that
## running the script for different pairs produces separate
## files without overwriting earlier results.
############################################################
dir.create("Plots", showWarnings = FALSE)

ggsave(paste0("Plots/Supp_", Cell1, "_", Cell2, "_bubble.pdf"),
       p1, width = 14, height = 8.5)
ggsave(paste0("Plots/Supp_", Cell1, "_", Cell2, "_pip_ranking.pdf"),
       p2, width = 11, height = 9)
ggsave(paste0("Plots/Supp_", Cell1, "_", Cell2, "_curves.pdf"),
       p3, width = 14, height = 9)

cat("Plots saved to Plots/ directory.\n")
