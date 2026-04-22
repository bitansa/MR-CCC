############################################################
## build_lr_database.R
##
## Builds the ligand-receptor (LR) pair database used by
## MR-CCC, drawing from two sources:
##
##   1. CellPhoneDB v5  -- curated ligand-receptor pairs
##      with gene symbols, UniProt IDs, pathway annotations,
##      and complex membership.
##
##   2. MSigDB Reactome (C2/CP:REACTOME) via msigdbr -- used
##      to map pathway names to gene sets for downstream
##      pathway activity scoring.
##
## OUTPUT:
##   lr_matrix  -- data frame of gene-level LR pairs with
##                 Ensembl IDs and pathway classification
##   key_cols   -- full MSigDB Reactome gene-set table
##
## USAGE:
##   Set DATA_DIR to the folder containing your CellPhoneDB
##   input files, then source this script.
############################################################
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(msigdbr)


############################################################
## USER SETTINGS
############################################################

# CellPhoneDB v5 database files can be downloaded from:
#   https://github.com/ventolab/cellphonedb-data/tree/master/data
#
# Download these five CSV files (click each file, then click "Raw"):
#   interaction_input.csv
#   protein_input.csv
#   complex_input.csv
#   gene_input.csv
#   transcription_factor_input.csv
#
# Save all five to a single local folder, then set DATA_DIR below.

DATA_DIR <- file.path("~", "data", "cellphonedb")
# Examples:
#   macOS / Linux : "~/data/cellphonedb"
#   Windows       : file.path(Sys.getenv("USERPROFILE"), "data", "cellphonedb")


############################################################
## SECTION 1: CellPhoneDB
##
## CellPhoneDB provides a curated set of ligand-receptor
## interactions with UniProt identifiers, gene symbols,
## complex annotations, and pathway classifications.
##
## Input files (all from the CellPhoneDB v5 download):
##   interaction_input.csv -- interaction records
##   protein_input.csv     -- protein metadata
##   complex_input.csv     -- multi-subunit complex definitions
##   gene_input.csv        -- gene symbol <-> Ensembl mapping
##   transcription_factor_input.csv -- TF annotations (optional)
############################################################

# ---- Load raw CellPhoneDB tables ----
inter <- read.csv(file.path(DATA_DIR, "interaction_input.csv"))
prot  <- read.csv(file.path(DATA_DIR, "protein_input.csv"))
comp  <- read.csv(file.path(DATA_DIR, "complex_input.csv"))
genes <- read.csv(file.path(DATA_DIR, "gene_input.csv"))
tfs   <- read.csv(file.path(DATA_DIR, "transcription_factor_input.csv"))  # optional; not used below

# ---- Step 1: Filter to Ligand-Receptor interactions ----
# CellPhoneDB records many interaction types; we keep only
# directional Ligand-Receptor pairs for the CCC analysis.
cols_keep <- c(
  "id_cp_interaction", "partner_a", "partner_b",
  "protein_name_a", "protein_name_b",
  "interactors",
  "classification", "directionality", "modulatory_effect",
  "is_ppi", "source", "curator",
  "reactome_complex", "reactome_reaction", "reactome_pathway",
  "comments"
)

ligrec_raw <- inter %>%
  filter(directionality == "Ligand-Receptor") %>%
  select(any_of(cols_keep))

# ---- Step 2: Rename columns and strip species suffix ----
# partner_a / protein_name_a = ligand side (sender)
# partner_b / protein_name_b = receptor side (receiver)
# _HUMAN suffix is present in CellPhoneDB v5 protein names.
ligrec_clean <- ligrec_raw %>%
  dplyr::rename(
    ligand_uniprot   = partner_a,
    receptor_uniprot = partner_b,
    ligand_name      = protein_name_a,
    receptor_name    = protein_name_b
  ) %>%
  mutate(
    ligand_name   = str_remove(ligand_name,   "_HUMAN$"),
    receptor_name = str_remove(receptor_name, "_HUMAN$")
  )

# ---- Step 3: Extract gene-level identifiers ----
# The `interactors` column encodes gene symbols as
# "GeneA-GeneB" (simple pairs) or "GeneA+GeneB-GeneC+GeneD"
# (multi-subunit complexes, separated by +).
# Fall back to protein names if `interactors` is absent.
ligrec_out <- ligrec_clean %>%
  mutate(
    gene_pair_str = case_when(
      !is.na(interactors) & str_detect(interactors, "-") ~ interactors,
      TRUE ~ paste0(ligand_name, "-", receptor_name)
    )
  ) %>%
  tidyr::separate(
    gene_pair_str,
    into   = c("lig_side", "rec_side"),
    sep    = "-",
    fill   = "right",
    remove = TRUE
  ) %>%
  mutate(
    # Split multi-subunit sides on "+" to get individual gene symbols
    ligand_genes   = str_split(coalesce(lig_side, ligand_name),   "\\+"),
    receptor_genes = str_split(coalesce(rec_side, receptor_name), "\\+")
  ) %>%
  select(
    id_cp_interaction,
    ligand_uniprot,   receptor_uniprot,
    ligand_name,      receptor_name,
    ligand_genes,     receptor_genes,
    classification,   directionality, modulatory_effect,
    is_ppi,           source, curator,
    reactome_complex, reactome_reaction, reactome_pathway,
    comments
  )

# ---- Step 4: Build gene-symbol <-> Ensembl lookup table ----
gene_map <- genes %>%
  transmute(
    symbol  = as.character(hgnc_symbol),
    ensembl = as.character(ensembl)
  ) %>%
  distinct() %>%
  filter(!is.na(symbol), symbol != "")

# ---- Step 5: Parse gene vectors from list columns ----
# Normalises separators and removes empty/malformed tokens.
parse_gene_vec <- function(x) {
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  x   <- as.character(x)
  x   <- str_replace_all(x, "[\\+/:;\\|\\s]+", ",")
  toks <- unlist(str_split(x, ",", simplify = FALSE), use.names = FALSE)
  toks <- toks |> str_trim() |> str_to_upper()
  toks <- toks[toks != "" & toks != "..."]
  unique(toks)
}

# ---- Step 6: Expand to one row per (ligand gene, receptor gene) pair ----
# For multi-subunit complexes this produces all cross-combinations
# of individual subunits, e.g. {A, B} x {C} -> (A,C) and (B,C).
pairs_gene_level <- ligrec_out %>%
  mutate(
    lig_vec = map(ligand_genes,   parse_gene_vec),
    rec_vec = map(receptor_genes, parse_gene_vec)
  ) %>%
  mutate(
    gene_pairs = map2(lig_vec, rec_vec,
                      ~ expand_grid(ligand_symbol   = .x,
                                    receptor_symbol = .y))
  ) %>%
  select(classification, gene_pairs) %>%
  unnest(gene_pairs) %>%
  distinct()

# ---- Step 7: Add Ensembl IDs and produce final LR matrix ----
# Keep only pairs where both ligand and receptor have a
# mapped Ensembl ID (required for eQTL instrument lookup).
gene_map_lig <- gene_map %>% dplyr::rename(ligand_ensembl   = ensembl)
gene_map_rec <- gene_map %>% dplyr::rename(receptor_ensembl = ensembl)

lr_matrix <- pairs_gene_level %>%
  left_join(gene_map_lig,
            by           = c("ligand_symbol" = "symbol"),
            relationship = "many-to-many") %>%
  left_join(gene_map_rec,
            by           = c("receptor_symbol" = "symbol"),
            relationship = "many-to-many") %>%
  transmute(
    ligand_symbol,
    ligand_ensembl,
    receptor_symbol,
    receptor_ensembl,
    pathway_name = classification
  ) %>%
  filter(!is.na(ligand_ensembl), !is.na(receptor_ensembl)) %>%
  distinct()

cat("CellPhoneDB: retained", nrow(lr_matrix),
    "gene-level LR pairs with Ensembl IDs.\n")


############################################################
## SECTION 2: MSigDB Reactome Pathways
##
## msigdbr provides programmatic access to the Molecular
## Signatures Database. Here we load the Reactome collection
## (C2 / CP:REACTOME) to obtain gene-set memberships for
## pathway activity scoring.
##
## The resulting table key_cols contains one row per
## (pathway, gene) entry and is used downstream to define
## the Y variable (pathway activity scores) in MR-CCC.
############################################################

# ---- Load Reactome gene sets for Homo sapiens ----
msig_reactome <- msigdbr(
  species       = "Homo sapiens",
  category      = "C2",
  subcollection = "CP:REACTOME"
)

# Retain the columns needed for downstream pathway scoring
key_cols <- msig_reactome %>%
  select(gs_name, gene_symbol, ensembl_gene, gs_description) %>%
  mutate(gs_name = trimws(as.character(gs_name)))

cat("MSigDB Reactome: loaded", n_distinct(key_cols$gs_name),
    "pathways covering", n_distinct(key_cols$gene_symbol), "genes.\n")

# ---- Quick sanity check: WNT pathways ----
# Confirms that well-known pathway families are present.
wnt_hits <- key_cols %>%
  filter(str_detect(gs_name, regex("WNT", ignore_case = TRUE))) %>%
  distinct(gs_name)

cat("Reactome pathways containing 'WNT':", nrow(wnt_hits), "\n")
