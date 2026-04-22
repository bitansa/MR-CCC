# MR-CCC: Bayesian Mendelian Randomization for Causal Cell-Cell Communication

Code for the paper:

> **Bayesian Mendelian Randomization for Causal Cell-Cell Communication**  
> Bitan Sarkar and Yang Ni, 2026

---

## Repository Contents

### `mr_ccc_gibbs.cpp`
The core Gibbs sampler implementing the MR-CCC model. Written in C++ using RcppArmadillo. Compile in R with:
```r
Rcpp::sourceCpp("mr_ccc_gibbs.cpp")
```
This makes the function `mr_ccc_gibbs()` available in your R session.

### `build_lr_database.R`
Builds the ligand-receptor (LR) pair database by combining:
- **CellPhoneDB v5** — curated LR interactions with gene symbols and pathway annotations. Download the required input files from: https://github.com/ventolab/cellphonedb-data/tree/master/data
- **MSigDB Reactome (C2/CP:REACTOME)** — pathway gene sets, accessed via the `msigdbr` R package (no download needed).

Set `DATA_DIR` at the top of the script to the folder containing the CellPhoneDB CSV files, then `source()` the script. Produces two objects: `lr_matrix` and `key_cols`.

### `real_data_analysis.R`
Runs the MR-CCC analysis for one ordered sender → receiver cell-type pair, iterating over all ligand-receptor-pathway triplets. Also produces three publication figures (bubble plot, PIP lollipop, effect curves).

**To change the cell-type pair**, edit the two lines at the top of the USER SETTINGS block:
```r
Cell1 <- "NKCells"         # sender cell type
Cell2 <- "MonocytesCells"  # receiver cell type
```
Valid values for both Cell1 and Cell2:
| Value | Cell type |
|-------|-----------|
| `"BCells"` | B cells |
| `"CD4Cells"` | CD4⁺ T cells |
| `"CD8Cells"` | CD8⁺ T cells |
| `"NKCells"` | NK cells |
| `"MonocytesCells"` | Monocytes |

Any ordered combination is valid, e.g. `Cell1 = "BCells"`, `Cell2 = "CD4Cells"` runs the B cells → CD4⁺ T cells analysis.

**Prerequisites before sourcing:**
```r
Rcpp::sourceCpp("mr_ccc_gibbs.cpp")
source("build_lr_database.R")
```

### `simulation_mrccc.R`
Simulation study comparing four methods (OLS, MVMR, MR-BMA, MR-CCC) across three scenarios and four sample sizes. Produces summary tables and boxplot figures.

**Prerequisites before sourcing:**
```r
Rcpp::sourceCpp("mr_ccc_gibbs.cpp")
```

---

## Data

Preprocessed OneK1K single-cell RNA-seq and donor data used in the real-data analysis are available at:

> https://doi.org/10.5281/zenodo.19675075

Download `B_T_NK_monocytes.rda` and `donor.rda`, save them to a local folder, and set `DATA_DIR` at the top of `real_data_analysis.R` to that folder path.

---

## Dependencies

```r
# CRAN
install.packages(c(
  "Rcpp", "RcppArmadillo",
  "dplyr", "tidyr", "ggplot2", "forcats", "scales",
  "purrr", "tibble", "stringr", "ggrepel", "ggtext",
  "msigdbr", "parallel"
))

# Bioconductor
BiocManager::install(c("AUCell", "UCell", "GSVA", "GenomicRanges"))
```

A C++ compiler compatible with C++11 is required (Xcode Command Line Tools on macOS; Rtools on Windows).

---

## License

MIT License. See `LICENSE` for details.
