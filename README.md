# MR-CCC: reproducibility repository

[![DOI](https://zenodo.org/badge/1216664035.svg)](https://doi.org/10.5281/zenodo.23167575)

This repository contains the scripts that reproduce every analysis, table
and figure in

> Sarkar B, Ni Y (2026). *MR-CCC: Bayesian Mendelian Randomization for
> Causal Cell–Cell Communication.* arXiv:2604.23917.
> <https://arxiv.org/abs/2604.23917>

The method itself is distributed as the R package **MRCCC**
(<https://github.com/bitansa/MRCCC>), which provides a validated user-facing
interface, documentation, a vignette and a test suite. This repository is
deliberately separate: it holds the paper's analysis pipeline, its saved
results, and the code that draws its figures, so that the package can be
installed and used without any of that, and so that the paper's numbers can
be regenerated without installing the package.

## Layout

```
MR-CCC/
├── Github Codes/          all scripts (R, one C++ sampler, one shell sweep)
├── Results/               saved outputs of every analysis (.rds, .csv)
├── Plots/                 every figure, as PDF
├── Data/                  NOT tracked; see "Data" below
└── README.md, LICENSE
```

Every script is written to be run **from the repository root** with paths of
the form `"Github Codes/<script>.R"`. Do not `setwd()` into `Github Codes/`:
the scripts write to `Results/` and `Plots/` relative to the root, and
sourcing from inside the folder would place those directories in the wrong
location.

## Data

Two kinds of input are required, neither of which is tracked in this
repository.

**Processed expression and genotype data** are archived on Zenodo under DOI
[10.5281/zenodo.19675075](https://doi.org/10.5281/zenodo.19675075):

| File | Contents |
|---|---|
| `B_T_NK_monocytes.rda` | pseudo-bulk expression matrices (genes × donors) for five cell types |
| `donor.rda` | donor metadata, SNP genotype matrix, and `GRanges` objects for SNP and gene coordinates |

Save both to one folder and point the scripts at it with

```r
Sys.setenv(MRCCC_DATA_DIR = "/path/to/folder")      # this session only
```

or a permanent line `MRCCC_DATA_DIR=/path/to/folder` in `~/.Renviron`. When
the variable is unset the scripts look in `~/data/mrccc`.

**CellPhoneDB v5** ligand–receptor tables are downloaded from
<https://github.com/ventolab/cellphonedb-data/tree/master/data>. Four files
are needed — `interaction_input.csv`, `protein_input.csv`,
`complex_input.csv`, `gene_input.csv` — saved to one folder and located with
`CELLPHONEDB_DIR` in the same way. A fifth,
`transcription_factor_input.csv`, is read if present but is not used.

## Requirements

R (≥ 4.4) with a working C++ toolchain (Xcode command-line tools on macOS,
Rtools on Windows, `build-essential` on Debian-based Linux), and the
following packages.

From CRAN:
`Rcpp`, `RcppArmadillo`, `dplyr`, `tidyr`, `tibble`, `purrr`, `stringr`,
`forcats`, `ggplot2`, `scales`, `ggrepel`, `ggtext`, `Matrix`, `msigdbr`
(≥ 10.0.0, which provides the `collection`/`subcollection` arguments used by
`build_lr_database.R`). From version 10, msigdbr may obtain the full MSigDB
gene-set data, including the human Reactome sets used here, from the
companion package `msigdbdf`; if `msigdbr()` reports that `msigdbdf` is
missing, install it as the message directs (at the time of writing,
`install.packages("msigdbdf", repos = "https://igordot.r-universe.dev")`).

From Bioconductor:
`GenomicRanges`, `AUCell`, `UCell`, `GSVA`.

```r
install.packages(c("Rcpp", "RcppArmadillo", "dplyr", "tidyr", "tibble",
                   "purrr", "stringr", "forcats", "ggplot2", "scales",
                   "ggrepel", "ggtext", "Matrix", "msigdbr", "BiocManager"))
BiocManager::install(c("GenomicRanges", "AUCell", "UCell", "GSVA"))
```

## Scripts

| Script | Purpose |
|---|---|
| `mr_ccc_gibbs.cpp` | The blocked Gibbs sampler (RcppArmadillo). Identical to `MRCCC/src/mr_ccc_gibbs.cpp` in the package. |
| `build_lr_database.R` | Builds the ligand–receptor database from CellPhoneDB and loads the Reactome gene sets. |
| `build_triplet_inputs.R` | Shared definition of how one ligand–receptor–pathway triplet is assembled; sourced by the refit scripts so that they use the same input construction as the primary analysis. It also warns when the pathway representation chosen in a refit differs from the one recorded in the primary output. |
| `real_data_analysis.R` | **Primary analysis** for one ordered cell-type pair: all triplets, posterior summaries, convergence diagnostics, error-rate control, three figures. |
| `run_pair.R`, `run_all_pairs.sh` | Run the primary analysis for one pair in a fresh R process, or sweep all twenty ordered pairs. Resumable. |
| `convergence_diagnostics.R` | Trace figure (joint log-likelihood and running inclusion probability, all chains) for the declared discoveries, per-parameter R̂ table, and runtime benchmark. |
| `posterior_sign_probabilities.R` | Posterior sign probabilities of the interaction effect and intervals conditional on inclusion, for the discoveries. |
| `leave_one_instrument_out.R` | Refits each discovery dropping one instrument at a time. |
| `instrument_overlap_sensitivity.R` | Instrument-set diagnostics across all twenty pairs (number of instruments 5/10/20, genomic overlap between ligand and receptor windows, dosage correlation within and between instrument sets, a screen for instruments associated with the partner gene) and, for the main-text pair, the window sensitivity (±100/200/500 kb), the first-stage F profile under single-instrument removal, and a refit on disjoint instruments. |
| `pathway_representation_sensitivity.R` | Refits under each of the six pathway-activity representations against the adaptive rule, and calibrates the adaptive rule against a pre-specified one on permuted outcomes. |
| `simulation_mrccc.R` | Simulation scenarios S1–S3 (correctly specified model). |
| `simulation_misspecification.R` | Simulation scenarios S4–S7 (weak, correlated and pleiotropic instruments; nonlinear exposure effect). |
| `figure_style.R`, `figure_functions.R` | Single source of every figure's appearance (theme, palettes, notation) and one plot-building function per figure. |
| `regenerate_figures.R` | Rebuilds every figure from `Results/` with no model fitting, including the main-text copies (`sim_score_main.pdf`, `sim_betaXZ_main.pdf`, `Main_pip_ranking_NKCells_MonocytesCells.pdf`) sized for the two-column journal layout. |

## Running the analyses

All commands assume the working directory is the repository root.

**1. Compile the sampler and build the database** (once per session):

```r
Rcpp::sourceCpp("Github Codes/mr_ccc_gibbs.cpp")
source("Github Codes/build_lr_database.R")
```

**2. Primary analysis** for the main-text pair, interactively:

```r
.mrccc_cell1 <- "NKCells"; .mrccc_cell2 <- "MonocytesCells"
source("Github Codes/real_data_analysis.R")
```

or for all twenty ordered pairs, non-interactively (each pair in its own R
process; this writes `Results/<Cell1>_<Cell2>_MR_CCC.rds` but leaves no
objects in an R session):

```bash
bash "Github Codes/run_all_pairs.sh" 2>&1 | tee Results/run_all_pairs.log
```

**3. Post-hoc analyses** on the main-text pair. These four use the objects
that the interactive run of step 2 leaves in the session, so they must be
sourced after that run, in the same session (the shell sweep is not enough):

```r
source("Github Codes/convergence_diagnostics.R")
source("Github Codes/posterior_sign_probabilities.R")
source("Github Codes/leave_one_instrument_out.R")
source("Github Codes/pathway_representation_sensitivity.R")
```

`instrument_overlap_sensitivity.R` reloads all data itself and should be run
in a **fresh R session**, after step 1 and once
`Results/NKCells_MonocytesCells_MR_CCC.rds` exists (from either route in
step 2):

```r
Rcpp::sourceCpp("Github Codes/mr_ccc_gibbs.cpp")
source("Github Codes/build_lr_database.R")
source("Github Codes/instrument_overlap_sensitivity.R")
```

**4. Simulations**, independent of the real data:

```r
Rcpp::sourceCpp("Github Codes/mr_ccc_gibbs.cpp")
source("Github Codes/simulation_mrccc.R")
source("Github Codes/simulation_misspecification.R")
```

**5. Figures only**, from saved results, with no model fitting:

```r
source("Github Codes/regenerate_figures.R")
```

### Canonical sampler settings

| Analysis | Chains × iterations | Burn-in | Retained draws |
|---|---|---|---|
| Real data and all sensitivities | 4 × 100,000, dispersed starts | 2,000 per chain | 392,000 pooled |
| Simulations S4–S7 | 1 × 400,000 | 2,000 | 398,000 |
| Simulations S1–S3 | 1 × 20,000 | 2,000 | 18,000 |

Thinning is 1 throughout. The two large analyses rest on the same number of
draws; the real data uses several chains so that the Gelman–Rubin statistic
can be reported for every triplet, and the simulations use one because
between-chain diagnostics on 10,800 independent fits are not reported. The
rationale for each chain length is documented at the top of the
corresponding script. The Monte Carlo standard error of the posterior
inclusion probability is reported for the real-data analyses (every triplet)
and for S4–S7 (every replicate, column `MR_pip_mcse`); the S1–S3 output does
not store it. Both simulation studies use 100 replicates per cell and report
the Monte Carlo standard error of every rejection rate.

**Seeds.** The real-data scripts set the seed 20260921, naming the
Mersenne-Twister generator explicitly (R's default, so a fresh session is
unaffected). `donor.rda` contains a saved `.Random.seed`, which `load()`
restores, so in `real_data_analysis.R` and `instrument_overlap_sensitivity.R`
the draws made after the data are loaded follow that saved state; it is fixed,
so runs remain reproducible. `convergence_diagnostics.R`,
`posterior_sign_probabilities.R` and `leave_one_instrument_out.R` seed each
chain separately. S1–S3 seed each of
the 1,200 replicate jobs individually (seeds 1 to 1,200, passed through
`run_one_replicate()` to `generate_data()`), so they are reproducible
regardless of parallelisation. S4–S7 run in parallel with `mclapply` from
the seed 20260921 and use the parallel-safe L'Ecuyer-CMRG generator
(`RNGkind("L'Ecuyer-CMRG")` in `run_misspecification()`), so a rerun is
reproducible from the seed for a given number of cores. The shipped S4–S7
results reproduce up to Monte Carlo error, not digit for digit. The
generator setting persists for the rest of the R session, so run the
simulations in a session separate from the real-data analyses.

### Preprocessing

Cell-level counts are summed to donor-level pseudo-bulk profiles. Donors whose
library size in either cell type is more than three median absolute deviations
above the median, or below a quarter of it, are excluded; each profile is then
divided by its library-size factor (relative to the median). These normalised
counts are used for ligand and receptor expression, the eQTL screen and all six
pathway scores. Genotypes are coded as
dosages 0/1/2; SNPs carried by at most 5% of donors with a called genotype are
removed, and remaining missing calls are imputed by the SNP's mean dosage. The
donor covariates are sex, age in ten-year bins and three genotype principal
components; the same covariate matrix is used in the eQTL screen, the
first-stage F statistic and the model. Because neither first stage of the model
has an intercept, X, Z, Y and every column of the instrument and covariate
matrices are centred before fitting.

### Instrument selection and pathway representation

For each gene, candidate SNPs are taken from a ±200 kb window around the
promoter, ranked by their covariate-adjusted marginal cis-eQTL p-value, and
retained in that order — skipping any SNP whose squared dosage correlation
with one already retained exceeds 0.8 — until ten are kept (`LD_R2_MAX` in
`real_data_analysis.R`; the same rule is applied in
`instrument_overlap_sensitivity.R`, and the threshold matches the largest
correlation examined in simulation S5). Pathway activity in the receiver cell
type is summarised by six representations (PC1, AUCell, mean, UCell, ssGSEA,
GSVA); the primary analysis uses, for each triplet, the representation most
correlated in absolute value with ligand expression (PC1 is oriented so that
higher scores mean higher mean expression of the pathway genes), and
`pathway_representation_sensitivity.R` reports the results under each fixed
representation and the false-positive rate of the adaptive rule on permuted
outcomes.

### Approximate runtimes

| Step | Time |
|---|---|
| Primary analysis, one pair | ~1–2 h (pathway scoring dominates; ~41 of ~151 candidate triplets are analysed) |
| All twenty pairs | ~26 h |
| Convergence diagnostics | ~3–4 min per discovery |
| Posterior sign probabilities | ~3–4 min per discovery (one refit of 4 chains, ~1 min, plus pathway scoring) |
| Leave-one-instrument-out | ~20 min per discovery |
| Pathway-representation sensitivity, main-text pair | ~12 h: Part 1 refits ~41 triplets under 7 rules (~290 refits) and Part 2 runs 10 permutations × 20 triplets × 2 rules (400 refits), each refit ~1 min |
| Instrument-overlap sensitivity | ~5 min for Tier 1 (all twenty pairs, no MCMC); ~2–3 h with Tier 2 (~41 refits plus pathway scoring) |
| Simulations S1–S3 (1,200 fits × 20,000 sweeps, with comparators) | ~1.7 h |
| Simulations S4–S7 (10,800 fits × 400,000 sweeps) | ~21 h |

These are approximate. Sampling costs about 14 s per 100,000 sweeps at
n = 651, so one refit of four chains takes about a minute; rebuilding the six
pathway scores of a triplet takes a further 2–3 minutes.

## Outputs

Every analysis writes its results to `Results/` **before** drawing anything,
so figures can be restyled or re-exported from the saved objects.

| File | Written by |
|---|---|
| `Results/<Cell1>_<Cell2>_MR_CCC.rds` | `real_data_analysis.R` — one row per triplet with posterior means, credible intervals, sign-reversal threshold, MCSE, ESS, R̂ (β_X, β_XZ, γ, log-likelihood), instrument strength, discovery flags, chosen pathway representation, and the ligand and receptor Ensembl identifiers (`ligand_ensembl`, `receptor_ensembl`; labels only, not required by any downstream script) |
| `Results/convergence_diagnostics_<C1>_<C2>.csv`, `Results/convergence_trace_<C1>_<C2>.rds` | `convergence_diagnostics.R` |
| `Results/runtime_benchmark.csv` | `convergence_diagnostics.R` |
| `Results/posterior_sign_probabilities_<C1>_<C2>.csv` | `posterior_sign_probabilities.R` |
| `Results/leave_one_out_{summary,detail}_<C1>_<C2>.csv` | `leave_one_instrument_out.R` |
| `Results/overlap_{detail,summary}_all_pairs.csv`, `Results/instrument_count_sensitivity.csv`, `Results/selection_variant_{summary,detail}.csv`, `Results/leave_one_out_F_profile_<C1>_<C2>.csv`, `Results/sensitivity_<C1>_<C2>.csv` | `instrument_overlap_sensitivity.R` |
| `Results/representation_{fixed,concord,null}_<C1>_<C2>.csv` | `pathway_representation_sensitivity.R` |
| `Results/simulation_S1_S3_raw.csv`, `Results/simulation_table_S{1,2,3}.csv` | `simulation_mrccc.R` |
| `Results/simulation_S{4,5,6,7}_*.csv`, `Results/simulation_S4_S7_summary.csv` | `simulation_misspecification.R` |
| `Results/sessionInfo_<script>.txt` | `real_data_analysis.R`, the four post-hoc scripts of step 3, `instrument_overlap_sensitivity.R`, `simulation_mrccc.R`, `simulation_misspecification.R` |
| `Plots/*.pdf` | the analysis scripts and `regenerate_figures.R` (identical names, so regenerated figures replace the originals) |

**Session information.** Each of these scripts ends by writing
`utils::capture.output(sessionInfo())` to `Results/sessionInfo_<script>.txt`,
which records the R and package versions used for that run. Record the
CellPhoneDB release (v5) and the msigdbr version (with its MSigDB release)
alongside these files: the ligand–receptor table and the Reactome gene sets
depend on both, and the CellPhoneDB tables are read from CSV files that no
package version identifies.

## Smoke test

`real_data_analysis.R` accepts two overrides for a quick end-to-end check:

```r
.mrccc_n_iter       <- 6000   # short chains; burn-in scales down automatically
.mrccc_max_triplets <- 3      # truncate the triplet universe
source("Github Codes/real_data_analysis.R")
```

Both are announced loudly at the start of the run and again beside the
discovery counts, and the chain length is recorded in the output object, so a
check can never be mistaken for a reportable result. Remove both objects
(`rm(.mrccc_n_iter, .mrccc_max_triplets)`) before a real run.

A check writes its outputs under names prefixed `smoke_` (for example
`Results/smoke_NKCells_MonocytesCells_MR_CCC.rds`), so that it never
overwrites a reported result. The prefix is applied whenever `N_ITER` is not
100,000, the triplet universe is truncated, or the clumping threshold is not
0.8. The post-hoc scripts of step 3 inherit it from the session, and
`instrument_overlap_sensitivity.R` applies it when its clumping threshold is
not 0.8. `regenerate_figures.R` and `run_all_pairs.sh` ignore `smoke_` files.

### Overrides

Each override is an object assigned in the session **before** the script is
sourced; none requires editing a script.

| Object | Read by | Effect |
|---|---|---|
| `.mrccc_cell1`, `.mrccc_cell2` | `real_data_analysis.R` | Sender and receiver cell types (default NK cells → monocytes); set by `run_pair.R`. |
| `.mrccc_n_iter`, `.mrccc_max_triplets` | `real_data_analysis.R` | Smoke test, as above. |
| `.mrccc_ld_r2_max` | `real_data_analysis.R`, `instrument_overlap_sensitivity.R` | LD clumping threshold r² for instrument selection (default 0.8; 1 disables clumping). Any other value is treated as a check and prefixes the outputs `smoke_`. |
| `.mrccc_run_tier2` | `instrument_overlap_sensitivity.R` | `FALSE` runs the diagnostics (Tier 1 and Sections 4b–4e) without the Tier 2 refit (default `TRUE`). |
| `.mrccc_sim_cores` | `simulation_mrccc.R`, `simulation_misspecification.R` | Caps the number of cores used by `mclapply` (default: all cores but one). |
| `.run_s1_s3` | `simulation_mrccc.R` | `FALSE` loads the functions and settings without running S1–S3, e.g. before `simulation_misspecification.R` in a new session (default `TRUE`). |

## Citation

Archived releases of this repository are on Zenodo:
[10.5281/zenodo.23167575](https://doi.org/10.5281/zenodo.23167575)
(this DOI always resolves to the latest version).

```bibtex
@article{sarkar2026mrccc,
  title   = {{MR-CCC}: Bayesian Mendelian Randomization for Causal Cell--Cell Communication},
  author  = {Sarkar, Bitan and Ni, Yang},
  journal = {arXiv preprint arXiv:2604.23917},
  year    = {2026}
}
```

## License

MIT. See `LICENSE`.
