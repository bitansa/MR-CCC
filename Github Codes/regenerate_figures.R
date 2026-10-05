############################################################
## regenerate_figures.R
##
## REBUILD EVERY FIGURE FROM SAVED RESULTS -- NO MODEL FITTING
##
## The analyses that produce the figures are expensive: the real-data sweep
## takes roughly an hour per cell-type pair (about a day for all
## twenty ordered pairs), the simulation grids take hours, and the trace
## figures require re-running four chains of 100,000 iterations per
## triplet. Every one of those scripts writes its results to Results/
## before plotting. This script reads those saved results and rebuilds
## every figure through the same plot-building functions the analysis
## scripts use (figure_functions.R, sourced via figure_style.R), so that
## figures can be restyled, re-labelled or re-exported at any time without
## re-running an analysis.
##
## USAGE (from the project root):
##   source("Github Codes/regenerate_figures.R")
##
## INPUTS (all under Results/; any that are absent are skipped with a
## message rather than an error, so the script can be run at any stage of
## the pipeline):
##   <Cell1>_<Cell2>_MR_CCC.rds            -> p1, p2, p3 for every pair found
##   simulation_S1_S3_raw.csv              -> p_score, p_betaX, p_betaXZ
##   simulation_S4_S7_summary.csv          -> p_s4, p_s5, p_s6, p_s7
##   simulation_S{4,5,6,7}_*.csv            -> S4-S7 estimate boxplots
##   convergence_trace_<Cell1>_<Cell2>.rds -> p_trace for every pair found
##
## OUTPUTS (all under Plots/). Apart from the main-text copies, the names are
## identical to those the analysis scripts write, so regenerated figures
## replace the originals in place:
##   Supp_<Cell1>_<Cell2>_bubble.pdf, _pip_ranking.pdf, _curves.pdf
##   sim_score.pdf, sim_betaX.pdf, sim_betaXZ.pdf
##   sim_S4_weak.pdf, sim_S5_correlated.pdf, sim_S6_pleiotropy.pdf,
##   sim_S7_nonlinear.pdf
##   sim_S4_S7_score_box.pdf, sim_S4_S7_betaX_box.pdf,
##   sim_S4_S7_betaXZ_box.pdf
##   convergence_trace_<Cell1>_<Cell2>.pdf
## Main-text copies, written only by this script:
##   sim_score_main.pdf, sim_betaXZ_main.pdf
##   Main_pip_ranking_NKCells_MonocytesCells.pdf
##
## SETTINGS THAT MUST MATCH THE ANALYSIS SCRIPTS
##   The saved real-data results carry both discovery flags (disc_pip and
##   disc_bfdr); DISCOVERY_RULE selects which one drives the figures and
##   pip_thresh is the PIP threshold drawn on them, exactly as in
##   real_data_analysis.R. pip_thr is the rejection threshold used to
##   re-derive the MR-CCC rejection column in the simulation long format,
##   exactly as in simulation_mrccc.R.
############################################################
source("Github Codes/figure_style.R")

DISCOVERY_RULE <- "pip"
pip_thresh     <- 0.5     # real data: PIP threshold shown on the figures
pip_thr        <- 0.5     # simulation: MR-CCC rejection threshold

## Reports a missing input and returns FALSE, so each block below can be
## skipped cleanly.
input_present <- function(path) {
  if (file.exists(path)) return(TRUE)
  message("Skipping: ", path, " not found.")
  FALSE
}

############################################################
## 1. REAL DATA: p1 (bubble), p2 (PIP ranking), p3 (effect curves)
##
## One set of three figures per ordered cell-type pair. Every
## Results/<Cell1>_<Cell2>_MR_CCC.rds is processed, so a complete sweep
## regenerates all twenty pairs.
############################################################
## Results of smoke-test runs are named with a "smoke_" prefix by
## real_data_analysis.R and are excluded here: they come from a truncated
## triplet universe and a shortened chain, so figures built from them would
## not be reportable.
rds_files <- sort(list.files("Results", pattern = "_MR_CCC\\.rds$",
                             full.names = TRUE))
.is_smoke <- startsWith(basename(rds_files), "smoke_")
if (any(.is_smoke)) {
  message("Skipping ", sum(.is_smoke),
          " smoke-test result file(s): ",
          paste(basename(rds_files[.is_smoke]), collapse = ", "))
  rds_files <- rds_files[!.is_smoke]
}
if (length(rds_files) == 0L) {
  message("Skipping real-data figures: no Results/*_MR_CCC.rds found.")
} else {
  cat("\n=== Real-data figures:", length(rds_files), "pair(s) ===\n")
  for (f in rds_files) {
    cells <- strsplit(sub("_MR_CCC\\.rds$", "", basename(f)), "_")[[1]]
    if (length(cells) != 2L) {
      message("Skipping ", f, ": file name does not encode <Cell1>_<Cell2>.")
      next
    }
    Cell1 <- cells[1]; Cell2 <- cells[2]
    cat("--", Cell1, "->", Cell2, "\n")

    # A pair whose saved object cannot be plotted (for example one written
    # without the discovery columns) is reported and skipped, so that one
    # bad file does not stop the remaining pairs.
    tryCatch({
      out <- readRDS(f) %>% add_plot_columns(DISCOVERY_RULE)
      pair_label <- pair_label_for(Cell1, Cell2)

      p1 <- plot_bubble(out, pair_label, pip_thresh)
      p2 <- plot_pip_ranking(out, pair_label, pip_thresh)
      p3 <- plot_effect_curves(out, pair_label, pip_thresh)

      # The caption names the pair and explains the encoding, so the in-plot
      # title is dropped to give the 42 rows the height they need.
      save_fig(p1 + labs(title = NULL, subtitle = NULL),
               paste0("Supp_", Cell1, "_", Cell2, "_bubble"),      12, 16)
      save_fig(p2, paste0("Supp_", Cell1, "_", Cell2, "_pip_ranking"), 14, 13)
      save_fig(p3, paste0("Supp_", Cell1, "_", Cell2, "_curves"),      14, 9)

      # Main-text copy of the focal pair's PIP ranking. The main text prints
      # it at about five inches, so it is drawn on a smaller canvas (the same
      # aspect ratio, hence the same height on the page) and the title, which
      # the caption already carries, is dropped; both make the text print
      # larger.
      if (Cell1 == "NKCells" && Cell2 == "MonocytesCells") {
        p2_main <- p2 + labs(title = NULL, subtitle = NULL,
                             color = "Pathway") +
          theme(axis.text.y     = element_text(size = 14),
                axis.title      = element_text(size = 20, face = "bold"),
                legend.position = "right",
                legend.title    = element_text(size = 18, face = "bold"),
                legend.text     = element_text(size = 16)) +
          guides(color = guide_legend(ncol = 1,
                                      override.aes = list(size = 4)))
        save_fig(p2_main, "Main_pip_ranking_NKCells_MonocytesCells",
                 10, 9.3)
      }
    }, error = function(e) {
      message("Skipping ", f, ": ", conditionMessage(e))
    })
  }
}

############################################################
## 2. SIMULATION S1-S3: p_score, p_betaX, p_betaXZ
############################################################
sim_raw_path <- "Results/simulation_S1_S3_raw.csv"
if (input_present(sim_raw_path)) {
  cat("\n=== Simulation S1-S3 figures ===\n")
  sim_raw  <- read.csv(sim_raw_path, stringsAsFactors = FALSE)
  sim_long <- make_long_results(sim_raw, pip_thr = pip_thr)

  p_score  <- plot_sim_score(sim_long)
  p_betaX  <- plot_sim_betaX(sim_long)
  p_betaXZ <- plot_sim_betaXZ(sim_long)

  save_fig(p_score,  "sim_score",  14, 9)
  save_fig(p_betaX,  "sim_betaX",  14, 9)
  save_fig(p_betaXZ, "sim_betaXZ", 14, 9)

  # Main-text copies (Figure 2): two panels side by side, each about three
  # inches wide in print. Same aspect ratio on a smaller canvas, without the
  # in-plot titles that the caption already carries.
  # The facets are too narrow for four rotated method names, so the main
  # copies identify the methods by fill colour with one legend instead.
  main_sim <- function(p) {
    p + labs(title = NULL, subtitle = NULL, x = NULL) +
      scale_fill_manual(values = pal_methods, name = NULL) +
      theme(axis.text.x     = element_blank(),
            axis.ticks.x    = element_blank(),
            legend.position = "bottom",
            legend.text     = element_text(size = 22))
  }
  save_fig(main_sim(p_score),  "sim_score_main",  9.8, 6.3)
  save_fig(main_sim(p_betaXZ), "sim_betaXZ_main", 9.8, 6.3)
}

############################################################
## 3. SIMULATION S4-S7: p_s4, p_s5, p_s6, p_s7
############################################################
sum_mis_path <- "Results/simulation_S4_S7_summary.csv"
if (input_present(sum_mis_path)) {
  cat("\n=== Simulation S4-S7 figures ===\n")
  sum_mis <- read.csv(sum_mis_path, stringsAsFactors = FALSE)

  p_s4 <- plot_s4_weak(sum_mis)
  p_s5 <- plot_s5_correlated(sum_mis)
  p_s6 <- plot_s6_pleiotropy(sum_mis)
  p_s7 <- plot_s7_nonlinear(sum_mis)

  save_fig(p_s4, "sim_S4_weak",       15, 10)
  save_fig(p_s5, "sim_S5_correlated", 15, 10)
  save_fig(p_s6, "sim_S6_pleiotropy", 18, 10)
  save_fig(p_s7, "sim_S7_nonlinear",  14, 10)
}

############################################################
## 3b. SIMULATION S4-S7: estimate-level boxplots
##
## Built from the per-replicate CSVs rather than the summary, since the
## boxplots need the individual estimates.
############################################################
mis_raw_paths <- file.path("Results",
                           c("simulation_S4_weak.csv",
                             "simulation_S5_correlated.csv",
                             "simulation_S6_pleiotropy.csv",
                             "simulation_S7_nonlinear.csv"))
if (all(vapply(mis_raw_paths, input_present, logical(1)))) {
  cat("\n=== Simulation S4-S7 estimate boxplots ===\n")
  mis_raw <- dplyr::bind_rows(
    lapply(mis_raw_paths, read.csv, stringsAsFactors = FALSE))

  save_fig(plot_mis_estimate_box(mis_raw, "score"),
           "sim_S4_S7_score_box",  16, 11)
  save_fig(plot_mis_estimate_box(mis_raw, "beta_X"),
           "sim_S4_S7_betaX_box",  16, 11)
  save_fig(plot_mis_estimate_box(mis_raw, "beta_XZ"),
           "sim_S4_S7_betaXZ_box", 16, 11)
}

############################################################
## 4. CONVERGENCE: p_trace
##
## One trace figure per Results/convergence_trace_<Cell1>_<Cell2>.rds.
############################################################
trace_files <- sort(list.files("Results", pattern = "^convergence_trace_.*\\.rds$",
                               full.names = TRUE))
if (length(trace_files) == 0L) {
  message("Skipping trace plots: no Results/convergence_trace_*.rds found.")
} else {
  cat("\n=== Convergence trace figures:", length(trace_files), "pair(s) ===\n")
  for (f in trace_files) {
    pair_tag <- sub("^convergence_trace_(.*)\\.rds$", "\\1", basename(f))
    tryCatch({
      trace_tab <- readRDS(f)
      p_trace   <- plot_trace(trace_tab)
      save_fig(p_trace, paste0("convergence_trace_", pair_tag), 14, 12)
    }, error = function(e) {
      message("Skipping ", f, ": ", conditionMessage(e))
    })
  }
}

cat("\nFigure regeneration complete.\n")
