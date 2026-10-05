############################################################
## figure_functions.R
##
## PLOT-BUILDING FUNCTIONS -- ONE DEFINITION PER FIGURE
##
## Each function takes the data frame(s) a figure needs and returns a
## ggplot object. It performs no model fitting and reads no files. The
## analysis scripts call these functions on freshly computed results;
## regenerate_figures.R calls the same functions on results reloaded from
## Results/, so a figure is defined in exactly one place and can be restyled
## or re-exported without re-running any analysis.
##
## This file is sourced by figure_style.R and is not sourced directly:
##   source("Github Codes/figure_style.R")
##
## FUNCTIONS
##   Real data (real_data_analysis.R, regenerate_figures.R)
##     fmt_cell(x)                           display name of a cell type
##     pair_label_for(Cell1, Cell2)          "sender -> receiver" title fragment
##     add_plot_columns(out, discovery_rule) LR label and `discover` flag
##     plot_bubble(out, pair_label, pip_thresh)          p1
##     plot_pip_ranking(out, pair_label, pip_thresh)     p2
##     plot_effect_curves(out, pair_label, pip_thresh)   p3
##
##   Simulation S1-S3 (simulation_mrccc.R, regenerate_figures.R)
##     make_long_results(df_raw, pip_thr)    wide results -> one row per method
##     plot_sim_score(sim_long)              p_score
##     plot_sim_betaX(sim_long)              p_betaX
##     plot_sim_betaXZ(sim_long)             p_betaXZ
##
##   Simulation S4-S7 (simulation_misspecification.R, regenerate_figures.R)
##     plot_s4_weak(sum_mis)                 p_s4
##     plot_s5_correlated(sum_mis)           p_s5
##     plot_s6_pleiotropy(sum_mis)           p_s6
##     plot_s7_nonlinear(sum_mis)            p_s7
##
##   Convergence (convergence_diagnostics.R, regenerate_figures.R)
##     plot_trace(trace_tab)                 p_trace
##
## Every function applies theme_mrccc() and takes its notation from the
## lab_*() helpers in figure_style.R. Per-figure theme() additions are
## limited to layout necessities (axis-text angle for long categorical
## labels, legend placement for a long pathway legend) and are commented
## where they occur.
############################################################
library(dplyr)
library(tidyr)
library(forcats)
library(purrr)
library(tibble)
library(ggrepel)

############################################################
## REAL DATA: SHARED PREPARATION
############################################################

## HTML-formatted cell-type display names. The superscripts are why the
## real-data titles are rendered with element_markdown().
fmt_cell <- function(x) {
  dplyr::case_when(
    x == "BCells"         ~ "B cells",
    x == "CD4Cells"       ~ "CD4<sup>+</sup> T cells",
    x == "CD8Cells"       ~ "CD8<sup>+</sup> T cells",
    x == "NKCells"        ~ "NK cells",
    x == "MonocytesCells" ~ "monocytes",
    TRUE                  ~ x
  )
}

pair_label_for <- function(Cell1, Cell2) {
  paste0(fmt_cell(Cell1), " → ", fmt_cell(Cell2))
}

## Adds the two columns every real-data figure uses: the ligand-receptor
## label and the `discover` flag. `discover` follows the discovery rule of
## the analysis; both underlying rules remain available as disc_pip and
## disc_bfdr.
add_plot_columns <- function(out, discovery_rule = "pip") {
  out %>%
    mutate(
      LR       = paste(ligand_col_name, receptor_col_name, sep = "–"),
      discover = if (identical(discovery_rule, "bfdr")) disc_bfdr else disc_pip
    )
}

## Pathway colours, keyed by name. The colours themselves come from
## pal_pathways() in figure_style.R, which holds every palette and which
## errors rather than return NA when a level has no colour.
## The levels are fixed, not taken from the data of one pair, so that a
## pathway keeps the same colour in every figure even when a pair lacks
## some of the classes.
pathway_levels_all <- c(
  "Adhesion by ICAM", "Signaling by Chemokines", "Signaling by Cholesterol",
  "Signaling by GABA", "Signaling by Integrin", "Signaling by Interferon",
  "Signaling by Interleukin", "Signaling by Leptin",
  "Signaling by Prostaglandin", "Signaling by Semaphorin",
  "Signaling by Thromboxane", "Signaling by WNT")

pathway_palette <- function(out) {
  extra <- setdiff(sort(unique(out$pathway_name)), pathway_levels_all)
  levels_used <- c(pathway_levels_all, extra)
  setNames(pal_pathways(length(levels_used)), levels_used)
}

############################################################
## REAL DATA: p1 -- BUBBLE PLOT OF |beta_X| AND |beta_XZ|
##
## Each bubble: one (ligand-receptor, pathway) combination.
## Size   = absolute standardized posterior mean effect.
## Fill   = PIP (blue gradient: light = low, dark = high).
## Border = black for PIP > pip_thresh (discovered triplets).
############################################################
plot_bubble <- function(out, pair_label, pip_thresh) {
  lr_order <- out %>%
    arrange(gamma_mean) %>%
    pull(LR) %>%
    unique()

  df_long <- out %>%
    mutate(
      Pathway   = fct_reorder(pathway_name, gamma_mean, .fun = max),
      LR_factor = factor(LR, levels = lr_order),
      pip_group = if_else(gamma_mean > pip_thresh, "high", "low")
    ) %>%
    pivot_longer(
      cols      = c(Beta_X_mean, Beta_XZ_mean),
      names_to  = "effect",
      values_to = "beta"
    ) %>%
    mutate(
      effect   = recode(effect,
                        Beta_X_mean  = lab_strip_abs_beta_X_std(),
                        Beta_XZ_mean = lab_strip_abs_beta_XZ_std()),
      size_val = abs(beta)
    )

  ggplot(df_long, aes(x = Pathway, y = LR_factor)) +
    # Low-PIP: faint, grey border
    geom_point(
      data  = df_long %>% filter(pip_group == "low"),
      aes(size = size_val, fill = gamma_mean),
      shape = 21, color = "grey55", stroke = 0.3, alpha = 0.75
    ) +
    # High-PIP: vivid, black border
    geom_point(
      data  = df_long %>% filter(pip_group == "high"),
      aes(size = size_val, fill = gamma_mean),
      shape = 21, color = "black", stroke = 0.7, alpha = 0.95
    ) +
    facet_wrap(~effect, nrow = 1, labeller = label_parsed) +
    scale_fill_gradient(
      low = pip_gradient_low, high = pip_gradient_high, limits = c(0, 1),
      breaks = pip_breaks, labels = pip_labels,
      name = lab_pip(), guide = guide_pip_colourbar()
    ) +
    scale_size_continuous(
      range = c(1.5, 9),
      name  = "Standardized\neffect magnitude"
    ) +
    # One line per pathway, without the shared "Signaling by " prefix:
    # wrapped multi-line labels at an angle overlap their neighbours.
    scale_x_discrete(labels = short_pathway) +
    labs(
      title    = wrap_title(paste0("Communication effects: ", pair_label,
                                   "  (MR-CCC)"), 58, markdown = TRUE),
      subtitle = wrap_title(paste0(
        "Point size = |standardized posterior mean|;  color = PIP;  ",
        "black border = PIP > ", pip_thresh), 95),
      x = "Pathway (signaling class)",
      y = "Ligand–receptor"
    ) +
    theme_mrccc(markdown_title = TRUE) +
    # The y axis carries one row per triplet, about forty per cell-type
    # pair, so its text is set well below the theme size; at the theme's
    # 22 pt the rows collide into a solid block.
    # Legends at the bottom: the PIP colour bar is horizontal, and the two
    # pathway panels need the full width for eleven columns each.
    theme(legend.box         = "vertical",
          legend.position    = "bottom",
          axis.text.x        = element_text(angle = 50, hjust = 1, vjust = 1,
                                            size = 16),
          axis.text.y        = element_text(size = 14),
          panel.grid.major.y = element_blank())
}

############################################################
## REAL DATA: p2 -- RANKED LOLLIPOP OF ALL TRIPLETS BY PIP
##
## All (ligand-receptor-pathway) triplets ranked by PIP.
## Segments and points coloured by pathway.
## Black ring and PIP label on discovered triplets.
## Shaded region = discovery zone (PIP > pip_thresh).
############################################################
plot_pip_ranking <- function(out, pair_label, pip_thresh) {
  n_disc  <- sum(out$discover)
  n_total <- nrow(out)
  pathway_colors <- pathway_palette(out)

  df_lollipop <- out %>%
    mutate(LR_ranked = fct_reorder(LR, gamma_mean))

  ggplot(df_lollipop, aes(x = gamma_mean, y = LR_ranked)) +
    # Shaded discovery zone
    annotate("rect",
             xmin = pip_thresh, xmax = 1.02, ymin = -Inf, ymax = Inf,
             fill = pip_gradient_high, alpha = 0.04) +
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
      hjust = -0.30, size = 6, fontface = "bold", color = "grey15"
    ) +
    # Threshold line and annotation
    geom_vline(xintercept = pip_thresh,
               linetype = "dashed", color = "grey30", linewidth = 0.8) +
    annotate("text",
             x = pip_thresh - 0.02, y = 1.6,
             label = paste0("PIP = ", pip_thresh),
             hjust = 1, size = 6, color = "grey30", fontface = "italic") +
    scale_x_continuous(
      limits = c(0, 1.12),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0", "0.25", "0.50", "0.75", "1.00")
    ) +
    # Single-line legend entries: two-line entries in a two-column legend
    # overlap the entry below them.
    scale_color_manual(values = pathway_colors,
                       name = "Pathway (signaling class)",
                       labels = short_pathway) +
    labs(
      # Sentence case, and wrapped: this is the longest title in the
      # repository and at the figure width it would otherwise be clipped,
      # taking the receiver cell type with it.
      title    = wrap_title(paste0("PIP ranking across all triplets: ",
                                   pair_label, "  (MR-CCC)"), 58,
                            markdown = TRUE),
      subtitle = wrap_title(paste0(
        "Triplets ranked by PIP: ", n_total, " shown, ",
        n_disc,  " with PIP > ", pip_thresh,
        " (black ring);  shaded region = discovery zone"
      ), 78),
      x = "Posterior inclusion probability (PIP)",
      y = NULL
    ) +
    theme_mrccc(markdown_title = TRUE) +
    theme(
      # One row per triplet, about forty per cell-type pair, in a panel a
      # few inches tall: the y-axis text is set well below the theme size
      # so that adjacent labels do not touch.
      axis.text.y        = element_text(size = 13),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey88", linewidth = 0.4),
      # The pathway legend has twelve entries, so it is placed at the right
      # and given two columns; stacked in one column it is taller than the
      # figure and its last entry is cut off at the canvas edge.
      legend.position    = "right",
      legend.text        = element_text(size = 16),
      legend.title       = element_text(size = 20, face = "bold"),
      legend.key.size    = grid::unit(14, "pt"),
      legend.key.height  = grid::unit(24, "pt")
    ) +
    guides(color = guide_legend(ncol = 2, byrow = TRUE,
                                override.aes = list(linewidth = 2, size = 3)))
}

############################################################
## REAL DATA: p3 -- RECEPTOR-MODULATED LIGAND EFFECT CURVES
##           (discovery pairs only, PIP > pip_thresh)
##
## Each curve shows how the total standardized ligand effect
##   effect(z) = beta_X^(s) + beta_XZ^(s) * z
## varies with standardized receptor expression z = Z / SD(Z).
##
## Only pairs that exceed the PIP threshold are plotted.
## Only pathway panels that contain at least one such discovery
## are rendered (pathways with no discoveries are suppressed).
## Every shown curve is labelled at its rightmost data point.
############################################################
plot_effect_curves <- function(out, pair_label, pip_thresh) {
  n_disc <- sum(out$discover)

  # A cell-type pair with no triplet above the threshold has nothing to draw.
  # Returning an annotated empty panel, rather than proceeding, avoids an
  # error from unnesting an empty list-column, which would otherwise stop
  # the pipeline for that pair after its results had already been saved.
  if (sum(out$gamma_mean > pip_thresh, na.rm = TRUE) == 0L) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, size = 9,
                 family = FONT_FAMILY,
                 label = paste0("No triplet exceeds PIP = ", pip_thresh,
                                " for this cell-type pair")) +
        labs(title = wrap_title(paste0("Receptor-modulated ligand effects: ",
                                       pair_label, "  (MR-CCC)"), 58,
                                markdown = TRUE)) +
        # theme_void() for the empty panel, but the TITLE must still go
        # through the shared theme's markdown element: pair_label contains
        # HTML superscripts (CD4<sup>+</sup>), and a plain element_text
        # renders those tags literally in the figure.
        theme_void(base_size = 26, base_family = FONT_FAMILY) +
        theme(plot.title = ggtext::element_markdown(size = 30,
                                                    face = "bold", hjust = 0),
              plot.title.position = "plot",
              plot.margin = margin(t = 10, r = 18, b = 10, l = 10,
                                   unit = "pt"))
    )
  }

  # Identify pathways that contain at least one discovered pair.
  # Panels for pathways with zero discoveries are not drawn.
  sig_pathways <- out %>%
    filter(gamma_mean > pip_thresh) %>%
    pull(pathway_name) %>%
    unique()

  # Build long-format data for plotting:
  #   - restrict to discovered pairs only (PIP > pip_thresh)
  #   - restrict to panels that have at least one discovery
  #   - unnest Z_vec / Effect_vec into individual (Z, effect) rows
  #   - sort by Z within each triplet for correct line drawing
  #   - keep the central 98% of donors' receptor expression: normalised
  #     counts are right-skewed, and a handful of donors many SDs above the
  #     mean would otherwise stretch the axis and show the curve mostly
  #     where there are almost no data
  df_lines <- out %>%
    filter(pathway_name %in% sig_pathways,
           gamma_mean   > pip_thresh) %>%
    mutate(
      triplet_id = row_number(),
      paired     = purrr::pmap(
        list(Z_vec, Effect_vec, Effect_lo, Effect_hi),
        function(z, e, lo, hi) tibble::tibble(Z = z, effect = e,
                                              lo = lo, hi = hi))
    ) %>%
    select(triplet_id, pathway_name, LR, gamma_mean, paired) %>%
    tidyr::unnest(paired) %>%
    group_by(triplet_id) %>%
    filter(Z >= stats::quantile(Z, 0.01), Z <= stats::quantile(Z, 0.99)) %>%
    arrange(Z, .by_group = TRUE) %>%
    ungroup()

  # Label data: one label per curve, placed at the rightmost
  # observed receptor expression value for that triplet.
  # Because all displayed curves are discoveries, every curve
  # receives a label.
  # Labels sit just right of the curve ends and are spread vertically only,
  # so that several curves ending at the same receptor level (common within
  # one pathway) do not stack their labels on top of each other. The x axis
  # is extended on the right to make room for them.
  df_labels <- df_lines %>%
    group_by(pathway_name) %>%
    mutate(label_x = max(Z) + 0.03 * diff(range(Z))) %>%
    group_by(triplet_id, pathway_name, LR, gamma_mean, label_x) %>%
    slice_max(order_by = Z, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Subtitle built from the notation helpers so that it matches the axes.
  bX  <- lab_beta_X_std()[[1]]
  bXZ <- lab_beta_XZ_std()[[1]]
  # Two lines: as one line this subtitle is wider than the canvas and the
  # clause giving the number of curves shown is clipped. A plotmath
  # expression cannot be wrapped by strwrap, so the break is explicit via
  # atop(), which stacks two expressions.
  subtitle_expr <- bquote(
    atop(paste("Each curve: ", .(bX), " + ", .(bXZ),
               "  ·  (Z / SD(Z));  shaded band = pointwise 95% credible interval"),
         paste(.(n_disc), .(if (n_disc == 1) " triplet" else " triplets"),
               " with PIP > ", .(pip_thresh),
               " shown;  receptor range: central 98% of donors"))
  )

  ggplot() +
    # Pointwise 95% credible band, drawn first so the line sits on top.
    # Kept deliberately faint and uncoloured: the band conveys uncertainty,
    # the line conveys the estimate, and colouring both would make the PIP
    # scale harder to read.
    geom_ribbon(
      data = df_lines,
      aes(x = Z, ymin = lo, ymax = hi, group = triplet_id),
      fill = "grey55", alpha = 0.16
    ) +
    # Every displayed curve is a discovery (PIP > 0.5, stated in the
    # caption), so colour identifies the triplet instead of encoding PIP;
    # the label of each curve carries the same colour.
    geom_line(
      data      = df_lines,
      aes(x = Z, y = effect, group = triplet_id, color = LR),
      linewidth = 1.25, alpha = 0.95
    ) +
    # Label every curve at its rightmost point; leader segments connect
    # labels that had to be moved apart.
    geom_text_repel(
      data               = df_labels,
      aes(x = label_x, y = effect, label = LR, color = LR),
      size               = 5.5,
      family             = FONT_FAMILY,
      fontface           = "bold",
      hjust              = 0,
      direction          = "y",
      box.padding        = 0.35,
      force              = 2,
      min.segment.length = 0,
      max.overlaps       = Inf,
      show.legend        = FALSE
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.45))) +
    # One panel per pathway; only panels with discoveries are shown. Strip
    # names drop the shared "Signaling by " prefix so that they fit narrow
    # panels.
    facet_wrap(~pathway_name, scales = "free",
               labeller = as_labeller(short_pathway)) +
    scale_color_manual(
      values = setNames(
        rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
              "#56B4E9", "#920000", "#004949"),
            length.out = length(unique(df_lines$LR))),
        sort(unique(df_lines$LR))),
      guide = "none"
    ) +
    labs(
      title    = wrap_title(paste0("Receptor-modulated ligand effects: ",
                                   pair_label, "  (MR-CCC)"), 58,
                            markdown = TRUE),
      subtitle = subtitle_expr,
      x = "Standardized receptor expression  (Z / SD(Z))",
      # Shortened: the full parenthetical made this title taller than the
      # figure, and a rotated axis title that does not fit is clipped at
      # both ends. The units are given in the caption instead.
      y = "Standardized ligand effect"
    ) +
    theme_mrccc(markdown_title = TRUE)
}

############################################################
## SIMULATION S1-S3: LONG FORMAT
##
## Reshapes the wide results data frame (one row per replicate, as written
## to Results/simulation_S1_S3_raw.csv) into long format (one row per method
## per replicate) for plotting and summarisation. pip_thr re-derives the
## MR-CCC rejection here to keep the long format self-consistent.
############################################################
make_long_results <- function(df_raw, pip_thr = 0.5) {
  bind_rows(
    df_raw %>% transmute(
      scenario, n, gamma_true, beta_X_true, beta_XZ_true,
      method      = "OLS",
      score       = OLS_score,
      reject      = OLS_reject,
      beta_X_hat  = OLS_beta_X,
      beta_XZ_hat = OLS_beta_XZ
    ),
    df_raw %>% transmute(
      scenario, n, gamma_true, beta_X_true, beta_XZ_true,
      method      = "MVMR",
      score       = MVMR_score,
      reject      = MVMR_reject,
      beta_X_hat  = MVMR_beta_X,
      beta_XZ_hat = MVMR_beta_XZ
    ),
    df_raw %>% transmute(
      scenario, n, gamma_true, beta_X_true, beta_XZ_true,
      method      = "MR-BMA",
      score       = MRBMA_score,
      reject      = MRBMA_reject,
      beta_X_hat  = MRBMA_beta_X,
      beta_XZ_hat = MRBMA_beta_XZ  # NA: interaction not modeled
    ),
    df_raw %>% transmute(
      scenario, n, gamma_true, beta_X_true, beta_XZ_true,
      method      = "MR-CCC",
      score       = MR_score,
      reject      = as.integer(MR_score > pip_thr),
      beta_X_hat  = MR_beta_X,
      beta_XZ_hat = MR_beta_XZ
    )
  ) %>%
    mutate(method = factor(method, levels = method_levels))
}

############################################################
## SIMULATION S1-S3: p_score -- COMMUNICATION SCORE
##
## Score definition per method:
##   OLS   : 1 - p_F  (joint F-test, H0: beta_X = beta_XZ = 0)
##   MVMR  : 1 - p_t  (t-test, H0: beta_X = 0)
##   MR-BMA: MIP_X    (marginal inclusion probability for X)
##   MR-CCC: PIP      (posterior inclusion probability P(gamma=1|data))
## Red dashed line: true communication state (0 for S1, 1 for S2/S3)
############################################################
plot_sim_score <- function(sim_long) {
  gamma_df <- sim_long %>%
    distinct(scenario) %>%
    mutate(gamma_true = case_when(scenario == "S1" ~ 0, TRUE ~ 1))

  ggplot(sim_long, aes(x = method, y = score, fill = method)) +
    geom_boxplot(outlier.alpha = 0.35, alpha = 0.85) +
    scale_fill_manual(values = pal_methods, guide = "none") +
    geom_hline(
      data      = gamma_df,
      aes(yintercept = gamma_true),
      color     = "red",
      linetype  = "dashed",
      linewidth = 0.9
    ) +
    facet_grid(scenario ~ n, scales = "free_y") +
    labs(
      title    = "Communication scores across methods and scenarios",
      subtitle = "Rows = scenario; columns = sample size",
      x        = "Method",
      y        = "Communication score"
    ) +
    theme_mrccc() +
    # Method names on a four-level categorical axis; angled to avoid overlap.
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 20))
}

############################################################
## SIMULATION S1-S3: p_betaX -- beta_X ESTIMATES
##
## OLS / MVMR / MR-CCC: point estimate of main causal effect
## MR-BMA: model-averaged causal effect (MACE) for X
## Red dashed line: true beta_X value
############################################################
plot_sim_betaX <- function(sim_long) {
  betaX_df <- sim_long %>% distinct(scenario, n, beta_X_true)
  bX <- lab_beta_X()[[1]]

  ggplot(sim_long, aes(x = method, y = beta_X_hat, fill = method)) +
    geom_boxplot(outlier.alpha = 0.35, alpha = 0.85) +
    scale_fill_manual(values = pal_methods, guide = "none") +
    geom_hline(
      data      = betaX_df,
      aes(yintercept = beta_X_true),
      colour    = "red",
      linetype  = "dashed",
      linewidth = 0.9
    ) +
    facet_grid(scenario ~ n, scales = "free_y") +
    labs(
      title    = bquote(paste("Estimated ", .(bX),
                              " across methods and scenarios")),
      subtitle = "Red dashed line = true value",
      x        = "Method",
      y        = lab_beta_X()
    ) +
    theme_mrccc() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 20))
}

############################################################
## SIMULATION S1-S3: p_betaXZ -- beta_XZ ESTIMATES
##
## Only OLS and MR-CCC estimate the interaction term.
## MR-BMA and MVMR are excluded (neither models beta_XZ).
## Red dashed line: true beta_XZ value.
############################################################
plot_sim_betaXZ <- function(sim_long) {
  sim_long_xz <- sim_long %>%
    filter(method %in% c("OLS", "MR-CCC")) %>%
    mutate(method = factor(method, levels = c("OLS", "MR-CCC")))
  betaXZ_df_xz <- sim_long_xz %>%
    distinct(scenario, n, beta_XZ_true)
  bXZ <- lab_beta_XZ()[[1]]

  ggplot(sim_long_xz, aes(x = method, y = beta_XZ_hat, fill = method)) +
    geom_boxplot(outlier.alpha = 0.35, alpha = 0.85) +
    scale_fill_manual(values = pal_methods, guide = "none") +
    geom_hline(
      data      = betaXZ_df_xz,
      aes(yintercept = beta_XZ_true),
      colour    = "red",
      linetype  = "dashed",
      linewidth = 0.9
    ) +
    facet_grid(scenario ~ n, scales = "free_y") +
    labs(
      title    = bquote(paste("Estimated ", .(bXZ),
                              " across methods and scenarios")),
      subtitle = "Red dashed line = true value  (MR-BMA and MVMR excluded: interaction not modeled)",
      x        = "Method",
      y        = lab_beta_XZ()
    ) +
    theme_mrccc() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 20))
}

############################################################
## SIMULATION S4-S7: MISSPECIFICATION FIGURES
##
## One figure per scenario, each plotting the quantity that scenario is
## designed to interrogate against the departure being varied.
##
## The NULL arm and the SIGNAL arm are shown together because they answer
## different questions from the same run: under the null the rejection rate
## is the false positive rate, so calibration is the issue; under signal it
## is power. A method can hold calibration and lose all power, or keep power
## by over-rejecting, and only the two panels together distinguish those.
############################################################

## Facet labels shared by S4-S7. The three arms mirror scenarios S1-S3 and
## are displayed in that order; mis_arm_order() imposes it, since the
## default alphabetical ordering would place the main-effect arm first.
## Unused levels are dropped, so data produced with fewer arms still plot.
## These label the ROWS of facet_grid, whose strips are rotated, so the
## available space is the panel HEIGHT rather than its width. The longer
## forms ("Main-effect arm (power; no true interaction)") do not fit and are
## clipped to an unreadable fragment, so each is stated in the fewest words
## that still say which quantity the row reports. The full definitions of
## the three arms belong in the caption.
mis_arm_labels <- c(
  null   = "Null\n(size)",
  signal = "Signal\n(power)",
  main   = "Main\n(power)")
mis_arm_order  <- function(df) {
  df$arm <- factor(df$arm, levels = c("null", "signal", "main"))
  df[!is.na(df$arm), , drop = FALSE]
}
mis_n_label    <- function(v) paste0("n = ", v)

## Methods are drawn in the display order of method_levels. Without this the
## legend of each S4-S7 figure falls back to character order (MR-BMA,
## MR-CCC, MVMR, OLS) while the S1-S3 figures and the boxplots use the
## display order, so the same supplementary section carries two different
## method orderings.
mis_method_order <- function(df) {
  df$method <- factor(df$method, levels = method_levels)
  df[!is.na(df$method), , drop = FALSE]
}

## Common layers for S4 and S5. `xvar` is the departure being varied.
## Plus or minus one Monte Carlo standard error on each rejection rate.
## Drawn rather than left to the caption, because the whole reading of these
## panels is whether a rate departs from the nominal 0.05, and that cannot
## be judged without the precision of the estimate. Returns NULL -- a no-op
## when added to a plot -- if the MCSE column is absent from the summary.
mis_mcse_layer <- function(df) {
  if (!"reject_mcse" %in% names(df)) return(NULL)
  geom_errorbar(aes(ymin = pmax(reject_rate - reject_mcse, 0),
                    ymax = pmin(reject_rate + reject_mcse, 1)),
                width = 0, linewidth = 0.7, alpha = 0.75,
                show.legend = FALSE)
}

mis_layers <- function(df, xvar, xlab, logx = FALSE) {
  df <- mis_method_order(mis_arm_order(df))
  p <- ggplot(df, aes(x = .data[[xvar]], y = reject_rate,
                      colour = method, shape = method, group = method)) +
    # 0.05 marks nominal size; it is meaningful in the null panel only.
    geom_hline(yintercept = 0.05, linetype = "dashed",
               colour = "grey35", linewidth = 0.8) +
    mis_mcse_layer(df) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 4.5) +
    facet_grid(arm ~ n, labeller = labeller(arm = mis_arm_labels,
                                            n   = mis_n_label)) +
    scale_colour_manual(values = pal_methods, name = "Method") +
    scale_shape_manual(values  = shp_methods, name = "Method") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = xlab, y = "Rejection rate") +
    theme_mrccc()
  if (logx) p <- p + scale_x_log10(breaks = c(2, 3, 10, 30, 100)) else p
}

## S4: weak instruments. Instrument strength is swept from the regime
## observed in the real analysis (F ~ 1.5) to that of scenarios S1-S3
## (F ~ 70), so both sit on one continuous axis. Plotted against the
## REALISED median F rather than the target, since the realised value is
## what a reader can compare with the F statistics reported for the data.
plot_s4_weak <- function(sum_mis) {
  mis_layers(sum_mis %>% filter(scenario == "S4"),
             "F_obs_median", "Median first-stage F (log scale)",
             logx = TRUE) +
    labs(title = "S4: behavior under weak instruments",
         subtitle = wrap_title(paste(
           "Dashed line = nominal 0.05 size.",
           "Left of the plot is the regime observed in the data."), 88))
}

## S5: correlated instruments.
plot_s5_correlated <- function(sum_mis) {
  mis_layers(sum_mis %>% filter(scenario == "S5"),
             "iv_corr", "AR(1) correlation between instruments") +
    scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9)) +
    labs(title = "S5: behavior under correlated instruments",
         subtitle = wrap_title(paste(
           "The largest correlation shown sets the clumping threshold",
           "used on the real data."), 88))
}

## S6: invalid / pleiotropic instruments. Faceted additionally by the
## FRACTION of invalid instruments, so that the breakdown point can be read
## in two directions: how many instruments are invalid, and how large their
## direct effects are.
## S6 differs from S4, S5 and S7 only in carrying a second column variable,
## the fraction of invalid instruments, so it adds that facet to the shared
## mis_layers() rather than restating the whole plot. Adding a facet_grid
## replaces the one mis_layers() set, which is what makes this work; writing
## the layers out again would leave two copies to keep in step.
plot_s6_pleiotropy <- function(sum_mis) {
  mis_layers(sum_mis %>% filter(scenario == "S6"),
             "pleio_delta", "Magnitude of direct effect") +
    facet_grid(arm ~ n + pleio_frac, labeller = labeller(
      arm        = mis_arm_labels,
      n          = mis_n_label,
      pleio_frac = function(v) paste0("Invalid = ", v))) +
    scale_x_continuous(breaks = c(0.1, 0.3)) +
    theme(panel.spacing.x = grid::unit(1.6, "lines")) +
    labs(title = "S6: behavior under invalid (pleiotropic) instruments",
         subtitle = wrap_title(paste(
           "Exclusion restriction violated by a direct instrument effect",
           "on the outcome."), 88))
}

## S7: nonlinear ligand effect.
## S7 sweeps the STRENGTH of the nonlinearity rather than fixing one value,
## so the figure is a dose-response in nl_coef on the same line-and-point
## form as S4-S6, not a bar chart over methods. nl_coef is on the same scale
## as beta_X = 0.3, so the right-hand end of the axis is curvature as large
## as the linear effect itself.
plot_s7_nonlinear <- function(sum_mis) {
  mis_layers(sum_mis %>% filter(scenario == "S7"),
             "nl_coef", "Strength of the quadratic term") +
    scale_x_continuous(breaks = c(0.05, 0.1, 0.2, 0.3)) +
    labs(title = "S7: behavior under a nonlinear ligand effect",
         subtitle = wrap_title(paste(
           "Outcome generated with a quadratic ligand effect;",
           "working model is linear plus interaction.",
           "In this scenario the null arm keeps the quadratic effect."), 88))
}

############################################################
## MISSPECIFICATION: ESTIMATE-LEVEL BOXPLOTS
##
## The line plots above show operating characteristics along each swept
## departure. These boxplots complement them with the per-replicate
## ESTIMATES, in the same visual language as the S1-S3 figures, so that the
## two halves of the simulation study read as one: where a score box sits
## (near 0 in the null arm, near 1 otherwise), whether beta_X is attenuated,
## and -- in the main-effect arm, where the true interaction is zero --
## whether the misspecification manufactures a spurious beta_XZ.
##
## Each scenario is shown at TWO representative settings of its departure
## rather than at every grid point, because a box per setting would be
## unreadable and the trends are already in the line plots: S4 at target
## F = 70 and 1, S5 at rho = 0 and 0.9, S6 at (fraction, effect) = (0.2,
## 0.1) and (0.4, 0.3), and S7 at the smallest and largest nl_coef. Only
## the n = 500 replicates are plotted (n_keep below).
############################################################
## Representative settings for the estimate boxplots. Each scenario
## contributes TWO columns wherever two informative settings exist, ordered
## from the milder to the more severe departure, so that the boxplot carries
## a within-scenario contrast rather than a single slice:
##
##   S4  target F = 70 and target F = 1, realised medians of about 85 and 2.
##       The first is effectively no misspecification -- the strong-instrument
##       regime of the correctly specified study -- and the second is the grid
##       point nearest the first-stage strength observed in the real analysis,
##       so the pair brackets the range the paper reasons about. Strips show
##       the realised values.
##   S5  rho = 0 and rho = 0.9. rho = 0 is exactly the S1-S3 design with
##       independent instruments, giving a second reference column.
##   S6  fraction 0.2 with effect 0.1, and fraction 0.4 with effect 0.3.
##       The grid contains no zero-pleiotropy cell, so the contrast is mild
##       against severe rather than absent against present.
##   S7  nl_coef 0.05 and 0.3, the mild and severe ends of the swept
##       curvature. Pooling the four strengths into one box would make S7
##       the only column that mixes settings, and would hide the very
##       contrast the scenario exists to show.
##
## The `panel` column separates settings that share a scenario label.
mis_representative <- function(mis_raw) {
  dplyr::bind_rows(
    mis_raw %>% filter(scenario == "S4", target_F == 70) %>%
      mutate(panel = "S4_strong"),
    mis_raw %>% filter(scenario == "S4", target_F == 1) %>%
      mutate(panel = "S4_weak"),
    mis_raw %>% filter(scenario == "S5", iv_corr == 0) %>%
      mutate(panel = "S5_indep"),
    mis_raw %>% filter(scenario == "S5", iv_corr == 0.9) %>%
      mutate(panel = "S5_corr"),
    mis_raw %>% filter(scenario == "S6", pleio_frac == 0.2,
                       pleio_delta == 0.1) %>%
      mutate(panel = "S6_mild"),
    mis_raw %>% filter(scenario == "S6", pleio_frac == 0.4,
                       pleio_delta == 0.3) %>%
      mutate(panel = "S6_severe"),
    # The extremes are taken over the S7 rows only: S4-S6 carry
    # nl_coef = 0, so a minimum over all rows would select nothing.
    mis_raw %>% filter(scenario == "S7") %>%
      filter(nl_coef == min(nl_coef)) %>%
      mutate(panel = "S7_mild"),
    mis_raw %>% filter(scenario == "S7") %>%
      filter(nl_coef == max(nl_coef)) %>%
      mutate(panel = "S7_severe")
  ) %>%
    mutate(panel = factor(panel,
                          levels = c("S4_strong", "S4_weak",
                                     "S5_indep", "S5_corr",
                                     "S6_mild", "S6_severe",
                                     "S7_mild", "S7_severe")))
}

## Column strips of the boxplots. Kept short because facet_grid gives each
## strip only the width of its own column, and a longer label is clipped:
## the scenario prefix and the parameter value are exactly what a reader
## needs to tell the columns apart, so neither can be the part that is lost.
## The parameter being varied is named in the axis titles and the caption.
mis_scenario_labels <- c(
  S4_strong = "S4\nF ≈ 85",
  S4_weak   = "S4\nF ≈ 2",
  S5_indep  = "S5\nρ = 0",
  S5_corr   = "S5\nρ = 0.9",
  S6_mild   = "S6\n20%, 0.1",
  S6_severe = "S6\n40%, 0.3",
  S7_mild   = "S7\nmild",
  S7_severe = "S7\nsevere")

## Long format for the estimate boxplots: one row per replicate x method,
## with the score and both coefficient estimates as columns.
mis_estimates_long <- function(mis_raw) {
  mis_representative(mis_raw) %>%
    tidyr::pivot_longer(
      cols = dplyr::matches("^(OLS|MVMR|MRBMA|MR)_(score|beta_X|beta_XZ)$"),
      names_to  = c("method", ".value"),
      names_pattern = "^(OLS|MVMR|MRBMA|MR)_(.*)$") %>%
    mutate(method = dplyr::recode(method, OLS = "OLS", MVMR = "MVMR",
                                  MRBMA = "MR-BMA", MR = "MR-CCC"),
           method = factor(method, levels = method_levels))
}

## One boxplot figure per quantity. Truth reference lines are drawn per arm
## from the effect sizes shared with the S1-S3 study: the score should sit
## near 0 in the null arm and near 1 otherwise; beta_X is 0.3 except in the
## null arm; beta_XZ is 0.3 only in the signal arm, so in the main-effect
## arm any box away from zero is a spurious interaction.
plot_mis_estimate_box <- function(mis_raw, quantity = c("score", "beta_X",
                                                        "beta_XZ"),
                                  n_keep = 500,
                                  beta_X_true = 0.3, beta_XZ_true = 0.3) {
  quantity <- match.arg(quantity)
  df <- mis_estimates_long(mis_raw) %>% filter(n == n_keep)
  df$value <- df[[quantity]]
  # The interaction is modelled by OLS and MR-CCC only; the other two
  # methods have no beta_XZ and are dropped from that panel.
  if (quantity == "beta_XZ") {
    df <- df %>% filter(method %in% c("OLS", "MR-CCC"), !is.na(value))
  }
  df <- mis_arm_order(df)
  # Outliers are not drawn, so they must not set the axis either: the
  # limits cover the central 99% of the plotted estimates.
  ylims <- stats::quantile(df$value, c(0.005, 0.995), na.rm = TRUE)
  if (quantity == "score") ylims <- c(0, 1)

  truth <- data.frame(arm = factor(c("null", "signal", "main"),
                                   levels = c("null", "signal", "main")))
  truth$yint <- switch(quantity,
    score   = c(0, 1, 1),
    beta_X  = c(0, beta_X_true, beta_X_true),
    beta_XZ = c(0, beta_XZ_true, 0))

  y_lab <- switch(quantity,
    score   = "Communication score",
    beta_X  = lab_beta_X(),
    beta_XZ = lab_beta_XZ())
  # The notation helpers in figure_style.R are the only place these symbols
  # are spelled out, so the title reuses them rather than restating the
  # plotmath; a second spelling is how two figures end up disagreeing.
  ttl <- switch(quantity,
    score   = "Communication scores under misspecification",
    beta_X  = bquote(paste("Estimated ", .(lab_beta_X()[[1]]),
                           " under misspecification")),
    beta_XZ = bquote(paste("Estimated ", .(lab_beta_XZ()[[1]]),
                           " under misspecification")))

  ggplot(df, aes(x = method, y = value, fill = method)) +
    geom_hline(data = truth, aes(yintercept = yint),
               linetype = "dashed", colour = "red3", linewidth = 0.8) +
    # Shape is carried as well as fill so that the methods remain
    # separable in greyscale, as in every other method-coloured figure:
    # four Okabe-Ito fills become four similar greys when printed in black
    # and white, so fill alone would not separate them.
    geom_boxplot(aes(shape = method), width = 0.65, outlier.size = 1.6,
                 outlier.shape = NA, alpha = 0.9, linewidth = 0.35) +
    coord_cartesian(ylim = ylims) +
    geom_point(aes(shape = method), stat = "summary", fun = median,
               size = 2.6, colour = "black", show.legend = FALSE) +
    facet_grid(arm ~ panel, labeller = labeller(
      arm = mis_arm_labels, panel = mis_scenario_labels)) +
    scale_fill_manual(values = pal_methods, name = "Method") +
    scale_shape_manual(values = shp_methods, name = "Method") +
    labs(title = ttl,
         subtitle = wrap_title(paste0(
           "Two settings per scenario, milder then more severe; n = ",
           n_keep, ". Red dashed line = target value."), 88),
         x = "Method", y = y_lab) +
    # These figures carry eight scenario columns and are therefore the
    # widest in the repository, so they are reduced further in print; the
    # scale keeps their text in the same range as the rest.
    theme_mrccc() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 18)) +
    # The x axis already names the methods and is titled "Method", so a
    # legend would only repeat it, in greyscale as well as in colour.
    guides(fill = "none", shape = "none")
}

############################################################
## CONVERGENCE: p_trace -- MCMC TRACE PLOTS
##
## trace_tab holds one row per displayed draw with columns
##   triplet    facet column label, ligand on line 1 and receptor on line 2
##   parameter  "log-likelihood", "Beta_X", "Beta_XZ" or
##              "gamma (running mean)"
##   chain      factor, chain index
##   iter       retained iteration
##   value      sampled value (running mean for the inclusion indicator)
## and carries the run settings as attributes n_chains, n_iter and burn_in,
## which convergence_diagnostics.R attaches before saving. The subtitle is
## built from those attributes; if they are absent, only the chain count,
## taken from the data, is stated.
##
## Facet rows are ordered by trace_parameter_levels() so that the
## log-likelihood sits above the inclusion indicator whichever subset of
## series was recorded.
############################################################
plot_trace <- function(trace_tab,
                       n_chains = attr(trace_tab, "n_chains"),
                       n_iter   = attr(trace_tab, "n_iter"),
                       burn_in  = attr(trace_tab, "burn_in")) {
  if (is.null(n_chains)) n_chains <- nlevels(factor(trace_tab$chain))

  # Order the facet rows. Levels present in the data are taken in the
  # canonical order; anything unrecognised is appended rather than dropped.
  lv_known   <- trace_parameter_levels()
  lv_present <- unique(as.character(trace_tab$parameter))
  trace_tab$parameter <- factor(
    trace_tab$parameter,
    levels = c(intersect(lv_known, lv_present), setdiff(lv_present, lv_known)))

  has_gamma <- "gamma (running mean)" %in% lv_present
  gamma_note <- if (has_gamma) {
    paste0(" The inclusion indicator is shown as a running mean, since it ",
           "takes only the values 0 and 1.")
  } else ""

  subtitle <- if (!is.null(n_iter) && !is.null(burn_in)) {
    paste0(
      n_chains, " chains, ", format(n_iter, big.mark = ","),
      " iterations, burn-in ", format(burn_in, big.mark = ","),
      "; thinned for display only.", gamma_note)
  } else {
    paste0(
      n_chains, " chains from dispersed starts; thinned for display only.",
      gamma_note)
  }

  # Short, because a rotated axis title longer than the figure is tall is
  # clipped at both ends; the running-mean qualification is in the subtitle.
  y_lab <- "Sampled value"

  ggplot(trace_tab, aes(x = iter, y = value, colour = chain)) +
    geom_line(linewidth = 0.35, alpha = 0.75) +
    # facet_wrap rather than facet_grid so that every panel has its own
    # y-axis: the triplets' log-likelihoods differ by hundreds of units, and a
    # shared axis flattens each trace into a line.
    facet_wrap(vars(parameter, triplet), scales = "free_y",
               ncol = max(length(unique(trace_tab$triplet)), 1L),
               labeller = labeller(
                 parameter = as_labeller(lab_strip_trace_parameters(),
                                         default = label_parsed))) +
    # pal_chains is named, and the chain factor has bare integers as levels,
    # so the values are indexed by position rather than by name. Taking only
    # as many colours as there are chains keeps the scale valid if the chain
    # count is ever changed, instead of erroring on a length mismatch.
    scale_colour_manual(values = unname(pal_chains)[seq_len(
                          max(nlevels(factor(trace_tab$chain)), 1L))],
                        labels = function(x) paste("Chain", x),
                        name   = "Chain") +
    # Three y breaks per panel: the panels are short, and more labels run
    # into each other.
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.6,
                                                     alpha = 1))) +
    # At most THREE breaks, in thousands. With four, the last label of one
    # panel and the first of the next meet at the panel boundary and render
    # as a single run of digits.
    scale_x_continuous(limits = c(0, NA),
                       breaks = breaks_at_most(3L),
                       labels = label_thousands,
                       expand = expansion(mult = 0.06)) +
    labs(
      title    = "MCMC trace plots across dispersed chains",
      subtitle = wrap_title(subtitle, 88),
      x = "Retained iteration",
      y = y_lab
    ) +
    theme_mrccc() +
    # Extra space between columns so that the last iteration label of one
    # panel and the first of the next do not run together.
    theme(panel.spacing.x = grid::unit(2.2, "lines"))
}
