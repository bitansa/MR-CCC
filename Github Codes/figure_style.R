############################################################
## figure_style.R
##
## SINGLE SOURCE OF TRUTH FOR FIGURE APPEARANCE
##
## Every figure produced by this repository -- the real-data figures
## (real_data_analysis.R), the simulation figures (simulation_mrccc.R,
## simulation_misspecification.R) and the convergence trace plot
## (convergence_diagnostics.R) -- takes its theme, colours, point shapes,
## mathematical notation and export settings from this file. No other
## script defines a theme, a palette, a plotmath label for a model quantity,
## or a call to ggsave(); they source this file and use what it provides.
## The plot-building functions themselves live in figure_functions.R, which
## this file sources, so that there is exactly one definition of each plot.
##
## USAGE (from the project root):
##   source("Github Codes/figure_style.R")
##
## PROVIDES:
##   theme_mrccc(markdown_title = FALSE)
##       The one theme. markdown_title = TRUE renders the plot title with
##       ggtext::element_markdown(), which the real-data figures need because
##       their titles contain HTML superscripts (CD4<sup>+</sup>).
##   pal_okabe_ito, pal_methods, pal_chains, shp_methods, method_levels,
##   pal_pathways(n), pip_gradient_low, pip_gradient_high
##       Colour-vision-deficiency-safe palette and its named sub-palettes,
##       the categorical pathway palette, and the two ends of the PIP
##       gradient used by every figure that colours by PIP.
##   wrap_title(x, width), wrap_label(width)
##       Text wrapping for titles, subtitles and facet strips. Every string
##       long enough to leave the canvas passes through one of these.
##   guide_pip_colourbar()
##       The one continuous-colour guide, sized so that its tick labels do
##       not collide at the theme's font size.
##   lab_beta_X_std(), lab_beta_XZ_std(), lab_beta_X(), lab_beta_XZ(),
##   lab_pip()
##       Plotmath expressions for every model quantity that appears on an
##       axis, in a title or in a legend, so that notation cannot drift
##       between figures.
##   lab_strip_abs_beta_X_std(), lab_strip_abs_beta_XZ_std(),
##   lab_strip_trace_parameters()
##       The same notation as parseable strings, for facet strips.
##   breaks_at_most(n_max), label_thousands
##       Axis helpers for long MCMC iteration axes.
##   save_fig(plot, name, width, height)
##       The only export path: Plots/<name>.pdf via ggsave(device = cairo_pdf).
##
## NOTATION RULE
##   Real-data figures report STANDARDIZED posterior means and therefore use
##   the superscripted forms hat(beta)[X]^{(s)} and hat(beta)[XZ]^{(s)}.
##   Simulation figures report estimates on the RAW scale, against known
##   truths, and use the un-superscripted forms hat(beta)[X] and
##   hat(beta)[XZ]. Any title or label that refers to an ESTIMATE carries the
##   hat, exactly as the axis it describes does.
##
## CAPITALISATION RULE: SENTENCE CASE
##   Every axis title, plot title, subtitle, legend title and facet label
##   begins with a capital letter; the remaining words are lower case except
##   proper nouns, gene symbols and acronyms (PIP, MCMC, MR-CCC, MVMR,
##   MR-BMA, OLS, SD, FDR, LD, AR(1)). Examples:
##     "Rejection rate"                      not "Rejection Rate"
##     "Median first-stage F (log scale)"    not "Median First-Stage F (Log Scale)"
##     "Null arm (false positive rate)"      not "Null Arm (False Positive Rate)"
##     "Signal arm (power)"                  not "Signal Arm (Power)"
##     "Invalid fraction = 0.2"              not "Invalid Fraction = 0.2"
##
## FONT SIZES AND CANVAS WIDTH
##   Figures are exported on a large canvas and reproduced in the article at
##   a single-column width of about seven inches, so the effective reduction
##   is 7/width and the printed size of any text is (source size) * 7/width.
##   The sizes below put axis text at 11 pt and axis titles at 14 pt for a
##   14-inch figure, which is the width used by almost every figure here.
##
##   A figure must therefore not be made wider without reducing its font
##   sizes to match: at 26 inches the same theme prints axis text at 5.9 pt,
##   below the 7 pt floor for a Bioinformatics figure. No figure in this
##   repository exceeds 18 inches, and the widest ones carry
##   theme_mrccc(scale = ...) to compensate.
##
## TEXT MUST FIT THE CANVAS
##   A title, subtitle, axis title or facet strip that is wider than the
##   figure is silently clipped by the graphics device: the text is drawn
##   and the part beyond the edge is simply not shown, with no warning. The
##   information lost is usually the end of the string, which is where
##   thresholds, counts and parameter values sit. Every such string is
##   therefore wrapped with wrap_title() or wrap_label(), and the theme sets
##   a plot margin so that wrapped text has somewhere to go.
##
## EXPORT
##   All figures are vector PDF written through cairo_pdf, which embeds the
##   fonts. A vector figure has no fixed resolution and satisfies any
##   dots-per-inch requirement exactly. save_fig() checks that the running R
##   build has cairo and stops with an actionable message if it does not,
##   rather than failing inside ggsave.
############################################################
library(ggplot2)
library(scales)

if (!requireNamespace("ggtext", quietly = TRUE)) {
  stop("Package 'ggtext' is required for the real-data figure titles. ",
       "Install with install.packages(\"ggtext\").")
}

############################################################
## THEME
############################################################
## One family for every glyph. Without this, cairo falls back per character
## and a single figure can mix three typefaces: the body text in one, and
## the arrow, en-dash and Greek letters of a title in another. Helvetica is
## present on macOS and Windows; on a Linux build without it, cairo falls
## back to its default sans and the figure is still internally consistent.
FONT_FAMILY <- "Helvetica"

## scale multiplies every font size. A figure wider than the 14-inch
## reference is reduced further in print, so it passes scale < 1 to keep the
## printed point size in range; see the canvas-width note in the header.
theme_mrccc <- function(markdown_title = FALSE, scale = 1) {
  sz <- function(x) x * scale
  title_element <- if (isTRUE(markdown_title)) {
    ggtext::element_markdown(size = sz(30), face = "bold", hjust = 0)
  } else {
    element_text(size = sz(30), face = "bold", hjust = 0)
  }
  theme_bw(base_size = sz(26), base_family = FONT_FAMILY) +
    theme(
      plot.title       = title_element,
      plot.subtitle    = element_text(size = sz(24), colour = "grey40", hjust = 0),
      axis.title       = element_text(size = sz(28), face = "bold"),
      axis.text        = element_text(size = sz(22), colour = "black"),
      strip.text       = element_text(size = sz(24), face = "bold"),
      strip.background = element_rect(fill = "grey93", colour = NA),
      legend.title     = element_text(size = sz(26), face = "bold"),
      legend.text      = element_text(size = sz(22)),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      # Titles are left-aligned to the whole plot rather than to the panel,
      # so a long title starts at the far left and has the full width to
      # run in, and the margin keeps wrapped lines off the canvas edge.
      plot.title.position = "plot",
      plot.margin      = margin(t = 10, r = 18, b = 10, l = 10, unit = "pt")
    )
}

############################################################
## TEXT WRAPPING
##
## Clipping is silent, so every string that can exceed the canvas is
## wrapped here rather than trusted to fit. The widths are in characters
## and are chosen for the theme's sizes at the figure widths actually used:
## a 14-inch canvas holds roughly 70 characters of 30 pt title text and 95
## of 24 pt subtitle text.
############################################################
## markdown = TRUE joins the wrapped lines with "<br>" instead of a newline.
## This is required for any string rendered by ggtext::element_markdown(),
## which is how the real-data titles are drawn so that CD4<sup>+</sup>
## renders as a superscript: in markdown a single newline is not a line
## break, so a "\n" is silently collapsed and the title is drawn on one line
## and clipped. The HTML tags themselves are counted by strwrap but occupy
## almost no width when rendered, so a markdown string wraps slightly
## earlier than its character count suggests, which errs the safe way.
wrap_title <- function(x, width = 70, markdown = FALSE) {
  paste(strwrap(x, width = width),
        collapse = if (isTRUE(markdown)) "<br>" else "\n")
}

## For facet strips and discrete axis labels, as a labeller-friendly
## function: wrap_label(24) returns a function suitable for
## scale_*_discrete(labels = ) or as_labeller().
wrap_label <- function(width = 24) {
  function(x) vapply(as.character(x),
                     function(s) paste(strwrap(s, width = width),
                                       collapse = "\n"),
                     character(1), USE.NAMES = FALSE)
}

############################################################
## PALETTES AND SHAPES
##
## Okabe-Ito (Okabe & Ito, 2008) is distinguishable under the common forms
## of colour vision deficiency. Methods are additionally separated by point
## shape so that the simulation panels remain readable in greyscale.
############################################################
pal_okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#000000")

## Display order of the four methods on categorical axes.
method_levels <- c("OLS", "MVMR", "MR-BMA", "MR-CCC")

pal_methods <- c("MR-CCC" = "#0072B2",
                 "MVMR"   = "#D55E00",
                 "MR-BMA" = "#009E73",
                 "OLS"    = "#CC79A7")

shp_methods <- c("MR-CCC" = 16, "MVMR" = 17, "MR-BMA" = 15, "OLS" = 18)

pal_chains <- c("Chain 1" = "#0072B2",
                "Chain 2" = "#D55E00",
                "Chain 3" = "#009E73",
                "Chain 4" = "#CC79A7")

## Categorical palette for Reactome pathways. The data carry twelve distinct
## pathways on most cell-type axes, so a palette of fewer than twelve leaves
## one pathway with no colour: ggplot assigns NA, which draws its points and
## its legend key invisible, and the pathway disappears from the figure
## without any warning. pal_pathways() therefore errors rather than recycle
## or return NA if asked for more colours than it has.
##
## The twelve are the eight Okabe-Ito colours followed by four further hues
## chosen to stay separable under deuteranopia and protanopia: a dark blue,
## a light blue, a dark red and a mid grey. Categorical pathway identity is
## secondary information in these figures -- the quantity being read is the
## inclusion probability -- so twelve categories is at the limit of what
## colour alone can carry, and the figures that use it also order the points
## by probability so that a reader never depends on the colour alone.
pal_pathways_full <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                       "#0072B2", "#D55E00", "#CC79A7", "#000000",
                       "#004949", "#92DBF2", "#920000", "#888888")

pal_pathways <- function(n) {
  if (n > length(pal_pathways_full)) {
    stop("pal_pathways(): asked for ", n, " colours but only ",
         length(pal_pathways_full), " are defined. Add colours rather than ",
         "recycling: a pathway with no colour is drawn invisible.")
  }
  pal_pathways_full[seq_len(n)]
}

## The two ends of the PIP gradient. Defined once here because three figures
## colour or fill by PIP and a reader compares them side by side.
pip_gradient_low  <- "#deebf7"
pip_gradient_high <- "#08306b"

## The one continuous-colour guide. At the theme's legend text size the
## default colourbar is too narrow for its own tick labels, which then
## overlap into a single run of digits; this sets a width that fits five
## labels and reduces the label size for continuous guides only.
guide_pip_colourbar <- function() {
  # Wide enough for five labels at 16 pt to sit apart; on a much narrower
  # bar adjacent labels overlap and become unreadable.
  guide_colourbar(barwidth = grid::unit(260, "pt"),
                  barheight = grid::unit(14, "pt"),
                  title.position = "top",
                  label.theme = element_text(size = 16,
                                             family = FONT_FAMILY))
}

## Breaks and labels for every PIP colour scale, written without trailing
## zeros so that adjacent labels stay short.
pip_breaks <- c(0, 0.25, 0.5, 0.75, 1)
pip_labels <- c("0", "0.25", "0.5", "0.75", "1")

## Short pathway names for crowded axes: the CellPhoneDB prefix
## "Signaling by " is common to almost every class and carries no
## information on an axis already titled "Pathway".
short_pathway <- function(x) sub("^Signaling by ", "", x)

############################################################
## NOTATION HELPERS
##
## Each returns a plotmath expression. They are the only place in the
## repository where these symbols are spelled out.
############################################################
lab_beta_X_std  <- function() expression(hat(beta)[X]^{(s)})
lab_beta_XZ_std <- function() expression(hat(beta)[XZ]^{(s)})
lab_beta_X      <- function() expression(hat(beta)[X])
lab_beta_XZ     <- function() expression(hat(beta)[XZ])
lab_pip         <- function() expression(paste("PIP  (", hat(gamma), ")"))

## Parseable strings for facet strips (label_parsed), same notation.
lab_strip_abs_beta_X_std  <- function() "abs(hat(beta)[X]^{(s)})"
lab_strip_abs_beta_XZ_std <- function() "abs(hat(beta)[XZ]^{(s)})"

## Lookup for the parameter strips of the trace plot: the sampled quantities
## are the parameters themselves, so they carry no hat.
lab_strip_trace_parameters <- function() {
  c("Beta_X"               = "beta[X]",
    "Beta_XZ"              = "beta[XZ]",
    # The running-mean qualification is stated in the trace subtitle; in the
    # rotated strip it does not fit the panel height.
    "gamma (running mean)" = "gamma",
    # Quoted, because these strips are drawn with label_parsed and an
    # unquoted "Log-likelihood" would parse as a subtraction and render the
    # hyphen as a minus sign.
    "log-likelihood"       = "atop('Log-', 'likelihood')")
}

## Facet-row order for the trace figure, so that the log-likelihood appears
## above the inclusion indicator whichever subset is drawn. Any level not
## listed here keeps its position after those that are.
trace_parameter_levels <- function() {
  c("log-likelihood", "Beta_X", "Beta_XZ", "gamma (running mean)")
}

############################################################
## AXIS HELPERS
############################################################
## Break function returning at most n_max pretty breaks inside the limits.
## Used on long iteration axes where the default break count produces
## overlapping labels at the theme's font size.
breaks_at_most <- function(n_max = 4L) {
  function(lims) {
    for (n_target in seq(n_max - 1L, 1L)) {
      b <- pretty(lims, n = n_target)
      b <- b[b >= lims[1] & b <= lims[2]]
      if (length(b) <= n_max) return(b)
    }
    b
  }
}

## Iteration counts in thousands. The retained range is n_iter - burn_in,
## 98,000 per chain at the current protocol (100,000 iterations, 2,000
## burn-in). plot_trace() expands the axis by 6% on each side and ggplot2
## passes the expanded range to the break function, so the breaks render as
## 0k, 50k and 100k; the 100k break lies just beyond the last retained
## iteration, inside the expanded panel.
label_thousands <- scales::label_number(scale = 1e-3, suffix = "k")

############################################################
## EXPORT
############################################################
save_fig <- function(plot, name, width, height) {
  if (!isTRUE(capabilities("cairo"))) {
    stop("This R build has no cairo device, so fonts cannot be embedded ",
         "at export. Install a cairo-enabled R, or change save_fig() to ",
         "device = pdf and embed afterwards with ",
         "extrafont::embed_fonts().")
  }
  dir.create("Plots", showWarnings = FALSE)
  path <- file.path("Plots", paste0(name, ".pdf"))
  # units = "in" is ggsave's default and is stated here because the printed
  # font size depends on the width: see the canvas-width note in the header.
  ggsave(path, plot, width = width, height = height, units = "in",
         device = cairo_pdf)
  cat("Figure written:", path, "\n")
  invisible(path)
}

############################################################
## PLOT-BUILDING FUNCTIONS
############################################################
source("Github Codes/figure_functions.R")
