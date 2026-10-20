#
# Determination of Dissociation Constant (Kd) — v5
#
# Requirements:
#   install.packages(c("shiny", "bslib", "ggplot2", "minpack.lm"))
#   bslib >= 0.4.0, ggplot2 >= 3.5.0
#
# Run modes:
#   Server mode (default) — ExpN > 1000 triggers a download-and-run-locally warning:
#     Rscript -e "shiny::runApp('app.R')"
#
#   Local mode — ExpN of 5000, 10000 or 100000 allowed:
#     Shell:     BINDING_APP_LOCAL=true Rscript -e "shiny::runApp('app.R')"
#     R console: Sys.setenv(BINDING_APP_LOCAL = "true"); shiny::runApp("app.R")
#
# Changes in v5 (relative to v4):
#   1. Normality test: the KS test against a normal whose mean and SD were
#      estimated from the same data (KSgeneral::cont_ks_test, invalid p-values)
#      is replaced by a parametric bootstrap KS test with B = 200. The bootstrap
#      sample size equals the number of values actually tested. D and n are
#      reported with p. KSgeneral is no longer required.
#   2. Every test is wrapped in tryCatch and the whole simulation is guarded,
#      so an error is reported in the app instead of disconnecting the session.
#   3. Formula 1 uses the cancellation-free root 2*Pt*Lt / (S + sqrt(S^2 - 4*Pt*Lt)),
#      algebraically identical to (S - sqrt(S^2 - 4*Pt*Lt)) / 2.
#   4. Inputs are validated on the server. [P]total/Kd > 10 gives a warning;
#      [P]total/Kd > 1000 aborts the run.
#   5. Fit outcomes are recorded (normal, Kd at the lower limit 0 nM, not
#      converged, error). Not-converged and error fits are excluded from all
#      statistics; fits with Kd at the lower limit are kept on the linear scale
#      and excluded on log10.
#   6. Optional serial-dilution design (top concentration and dilution factor).
#      The theoretical fraction bound is shown for every ligand concentration.
#   7. Concentration (pipetting) error is no longer simulated; only random error
#      in the [PL] readout is simulated.
#   8. Fitted [P]total is no longer displayed; results focus on Kd.
#   9. Numbers are formatted by significant digits; median and 2.5-97.5% range added.
#  10. Number of Simulated Experiments: 1000, 5000, 10000 or 100000; values above
#      1000 run in local mode only.
#  11. Fitted Kd distributions can be shown as histograms (default) or as normal
#      Q-Q plots with a pointwise 95% band. Both are drawn square (height follows
#      the column width); the Q-Q plotting panel has a 1:1 aspect ratio.
#

library(shiny)
library(bslib)
library(ggplot2)
library(minpack.lm)

# Run-mode flag
# TRUE  → local run, no ExpN restriction.
# FALSE → server run (default), ExpN > 1000 shows a warning modal.
IS_LOCAL <- identical(Sys.getenv("BINDING_APP_LOCAL"), "true")

# Limits and settings
KD_PT_MIN    <- 0.001     # nM, allowed range for Kd and [P]total
KD_PT_MAX    <- 50000     # nM
TOP_MIN      <- 0.001     # nM, allowed range for the top ligand concentration
TOP_MAX      <- 1e7       # nM
RATIO_WARN   <- 10        # [P]total/Kd above this: warning (titration regime)
RATIO_MAX    <- 1000      # [P]total/Kd above this: run is refused
KS_B         <- 200L      # bootstrap replicates for the KS test
BOUNDARY_TOL <- 1e-6      # fitted Kd <= BOUNDARY_TOL * input Kd counts as the Kd = 0 bound

# Section 1 — Mathematical helper functions -------------------------------

#' Compute binding saturation curve data
#' Derives [L]total from [PL]: [L]t = ([P]t + Kd - [PL])*[PL] / ([P]t - [PL])
#' @param pt_const       [P]total (nM)
#' @param kd_const       Kd (nM)
#' @param plot_points    Smoothness (number of [PL] steps)
#' @param shown_points_num  Points displayed ([PL] from 0 to 99.0% of [P]total)
compute_binding_curve <- function(pt_const,
                                   kd_const,
                                   plot_points      = 1000L,
                                   shown_points_num = 991L) {
  num_array <- seq(0L, plot_points - 1L) / plot_points
  pl <- pt_const * num_array[seq_len(shown_points_num)]
  pt <- rep(pt_const, shown_points_num)
  kd <- rep(kd_const, shown_points_num)

  lt <- (pt + kd - pl) * pl / (pt - pl)
  lf <- lt - pl
  lt[1L] <- 0; lf[1L] <- 0   # Force exact zero at origin

  # Index 501 (R, 1-based) = Python index 500 → 50% occupancy → apparent Kd
  # Index 951 → 95% occupancy
  app_kd      <- lt[501L]
  lt_95       <- lt[951L]
  lf_95       <- lf[951L]
  fraction_bp <- pl / pt

  list(pt = pt, pl = pl, kd = kd, lt = lt, lf = lf,
       fraction_bp = fraction_bp,
       app_kd = app_kd, lt_95 = lt_95, lf_95 = lf_95)
}

#' Exact binding formula — Formula 1 (smaller root of the quadratic)
#' [PL]^2 - S*[PL] + Pt*Lt = 0 with S = Pt + Lt + Kd.
#' Written as 2*Pt*Lt / (S + sqrt(S^2 - 4*Pt*Lt)) rather than
#' (S - sqrt(S^2 - 4*Pt*Lt)) / 2: the two are algebraically identical, but the
#' second subtracts two nearly equal numbers when Pt << Kd + Lt and loses
#' precision (catastrophic cancellation), which corrupts numerical derivatives.
formula1 <- function(Lt, Kd, Pt) {
  S <- Pt + Lt + Kd
  2 * Pt * Lt / (S + sqrt(pmax(S^2 - 4 * Pt * Lt, 0)))
}

#' Approximate binding formula — Formula 2 ([L]free ≈ [L]total)
formula2 <- function(Lt, Kd, Pt) (Pt * Lt) / (Kd + Lt)

#' Levenberg-Marquardt fit via nlsLM() — Formula 1 (Exact)
#' Returns c(Kd, Pt, conv): conv = 1 converged, 0 not converged, NA failed.
fit_f1 <- function(lt_data, pl_obs, init_kd, init_pt, paraRatio = 0.5) {
  df <- data.frame(lt = lt_data, pl = pl_obs)
  tryCatch({
    fit <- suppressWarnings(nlsLM(
      pl ~ 2 * Pt * lt / ((Pt + lt + Kd) + sqrt(pmax((Pt + lt + Kd)^2 - 4 * Pt * lt, 0))),
      data = df, start = list(Kd = init_kd * paraRatio, Pt = init_pt * paraRatio),
      lower = c(Kd = 0, Pt = 0)))
    c(Kd = unname(coef(fit)["Kd"]), Pt = unname(coef(fit)["Pt"]),
      conv = as.numeric(isTRUE(fit$convInfo$isConv)))
  }, error = function(e) c(Kd = NA_real_, Pt = NA_real_, conv = NA_real_))
}

#' Levenberg-Marquardt fit via nlsLM() — Formula 2 (Approximate)
#' Returns c(Kd, Pt, conv) as fit_f1(); here Pt is the plateau (Bmax).
fit_f2 <- function(lt_data, pl_obs, init_kd, init_pt, paraRatio = 0.5) {
  df <- data.frame(lt = lt_data, pl = pl_obs)
  tryCatch({
    fit <- suppressWarnings(nlsLM(
      pl ~ Pt * lt / (Kd + lt),
      data = df, start = list(Kd = init_kd * paraRatio, Pt = init_pt * paraRatio),
      lower = c(Kd = 0, Pt = 0)))
    c(Kd = unname(coef(fit)["Kd"]), Pt = unname(coef(fit)["Pt"]),
      conv = as.numeric(isTRUE(fit$convInfo$isConv)))
  }, error = function(e) c(Kd = NA_real_, Pt = NA_real_, conv = NA_real_))
}

#' Classify each fit: "ok", "boundary" (Kd at the lower bound 0),
#' "nonconv" (iteration limit reached) or "failed" (nlsLM error).
fit_status <- function(kd_hat, conv, kd_input) {
  ifelse(is.na(kd_hat) | is.na(conv), "failed",
    ifelse(conv == 0, "nonconv",
      ifelse(kd_hat <= BOUNDARY_TOL * kd_input, "boundary", "ok")))
}

#' Ligand concentration series (ascending)
#' "default":  log-spaced from 0.001*Kd (at least 1e-4 nM) to [L]total at 99% bound
#' "dilution": top_conc / dil_factor^k, k = conc_pts-1, ..., 0
build_lt_series <- function(mode, kd_const, pt_const, conc_pts,
                            top_conc = NA_real_, dil_factor = NA_real_) {
  if (identical(mode, "dilution")) {
    top_conc / dil_factor^((conc_pts - 1L):0L)
  } else {
    lt_99 <- 99 * kd_const + 0.99 * pt_const
    10^seq(log10(max(1e-4, 0.001 * kd_const)), log10(lt_99), length.out = conc_pts)
  }
}

#' KS distance between x and the normal distribution fitted to x
#' (mean and SD estimated from x; the Lilliefors statistic).
#' Identical to ks.test(x, "pnorm", mean(x), sd(x))$statistic. Written out because
#' only D is needed: the p-value that ks.test reports assumes a fully specified
#' null and is not valid here, and ks.test warns on tied values (e.g. fits at Kd = 0).
ks_stat_normal <- function(x) {
  n <- length(x)
  s <- sd(x)
  if (n < 2L || !is.finite(s) || s <= 0) return(NA_real_)
  z  <- sort((x - mean(x)) / s)
  Fz <- pnorm(z)
  i  <- seq_len(n)
  max(Fz - (i - 1) / n, i / n - Fz)
}

#' Parametric bootstrap KS test of normality.
#' Because the normal is a location-scale family, the null distribution of the
#' statistic depends only on n, so bootstrap samples are drawn from N(0, 1) with
#' the same n as the data tested. One null distribution per n is cached in
#' `null_cache` (an environment) and shared by tests with equal n.
#' p = (1 + #{D* >= D}) / (B + 1); the smallest attainable p is 1/(B + 1).
ks_boot_test <- function(x, null_cache, B = KS_B) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 5L)
    return(list(n = n, D = NA_real_, p = NA_real_, msg = "fewer than 5 values"))
  D <- ks_stat_normal(x)
  if (is.na(D))
    return(list(n = n, D = NA_real_, p = NA_real_, msg = "zero variance"))
  key <- as.character(n)
  if (is.null(null_cache[[key]]))
    null_cache[[key]] <- vapply(seq_len(B),
                                function(b) ks_stat_normal(rnorm(n)), numeric(1))
  p <- (1 + sum(null_cache[[key]] >= D)) / (B + 1)
  list(n = n, D = D, p = p, msg = NULL)
}

# Section 2 — Formatting helpers ------------------------------------------

#' Format a number to `digits` significant digits
fmt <- function(x, digits = 4) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  sci <- abs(x) >= 1e6 || (x != 0 && abs(x) < 1e-3)
  format(signif(x, digits), scientific = sci, trim = TRUE)
}

#' Format a bootstrap p-value as an HTML fragment "<i>p</i> = x"; the smallest
#' attainable value 1/(B+1) is shown as "<i>p</i> < 0.005" (for B = 200)
fmt_p <- function(p, B = KS_B) {
  if (length(p) != 1L || is.na(p)) return("<i>p</i> = NA")
  if (p <= 1 / (B + 1) + 1e-12) return(sprintf("<i>p</i> &lt; %.3f", 1 / (B + 1)))
  paste0("<i>p</i> = ", formatC(p, format = "f", digits = 3))
}

# Section 3 — Plot functions ----------------------------------------------

#' Binding saturation curve
make_conc_plot <- function(data, scatter_df = NULL, log_x = FALSE) {
  df <- data.frame(lf = data$lf, lt = data$lt, pl = data$pl)
  if (log_x) {
    df_lf <- df[df$lf > 0, ]; df_lt <- df[df$lt > 0, ]
  } else {
    df_lf <- df_lt <- df
  }
  p <- ggplot() +
    geom_line(data = df_lf, aes(x = lf, y = pl, colour = "[L]free"),  linewidth = 1.0) +
    geom_line(data = df_lt, aes(x = lt, y = pl, colour = "[L]total"), linewidth = 1.0) +
    scale_colour_manual(values = c("[L]free" = "#66c2a5", "[L]total" = "#fc8d62"), name = NULL) +
    labs(x = "Concentration of Free/Total Ligand (nM)",
         y = "Concentration of Bound Protein ([PL], nM)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.98, 0.05),
          legend.justification   = c(1, 0),
          legend.text            = element_text(size = 13),
          legend.background    = element_rect(fill = alpha("white", 0.7), colour = NA),
          legend.key             = element_blank(),
          legend.key.width     = unit(1.5, "cm"))
  if (log_x) p <- p + scale_x_log10()
  if (!is.null(scatter_df) && nrow(scatter_df) > 0) {
    sc <- if (log_x) scatter_df[scatter_df$lt > 0, ] else scatter_df
    if (nrow(sc) > 0)
      p <- p + geom_point(data = sc, aes(x = lt, y = pl),
                          colour = "#e6550d", alpha = 0.70, size = 2.5,
                          inherit.aes = FALSE)
  }
  p
}

#' Density histogram with fitted normal PDF overlay
make_histogram_plot <- function(vec, title, x_lab) {
  vec <- vec[is.finite(vec)]
  if (length(vec) < 2L) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "Too few values to plot",
                      colour = "grey55", size = 3.5) +
             theme_void() + labs(title = title))
  }
  m <- mean(vec); s <- sd(vec)
  p <- ggplot(data.frame(x = vec), aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30L,
                   fill = "#3182bd", colour = "white", alpha = 0.65) +
    labs(title = title, x = x_lab, y = "Density") +
    theme_bw(base_size = 12)
  if (is.finite(s) && s > 0)
    p <- p + stat_function(fun = dnorm, args = list(mean = m, sd = s),
                           colour = "#de2d26", linewidth = 0.9)
  p
}

#' Normal Q-Q plot
#' Red line: the fitted normal (mean and SD of the plotted values), the same
#' reference distribution as the KS test. Grey band: pointwise 95% range of the
#' i-th order statistic of a normal sample of the same size, from
#' U(i) ~ Beta(i, n - i + 1); computed directly, so no random numbers are used.
#' About 5% of points can fall outside the band by chance for normal data.
make_qq_plot <- function(vec, title, y_lab) {
  vec <- sort(vec[is.finite(vec)])
  n <- length(vec)
  s <- if (n >= 2L) sd(vec) else NA_real_
  if (n < 3L || !is.finite(s) || s <= 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "Too few distinct values to plot",
                      colour = "grey55", size = 3.5) +
             theme_void() + labs(title = title))
  }
  m <- mean(vec)
  i <- seq_len(n)
  df <- data.frame(z  = qnorm((i - 0.5) / n),
                   y  = vec,
                   lo = m + s * qnorm(qbeta(0.025, i, n - i + 1)),
                   hi = m + s * qnorm(qbeta(0.975, i, n - i + 1)))
  ggplot(df, aes(x = z)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
    geom_abline(intercept = m, slope = s, colour = "#de2d26", linewidth = 0.9) +
    geom_point(aes(y = y), colour = "#3182bd", size = 1, alpha = 0.6) +
    labs(title = title, x = "Theoretical normal quantile", y = y_lab) +
    theme_bw(base_size = 12) +
    theme(aspect.ratio = 1)   # square plotting panel
}

# Section 4 — UI ----------------------------------------------------------

ui <- fluidPage(
  title = "Binding Curve Viewer — Determination of Kd",
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),

  tags$head(tags$style(HTML("
    body { padding: 0 14px 20px 14px; }
    .p1-header {
      padding: 10px 2px 8px 2px;
      border-bottom: 2px solid #2c3e50;
      margin-bottom: 16px;
    }
    .p1-header h1 { display:inline; font-size:1.42rem; font-weight:650; margin:0; }
    .p1-header .home-link {
      margin-left:14px; font-size:0.95rem; font-weight:500; white-space:nowrap;
    }
    .nav-bar {
      display:flex; align-items:center; justify-content:center;
      margin-top:10px; gap:14px;
    }
    .nav-info {
      text-align:center; min-width:160px; font-size:13px;
      padding:5px 12px; border:1px solid #adb5bd; border-radius:5px;
      background:#f8f9fa; line-height:1.5;
    }
    .sim-console {
      background:#f8f9fa; border:1px solid #dee2e6; border-radius:6px;
      padding:14px 18px; margin:14px 0 10px 0;
    }
    .res-panel {
      border:1px solid #dee2e6; border-radius:6px; padding:12px 16px;
      background:#ffffff; font-size:13px; min-height:90px; margin-bottom:10px;
      line-height:1.6;
    }
    .res-warn { color:#a04000; }
    .hist-label {
      font-weight:650; color:#495057; margin:12px 0 4px 2px;
      font-size:0.82rem; text-transform:uppercase; letter-spacing:0.06em;
    }
    .conc-table td, .conc-table th { padding:1px 8px !important; }
    .dist-bar {
      display:flex; align-items:center; flex-wrap:wrap; gap:6px 22px;
      margin:12px 0 4px 0;
    }
    .dist-bar .hist-label { margin:0 0 0 2px; }
    .dist-bar .form-group, .dist-bar .shiny-input-container { margin-bottom:0; }
    .dist-note { font-size:12.5px; color:#555; margin:2px 0 0 2px; }
  "))),

  # p1_header
  div(class = "p1-header",
    h1(HTML("Determination of Dissociation Constant (<i>K</i><sub>d</sub>)")),
    tags$a(class = "home-link", href = "/binding-curve-viewer/", "[Back to Home]")
  ),

  # Main row: config_col | tab_conc + nav | data_table card
  # Plots and Binding Parameters are NOT updated reactively; they refresh only
  # when "Start Simulation" is clicked.
  fluidRow(

    # config_col: only Kd and [P]total inputs
    column(3,
      numericInput("kd_val",
        HTML("Equilibrium Dissociation Constant (<b><i>K</i><sub>d</sub></b>, nM)
              <br><small style='color:#666;font-weight:400'>0.001 – 50 000 nM</small>"),
        value = 100, min = KD_PT_MIN, max = KD_PT_MAX, step = 1),
      numericInput("pt_val",
        HTML("Total Fixed Component Conc. (<b>[P]<sub>total</sub></b>, nM)
              <br><small style='color:#666;font-weight:400'>0.001 – 50 000 nM;
              [P]<sub>total</sub>/<i>K</i><sub>d</sub> &le; 1000</small>"),
        value = 100, min = KD_PT_MIN, max = KD_PT_MAX, step = 1)
    ),

    # tab_conc: log and linear plots + navigation arrows
    column(5,
      tabsetPanel(id = "tab_conc",
        tabPanel("Logarithmic X-axis", plotOutput("plot_log",    height = "380px")),
        tabPanel("Linear X-axis",      plotOutput("plot_linear", height = "380px"))
      ),
      div(class = "nav-bar",
        actionButton("nav_left",  HTML("&#9664;"), class = "btn btn-primary btn-sm"),
        uiOutput("nav_text"),
        actionButton("nav_right", HTML("&#9654;"), class = "btn btn-primary btn-sm")
      )
    ),

    # data_table card: binding parameters + (after simulation) ligand series
    column(4,
      bslib::card(
        height = 500,
        bslib::card_header(HTML("<b>Binding Parameters</b>")),
        bslib::card_body(uiOutput("data_table"))
      )
    )
  ),

  # Monte Carlo Simulation Console
  fluidRow(column(12,
    div(class = "sim-console",
      tags$h5("Monte Carlo Simulation",
              style = "margin:0 0 12px 0; font-weight:700; color:#2b579a;"),

      # Row 1: ExpN | Ligand concentrations | Replicates | [PL] noise | Random seed
      fluidRow(
        column(2,
          selectInput("expn",
            "Number of Simulated Experiments",
            choices  = c("1000" = "1000", "5000" = "5000", "10,000" = "10000", "100,000" = "100000"),
            selected = "1000")
        ),
        column(2,
          selectInput("conc_pts",
            "Number of Ligand Concentrations per Experiment",
            choices  = as.character(5:10),
            selected = "9")
        ),
        column(2,
          selectInput("repn",
            "Replicates per Ligand Concentration",
            choices  = as.character(3:10),
            selected = "5")
        ),
        column(2,
          selectInput("ratio_sd",
            HTML("[PL] Noise SD as Fraction of [P]<sub>total</sub>"),
            choices  = c("0.0125", "0.0250", "0.0500"),
            selected = "0.0250")
        ),
        column(2,
          numericInput("seed_val",
            HTML("Random Seed <small style='color:#666;font-weight:400'>(positive integer)</small>"),
            value = 1, min = 1, step = 1)
        )
      ),

      # Row 2: ligand concentration design | Start + Clear buttons
      fluidRow(
        column(4,
          radioButtons("design_mode", "Ligand Concentration Series",
            choiceNames = list(
              HTML("Default: log-spaced, 0.001&times;<i>K</i><sub>d</sub> to
                    [L]<sub>total</sub> at 99% bound"),
              HTML("Serial dilution: top concentration and dilution factor")),
            choiceValues = list("default", "dilution"),
            selected = "default")
        ),
        column(2,
          conditionalPanel("input.design_mode == 'dilution'",
            numericInput("top_conc",
              HTML("Top [L]<sub>total</sub> (nM)"),
              value = 10000, min = TOP_MIN, max = TOP_MAX))
        ),
        column(2,
          conditionalPanel("input.design_mode == 'dilution'",
            selectInput("dil_factor", "Dilution Factor",
              choices  = c("2" = "2", "3" = "3",
                           "3.16 (half-log)" = as.character(sqrt(10)),
                           "10" = "10"),
              selected = "3"))
        ),
        column(3,
          tags$div(style = "margin-top:25px;",
            actionButton("start_sim", "Start Simulation",
              class = "btn btn-success btn-sm", style = "margin-right:6px;"),
            actionButton("clear_sim", "Clear Simulation",
              class = "btn btn-danger btn-sm")
          )
        )
      ),
      conditionalPanel("input.design_mode == 'dilution'", uiOutput("design_hint")),

      # Progress message (status text after completion)
      uiOutput("progress_ui")
    )
  )),

  # Simulation results
  fluidRow(
    column(6, div(class = "res-panel", uiOutput("res_f1"))),
    column(6, div(class = "res-panel", uiOutput("res_f2")))
  ),

  # Distribution grid — Kd distributions for both formulas, shown as
  # histograms (default) or normal Q-Q plots. Switching only changes the display.
  div(class = "dist-bar",
    div(class = "hist-label", HTML(
      "Fitted <span style='text-transform:none'><i>K</i><sub>d</sub></span> Distributions &mdash;
       Formula 1 (Exact) and Formula 2 (Approx)")),
    radioButtons("dist_plot", NULL, inline = TRUE,
      choices  = c("Histogram" = "hist", "Normal Q-Q plot" = "qq"),
      selected = "hist")
  ),
  # height = "auto": each plot is drawn as tall as it is wide (see make_dist_output)
  fluidRow(
    column(3, plotOutput("hist_kd1",  height = "auto")),
    column(3, plotOutput("hist_lkd1", height = "auto")),
    column(3, plotOutput("hist_kd2",  height = "auto")),
    column(3, plotOutput("hist_lkd2", height = "auto"))
  ),
  conditionalPanel("input.dist_plot == 'hist'",
    tags$p(class = "dist-note",
      "Red curve: normal density with the mean and SD of the plotted values.")),
  conditionalPanel("input.dist_plot == 'qq'",
    tags$p(class = "dist-note", paste(
      "Red line: normal distribution with the mean and SD of the plotted values",
      "(the same reference as the KS test). Grey band: pointwise 95% range for a",
      "normal sample of the same size; about 5% of points can fall outside it by",
      "chance even when the values are normal."))),
  tags$br()
)

# Section 5 — Server ------------------------------------------------------

server <- function(input, output, session) {

  # Reactive state
  sim <- reactiveValues(
    # Committed binding-curve data.
    # Initialised with defaults (Kd=100, Pt=100) so plots show on page load.
    # Updated ONLY when "Start Simulation" is clicked — not on every input change.
    bd           = compute_binding_curve(100, 100),

    # Raw scatter data for first ≤10 experiments (for navigation overlay)
    all_exps     = NULL,    # data.frame(exp_idx, lt, pl)

    # Fitted Kd from all ExpN experiments and the status of each fit
    kd1_all      = NULL,
    kd2_all      = NULL,
    st1          = NULL,    # "ok" / "boundary" / "nonconv" / "failed"
    st2          = NULL,

    # Bootstrap KS results (NULL means "simulation not run yet")
    ks_results   = NULL,

    # Ligand concentration series shown in the card after simulation
    lt_nominal   = NULL,

    # Navigation state
    current_idx  = 1L,
    total_expn   = 0L,
    status_msg   = "Simulation not started."
  )

  clear_results <- function() {
    sim$all_exps   <- NULL
    sim$kd1_all    <- NULL;  sim$kd2_all    <- NULL
    sim$st1        <- NULL;  sim$st2        <- NULL
    sim$ks_results <- NULL;  sim$lt_nominal <- NULL
    sim$total_expn <- 0L;    sim$current_idx <- 1L
  }

  # Current scatter points (reactive on navigation, not on inputs)
  current_sc <- reactive({
    if (is.null(sim$all_exps)) return(NULL)
    sim$all_exps[sim$all_exps$exp_idx == sim$current_idx, ]
  })

  # Binding curve plots
  # Both plots depend on sim$bd (committed values) and current_sc() (navigation).
  output$plot_log <- renderPlot({
    make_conc_plot(sim$bd, current_sc(), log_x = TRUE)
  })
  output$plot_linear <- renderPlot({
    make_conc_plot(sim$bd, current_sc(), log_x = FALSE)
  })

  # Binding Parameters card
  # Depends only on sim$bd (committed) and sim$lt_nominal.
  output$data_table <- renderUI({
    d <- sim$bd

    param_html <- paste0(
      "<table class='table table-sm table-borderless mb-0' style='font-size:13.5px;'>",
      "<tr>",
        "<td>Apparent <b><i>K</i><sub>d</sub></b></td>",
        "<td style='padding:2px 6px;'>=</td>",
        "<td><b>", formatC(d$app_kd, format = "f", digits = 5), " nM</b></td>",
      "</tr><tr>",
        "<td><b><i>K</i><sub>d</sub></b></td>",
        "<td style='padding:2px 6px;'>=</td>",
        "<td>", formatC(d$kd[1L], format = "f", digits = 5), " nM</td>",
      "</tr><tr>",
        "<td>App. <i>K</i><sub>d</sub> / <i>K</i><sub>d</sub></td>",
        "<td style='padding:2px 6px;'>=</td>",
        "<td>", formatC(d$app_kd / d$kd[1L], format = "f", digits = 2), "</td>",
      "</tr>",
      "<tr><td colspan='3'><hr style='margin:8px 0;'></td></tr>",
      "<tr><td colspan='3'>",
        "<b>[L]<sub>free</sub></b> = ", formatC(d$lf_95, format = "f", digits = 5), " nM",
        "<br><small style='color:#555'>when 95% protein is bound</small>",
      "</td></tr>",
      "<tr><td colspan='3'>",
        "<b>[L]<sub>total</sub></b> = ", formatC(d$lt_95, format = "f", digits = 5), " nM",
        "<br><small style='color:#555'>when 95% protein is bound</small>",
      "</td></tr>",
      "</table>"
    )

    # Ligand series with the theoretical fraction bound — shown after simulation
    conc_html <- if (!is.null(sim$lt_nominal)) {
      lt  <- sort(unique(sim$lt_nominal))
      occ <- formula1(lt, d$kd[1L], d$pt[1L]) / d$pt[1L]
      informative <- occ >= 0.05 & occ <= 0.95
      rows <- paste0(
        "<tr", ifelse(informative, " style='font-weight:700;'", ""), ">",
        "<td>", vapply(lt, fmt, character(1)), "</td>",
        "<td>", sprintf("%.1f%%", 100 * occ), "</td></tr>",
        collapse = "")
      paste0(
        "<hr style='margin:8px 0;'>",
        "<span style='display:block; font-size:13px; font-weight:bold; margin:0 0 4px 0;'>",
        "Ligand Concentrations Used</span>",
        "<table class='table table-sm mb-1 conc-table' style='font-size:12.5px;'>",
        "<thead><tr><th>[L]<sub>total</sub> (nM)</th><th>Fraction bound</th></tr></thead>",
        "<tbody>", rows, "</tbody></table>",
        "<small style='color:#555'>", sum(informative), " of ", length(lt),
        " concentrations lie between 5% and 95% bound (bold).</small>"
      )
    } else ""

    HTML(paste0(param_html, conc_html))
  })

  # Hint for the serial-dilution design (live, not committed)
  output$design_hint <- renderUI({
    kd  <- input$kd_val;   pt <- input$pt_val
    top <- input$top_conc; f  <- suppressWarnings(as.numeric(input$dil_factor))
    n   <- as.integer(input$conc_pts)
    ok  <- function(v) length(v) == 1L && is.numeric(v) && is.finite(v) && v > 0
    if (!ok(kd) || !ok(pt)) return(NULL)
    lt99 <- 99 * kd + 0.99 * pt
    series_txt <- if (ok(top) && ok(f))
      sprintf(" Series: %s &rarr; %s nM (%d points, %s-fold steps).",
              fmt(top / f^(n - 1L)), fmt(top), n, fmt(f, 3)) else ""
    tags$p(style = "font-size:12.5px; color:#555; margin:0 0 4px 0;",
      HTML(sprintf("For the current inputs, [L]<sub>total</sub> at 99%% bound = %s nM.%s",
                   fmt(lt99), series_txt)))
  })

  # Navigation display
  output$nav_text <- renderUI({
    max_shown <- min(10L, sim$total_expn)
    if (max_shown == 0L)
      div(class = "nav-info", "Current Exp No.", tags$br(), tags$strong("0 / 0"))
    else
      div(class = "nav-info", "Current Exp No.", tags$br(),
          tags$strong(paste0(sim$current_idx, " / ", max_shown)))
  })
  observeEvent(input$nav_left, {
    max_idx <- min(10L, sim$total_expn)
    if (max_idx > 0L)
      sim$current_idx <- (sim$current_idx - 2L) %% max_idx + 1L
  })
  observeEvent(input$nav_right, {
    max_idx <- min(10L, sim$total_expn)
    if (max_idx > 0L)
      sim$current_idx <- sim$current_idx %% max_idx + 1L
  })

  # Progress / status message
  output$progress_ui <- renderUI({
    tags$p(HTML(sim$status_msg),
           style = "font-weight:bold; color:#1a5276; margin:10px 0 0 0; font-size:13px;")
  })

  # Results panels
  # KS results are read from sim$ks_results, computed inside the simulation
  # observer (not reactively here).
  ks_line <- function(label, r) {
    if (is.null(r) || is.na(r$D))
      return(paste0(label, ": not available",
                    if (!is.null(r$msg)) paste0(" (", r$msg, ")") else ""))
    paste0(label, ": <i>D</i> = ", fmt(r$D, 3), ", ", fmt_p(r$p),
           " (<i>n</i> = ", r$n, ")")
  }

  make_results_html <- function(label, kd_vec, status, ks_lin, ks_log) {
    n_all <- length(status)
    cnt <- table(factor(status, levels = c("ok", "boundary", "nonconv", "failed")))
    kd_used <- kd_vec[status %in% c("ok", "boundary")]
    stats_html <- if (length(kd_used) >= 2L) {
      q <- quantile(kd_used, c(0.025, 0.5, 0.975), names = FALSE)
      paste0(
        "Fitted <i>K</i><sub>d</sub> (nM), <i>n</i> = ", length(kd_used), ": ",
        "mean &plusmn; SD = ", fmt(mean(kd_used)), " &plusmn; ", fmt(sd(kd_used)),
        "; median = ", fmt(q[2]),
        " (2.5&ndash;97.5%: ", fmt(q[1]), "&ndash;", fmt(q[3]), ")<br>")
    } else "Fitted <i>K</i><sub>d</sub>: too few usable fits<br>"

    warn <- character(0)
    if (cnt[["boundary"]] > 0L)
      warn <- c(warn, sprintf(paste0(
        "%s%% of fits have <i>K</i><sub>d</sub> at the lower limit (0 nM): ",
        "kept on the linear scale, excluded on log10"),
        fmt(100 * cnt[["boundary"]] / n_all, 3)))
    if (cnt[["nonconv"]] + cnt[["failed"]] > 0L)
      warn <- c(warn, sprintf(paste0(
        "%d fits that did not converge and %d fits that ended in an error ",
        "were excluded from all statistics"),
        cnt[["nonconv"]], cnt[["failed"]]))
    warn_html <- if (length(warn))
      paste0("<span class='res-warn'>&#9888; ", paste(warn, collapse = "; "), ".</span>") else ""

    # Fit outcomes: counts first, so "0 nM" never sits next to a count
    outcome_html <- if (cnt[["ok"]] == n_all) {
      sprintf("Fit outcomes: all %d fits converged normally.<br>", n_all)
    } else {
      sprintf(paste0(
        "Fit outcomes (of %d): %d normal &middot; %d with <i>K</i><sub>d</sub> ",
        "at the lower limit (0 nM) &middot; %d not converged &middot; %d errors<br>"),
        n_all, cnt[["ok"]], cnt[["boundary"]], cnt[["nonconv"]], cnt[["failed"]])
    }

    HTML(paste0(
      "<b>", label, " Results:</b><br>",
      outcome_html,
      stats_html,
      "Normality test (Kolmogorov&ndash;Smirnov, mean and SD estimated from the data; ",
      "parametric bootstrap, <i>B</i> = ", KS_B, "):<br>",
      "&nbsp;&nbsp;", ks_line("<i>K</i><sub>d</sub>", ks_lin), "<br>",
      "&nbsp;&nbsp;", ks_line("log10(<i>K</i><sub>d</sub>)", ks_log), "<br>",
      warn_html
    ))
  }

  output$res_f1 <- renderUI({
    if (is.null(sim$kd1_all) || is.null(sim$ks_results))
      return(HTML("<b>Formula 1 (Exact) Results:</b><br>Simulation not started."))
    ks <- sim$ks_results
    make_results_html("Formula 1 (Exact)", sim$kd1_all, sim$st1, ks$f1_lin, ks$f1_log)
  })
  output$res_f2 <- renderUI({
    if (is.null(sim$kd2_all) || is.null(sim$ks_results))
      return(HTML("<b>Formula 2 (Approx) Results:</b><br>Simulation not started."))
    ks <- sim$ks_results
    make_results_html("Formula 2 (Approx)", sim$kd2_all, sim$st2, ks$f2_lin, ks$f2_log)
  })

  # Distribution outputs (histogram or normal Q-Q plot, per input$dist_plot)
  # Linear scale: fits with Kd > 0 plus those at the lower limit (Kd = 0).
  # log10 scale:  fits with Kd > 0 only.
  make_dist_output <- function(out_id, which, log_scale) {
    output[[out_id]] <- renderPlot({
      kd <- if (which == 1L) sim$kd1_all else sim$kd2_all
      st <- if (which == 1L) sim$st1     else sim$st2
      if (is.null(kd)) return(NULL)
      tag <- paste0("F", which)
      if (log_scale) {
        vec   <- log10(kd[st == "ok"])
        title <- sprintf("%s: log10(Kd)  (n = %d)", tag, length(vec))
        lab   <- "log10(Kd (nM))"
      } else {
        vec   <- kd[st %in% c("ok", "boundary")]
        title <- sprintf("%s: Kd  (n = %d)", tag, length(vec))
        lab   <- "Kd (nM)"
      }
      if (identical(input$dist_plot, "qq")) make_qq_plot(vec, title, lab)
      else make_histogram_plot(vec, title, lab)
    },
    # Square image: height follows the current width of the output, so the
    # plots stay square at any window size (histograms and Q-Q plots alike)
    height = function() {
      w <- session$clientData[[paste0("output_", out_id, "_width")]]
      if (is.null(w) || !is.finite(w) || w <= 0) 300 else w
    })
  }
  make_dist_output("hist_kd1",  1L, FALSE)
  make_dist_output("hist_lkd1", 1L, TRUE)
  make_dist_output("hist_kd2",  2L, FALSE)
  make_dist_output("hist_lkd2", 2L, TRUE)

  # Input validation; returns NULL when all inputs are valid
  check_inputs <- function(kd, pt, seed, mode, top) {
    is_num <- function(v) length(v) == 1L && is.numeric(v) && is.finite(v)
    msgs <- character(0)
    if (!is_num(kd) || kd < KD_PT_MIN || kd > KD_PT_MAX)
      msgs <- c(msgs, "<i>K</i><sub>d</sub> must be between 0.001 and 50 000 nM.")
    if (!is_num(pt) || pt < KD_PT_MIN || pt > KD_PT_MAX)
      msgs <- c(msgs, "[P]<sub>total</sub> must be between 0.001 and 50 000 nM.")
    if (!is_num(seed) || seed < 1 || seed > .Machine$integer.max || seed != round(seed))
      msgs <- c(msgs, "Random seed must be a positive integer.")
    if (identical(mode, "dilution") &&
        (!is_num(top) || top < TOP_MIN || top > TOP_MAX))
      msgs <- c(msgs, "Top [L]<sub>total</sub> must be between 0.001 and 10 000 000 nM.")
    if (length(msgs)) paste(msgs, collapse = "<br>") else NULL
  }

  # Start Simulation
  observeEvent(input$start_sim, {

    expN <- as.integer(input$expn)

    # Server-mode guard: warn when ExpN > 1000
    # In local mode (IS_LOCAL = TRUE) this check is skipped entirely.
    if (!IS_LOCAL && expN > 1000L) {
      showModal(modalDialog(
        title = tags$span(
          style = "color:#c0392b; font-weight:700;",
          "\u26a0\ufe0f  Large Simulation \u2014 Server Load Warning"
        ),
        HTML(paste0(
          "Running <b>Number of Simulated Experiments &gt; 1000</b> on the web server ",
          "may cause heavy load and slow the application for all users.<br><br>",
          "Please <b>download and run the app locally</b>:<br>",
          "<a href='https://github.com/ydu-sci/Binding_Curve_Viewer' target='_blank' ",
          "style='font-weight:600;'>github.com/ydu-sci/Binding_Curve_Viewer</a><br>(app file: <code>src/shiny-apps/kd-sim/app.R</code>)<br><br>",
          "Once downloaded, launch it from that folder in <b>local mode</b> (no restriction):<br>",
          "<code style='background:#f1f1f1; padding:3px 6px; border-radius:3px;'>",
          "BINDING_APP_LOCAL=true Rscript -e \"shiny::runApp('app.R')\"</code><br>",
          "<small style='color:#555;'>Or from R console: ",
          "<code>Sys.setenv(BINDING_APP_LOCAL=\"true\"); shiny::runApp(\"app.R\")</code></small>"
        )),
        footer = modalButton("OK"),
        easyClose = TRUE
      ))
      return()   # Abort — do not start simulation
    }

    # Validate inputs
    design_mode <- input$design_mode
    err <- check_inputs(input$kd_val, input$pt_val, input$seed_val,
                        design_mode, input$top_conc)
    if (!is.null(err)) {
      showNotification(HTML(err), type = "error", duration = 8)
      return()
    }

    repN     <- as.integer(input$repn)
    ratio    <- as.numeric(input$ratio_sd)
    conc_pts <- as.integer(input$conc_pts)
    kd_const <- as.numeric(input$kd_val)
    pt_const <- as.numeric(input$pt_val)
    seed     <- as.integer(input$seed_val)
    top_conc <- if (identical(design_mode, "dilution")) as.numeric(input$top_conc) else NA_real_
    dil_fac  <- if (identical(design_mode, "dilution")) as.numeric(input$dil_factor) else NA_real_

    # [P]total/Kd limits: > RATIO_MAX refused, > RATIO_WARN warned
    ratio_pk <- pt_const / kd_const
    if (ratio_pk > RATIO_MAX) {
      showModal(modalDialog(
        title = tags$span(style = "color:#c0392b; font-weight:700;",
                          HTML("[P]<sub>total</sub> / <i>K</i><sub>d</sub> too high")),
        HTML(sprintf(paste0(
          "[P]<sub>total</sub> / <i>K</i><sub>d</sub> = %s exceeds %d.<br><br>",
          "Under such strong titration conditions the binding curve reflects ",
          "[P]<sub>total</sub> rather than <i>K</i><sub>d</sub>, and <i>K</i><sub>d</sub> ",
          "cannot be estimated. Please lower [P]<sub>total</sub> or raise ",
          "<i>K</i><sub>d</sub>."), fmt(ratio_pk, 3), RATIO_MAX)),
        footer = modalButton("OK"),
        easyClose = TRUE
      ))
      return()
    }
    ratio_note <- NULL
    if (ratio_pk > RATIO_WARN) {
      ratio_note <- sprintf(paste0(
        "[P]<sub>total</sub> / <i>K</i><sub>d</sub> = %s (&gt; %d): titration regime. ",
        "Fitted <i>K</i><sub>d</sub> values are imprecise and some fits may reach ",
        "<i>K</i><sub>d</sub> = 0; consider lowering [P]<sub>total</sub> below ",
        "<i>K</i><sub>d</sub>."), fmt(ratio_pk, 3), RATIO_WARN)
      showNotification(HTML(ratio_note), type = "warning", duration = 12)
    }

    # Commit binding-curve parameters and reset previous results
    sim$bd <- compute_binding_curve(pt_const, kd_const)
    clear_results()
    sim$status_msg <- sprintf("Running (ExpN=%d, seed=%d) ...", expN, seed)

    tryCatch(
      withProgress(
        message = sprintf("Monte Carlo Simulation (ExpN = %d)", expN),
        value   = 0,
        {
          # Set random seed
          set.seed(seed)

          # Ligand concentration series (nominal = true; no concentration error)
          lt_nominal <- build_lt_series(design_mode, kd_const, pt_const, conc_pts,
                                        top_conc, dil_fac)
          lt_array <- rep(lt_nominal, each = repN)
          n_pts    <- length(lt_array)
          noise_sd <- ratio * pt_const
          n_total  <- n_pts * expN   # repN × conc_pts × expN

          # Noise-free [PL] (identical in every experiment) + additive [PL]
          # readout noise. The only random draws are the n_total noise values,
          # in the same order as in v4 with zero concentration error.
          pl_true     <- formula1(lt_array, kd_const, pt_const)
          noisy_array <- array(rep(pl_true, expN) + rnorm(n_total, sd = noise_sd),
                               dim = c(repN, conc_pts, expN))

          # Storage
          store_idx <- integer(0); store_lt <- numeric(0); store_pl <- numeric(0)
          kd1_vec <- numeric(expN); conv1 <- numeric(expN)
          kd2_vec <- numeric(expN); conv2 <- numeric(expN)
          upd_every <- max(1L, expN %/% 100L)

          # Main simulation loop (progress: 0 → 88%)
          for (i in seq_len(expN)) {
            noisy_pl <- as.vector(noisy_array[, , i])

            # Store raw data only for first 10 experiments (scatter overlay)
            if (i <= 10L) {
              store_idx <- c(store_idx, rep(i, n_pts))
              store_lt  <- c(store_lt,  lt_array)
              store_pl  <- c(store_pl,  noisy_pl)
            }

            r1 <- fit_f1(lt_array, noisy_pl, kd_const, pt_const)
            kd1_vec[i] <- r1[["Kd"]]; conv1[i] <- r1[["conv"]]

            r2 <- fit_f2(lt_array, noisy_pl, kd_const, pt_const)
            kd2_vec[i] <- r2[["Kd"]]; conv2[i] <- r2[["conv"]]

            if (i %% upd_every == 0L || i == expN)
              incProgress(0.88 * upd_every / expN,
                          detail = sprintf("Exp %d / %d  (%.0f%%)",
                                           i, expN, 100 * i / expN))
          }

          st1 <- fit_status(kd1_vec, conv1, kd_const)
          st2 <- fit_status(kd2_vec, conv2, kd_const)

          # Bootstrap KS tests (progress: 88% → 100%)
          # Linear scale: converged fits incl. those at the Kd = 0 bound.
          # log10 scale:  converged fits excl. those at the bound.
          null_cache <- new.env(parent = emptyenv())
          ks_specs <- list(
            f1_lin = list(x = kd1_vec[st1 %in% c("ok", "boundary")], name = "F1: Kd"),
            f1_log = list(x = log10(kd1_vec[st1 == "ok"]),            name = "F1: log10(Kd)"),
            f2_lin = list(x = kd2_vec[st2 %in% c("ok", "boundary")], name = "F2: Kd"),
            f2_log = list(x = log10(kd2_vec[st2 == "ok"]),            name = "F2: log10(Kd)")
          )
          ks_res <- list()
          for (k in names(ks_specs)) {
            incProgress(0.03, detail = paste0("KS test: ", ks_specs[[k]]$name, " ..."))
            ks_res[[k]] <- tryCatch(
              ks_boot_test(ks_specs[[k]]$x, null_cache),
              error = function(e) list(n = length(ks_specs[[k]]$x), D = NA_real_,
                                       p = NA_real_, msg = conditionMessage(e)))
          }

          # Commit all results
          sim$all_exps   <- data.frame(exp_idx = store_idx,
                                       lt      = store_lt,
                                       pl      = store_pl)
          sim$kd1_all    <- kd1_vec;  sim$st1 <- st1
          sim$kd2_all    <- kd2_vec;  sim$st2 <- st2
          sim$lt_nominal <- lt_nominal
          sim$ks_results <- ks_res
          sim$total_expn  <- expN
          sim$current_idx <- 1L
          sim$status_msg  <- paste0(
            sprintf("Completed: %d runs, seed = %d.", expN, seed),
            if (!is.null(ratio_note))
              paste0("<br><span style='color:#a04000; font-weight:400;'>&#9888; ",
                     ratio_note, "</span>") else "")
        }
      ),
      error = function(e) {
        clear_results()
        sim$status_msg <- paste0(
          "<span style='color:#c0392b;'>Simulation stopped with an error: ",
          htmltools::htmlEscape(conditionMessage(e)), "</span>")
        showNotification(paste("Simulation error:", conditionMessage(e)),
                         type = "error", duration = 10)
      }
    )
  })

  # Clear Simulation
  # Clears all simulation results including the ligand list in the card.
  # sim$bd (committed binding curve) is preserved.
  observeEvent(input$clear_sim, {
    clear_results()
    sim$status_msg <- "Simulation cleared."
  })
}

# Section 6 — Launch ------------------------------------------------------

shinyApp(ui = ui, server = server)
