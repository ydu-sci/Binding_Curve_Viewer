#
# Competitive Binding — Inhibition Curve + MC Simulation — v5
#
# Pure R Shiny implementation of the competitive binding equilibrium IC50 curve
# and Monte Carlo simulation.
#
# Key mathematical references:
#   formula_competitive_pl : FEBS Lett. 1995, 360(2), 111-4
#   inhibition_calc        : paper SI, exact IC50 at 50% inhibition
#   calc_ki_app1           : Cheng-Prusoff, Biochem. Pharmacol. 1973, 22, 3099
#   calc_ki_exact          : Anal. Biochem. 2004, 332(2), 261-73
#   logistic_func          : 4-parameter logistic for apparent IC50 fitting
#
# Requirements:
#   install.packages(c("shiny", "bslib", "ggplot2", "minpack.lm"))
#   bslib >= 0.4.0, ggplot2 >= 3.5.0
#
# Run modes (same convention as the other apps):
#   Server mode (ExpN > 1000 shows a warning):
#     Rscript -e "shiny::runApp('app.R')"
#   Local mode (ExpN of 5000, 10000 or 100000):
#     BINDING_APP_LOCAL=true Rscript -e "shiny::runApp('app.R')"
#
# Changes in v5 (relative to v4):
#   1. Equilibrium [P]free: the trigonometric solution of the Wang (1995) cubic is
#      used as a starting value and refined by bounded Newton steps on the mass
#      balance. In v4 the trigonometric solution alone lost all precision when
#      [I]0 or [L]0 >> [P]0 (wrong simulated data, and mechanistic fits that
#      stayed at their starting value).
#   2. The [PL] quadratic in inhibition_calc() and calc_ki_exact() uses the
#      cancellation-free root 2*P0*L0 / (S + sqrt(S^2 - 4*P0*L0)).
#   3. Default inhibitor series is log-spaced from apparent IC50/1000 to apparent
#      IC50*1000 (v4 used the free-inhibitor IC50, which missed 50% inhibition in
#      tight binding). Optional serial-dilution series (top [I]0, dilution
#      factor). The theoretical inhibition is shown for every [I]0.
#   4. [P]0/Ki > 10 gives a tight-binding warning; [P]0/Ki > 1000 aborts the run.
#   5. Fit outcomes are recorded (normal, at the lower bound 0, not converged,
#      error). Not-converged and error fits are excluded from all statistics;
#      fits at the lower bound are kept on the linear scale, excluded on log10.
#   6. Logistic starting values: apparent IC50 at 0.5 x the true value (as for
#      all other fits), lower asymptote from the data at the highest [I]0.
#   7. Ki from each experiment's fitted apparent IC50 by Wang's group equation
#      (Nikolovska-Coleska et al. 2004), besides Cheng-Prusoff.
#   8. Normality test: parametric bootstrap KS test (B = 200) replaces
#      KSgeneral::cont_ks_test; KSgeneral is no longer required.
#   9. Concentration (pipetting) errors removed; default [PL] noise 0.025.
#  10. Inputs validated on the server; tryCatch guards the simulation; ExpN
#      1000, 5000, 10000 or 100000 (above 1000 in local mode only); results with fit outcomes, median and
#      2.5-97.5% range; histograms or normal Q-Q plots, drawn square;
#      significant-digit formatting; "Back to Home" link.
#

library(shiny)
library(bslib)
library(ggplot2)
library(minpack.lm)

IS_LOCAL <- identical(Sys.getenv("BINDING_APP_LOCAL"), "true")

# Limits and settings
CONC_MIN <- 0.001;  CONC_MAX <- 50000   # nM, for [P]0, [L]0, Kd, Ki
TOP_MIN  <- 0.001;  TOP_MAX  <- 1e9     # nM, top [I]0 of a serial dilution
RATIO_WARN   <- 10      # [P]0/Ki above this: tight-binding warning
RATIO_MAX    <- 1000    # [P]0/Ki above this: run is refused
KS_B         <- 200L    # bootstrap replicates for the KS test
BOUNDARY_TOL <- 1e-6    # estimate <= BOUNDARY_TOL * true value counts as the lower bound 0
DEFAULTS <- list(p0 = 100, l0 = 75, kd = 100, ki = 100)

# Section 1 — Mathematical helper functions -------------------------------

#' Free protein [P] at equilibrium in the three-component competitive system.
#' Starting value: trigonometric solution of the cubic (FEBS Lett. 1995, 360(2),
#' 111-4). It is then refined by bounded Newton steps on the mass balance
#'   f(P) = P + P*L0/(Kd + P) + P*I0/(Ki + P) - P0 = 0,
#' which is increasing and concave in P, so the root lies in
#' [P0/(1 + L0/Kd + I0/Ki), P0]. The trigonometric form alone loses precision
#' when [I]0 or [L]0 >> [P]0 (subtraction of nearly equal numbers, acos near +-1).
#' Vectorised over all arguments.
free_protein_competitive <- function(p0, l0, i0, kd, ki, iters = 5L) {
  kd <- pmax(kd, 1e-300); ki <- pmax(ki, 1e-300)
  a <- kd + ki + l0 + i0 - p0
  b <- ki * (l0 - p0) + kd * (i0 - p0) + kd * ki
  c <- -kd * ki * p0
  inner <- pmax(a^2 - 3 * b, 0)
  arg   <- pmin(pmax((-2 * a^3 + 9 * a * b - 27 * c) / (2 * sqrt(inner^3)), -1), 1)
  P     <- -a / 3 + (2 / 3) * sqrt(inner) * cos(acos(arg) / 3)
  lb    <- p0 / (1 + l0 / kd + i0 / ki)
  P     <- ifelse(is.finite(P), pmin(pmax(P, lb), p0), lb)
  for (k in seq_len(iters)) {
    f  <- P + P * l0 / (kd + P) + P * i0 / (ki + P) - p0
    fp <- 1 + l0 * kd / (kd + P)^2 + i0 * ki / (ki + P)^2
    P  <- pmin(pmax(P - f / fp, lb), p0)
  }
  P
}

#' Equilibrium [PL] in the competitive system (Wang equation). Vectorised.
formula_competitive_pl <- function(p0, l0, i0, kd, ki) {
  P <- free_protein_competitive(p0, l0, i0, kd, ki)
  P * l0 / (kd + P)
}

#' Full equilibrium species: [P]free, [PL]eq, [PI]eq.
competitive_binding_species <- function(p0, l0, i0, kd, ki) {
  P <- free_protein_competitive(p0, l0, i0, kd, ki)
  list(p_eq = P, pl_eq = P * l0 / (kd + P), pi_eq = P * i0 / (ki + P))
}

#' [PL] without inhibitor: smaller root of the quadratic, cancellation-free.
pl_no_inhibitor <- function(p0, l0, kd) {
  S <- p0 + l0 + kd
  2 * p0 * l0 / (S + sqrt(pmax(S^2 - 4 * p0 * l0, 0)))
}

#' Exact Apparent IC50 (total [I]) and IC50 (free [I]) at 50% inhibition.
inhibition_calc <- function(p0, l0, kd, ki) {
  pl_eq_0  <- pl_no_inhibitor(p0, l0, kd)
  pl_eq_50 <- pl_eq_0 / 2
  p_free   <- pl_eq_50 * kd / (l0 - pl_eq_50)
  pi_eq    <- p0 - p_free - pl_eq_50
  app_IC50 <- pi_eq * (ki + p_free) / p_free
  IC50     <- app_IC50 - pi_eq
  list(pl_eq_0 = pl_eq_0, app_IC50 = app_IC50, IC50 = IC50)
}

#' Approximate Ki via the Cheng-Prusoff equation.
calc_ki_app1 <- function(l0, app_ic50, kd) app_ic50 / (1 + l0 / kd)

#' Exact Ki from Apparent IC50 (Wang's group equation, Nikolovska-Coleska et
#' al., Anal. Biochem. 2004). Requires [P]0. Vectorised over app_ic50.
calc_ki_exact <- function(p0, l0, app_ic50, kd) {
  pl_eq0 <- pl_no_inhibitor(p0, l0, kd)
  PL50   <- pl_eq0 / 2
  L50    <- l0 - PL50
  I50    <- app_ic50 - p0 + kd * PL50 / L50 + PL50
  I50 / (L50 / kd + (p0 - pl_eq0) / kd + 1)
}

#' 4-parameter logistic function for apparent IC50 fitting.
#'   a = lower asymptote (PL at high inhibitor, approaches 0)
#'   b = Hill slope
#'   c = apparent IC50
#'   d = upper asymptote (PL at zero inhibitor, approaches PL_eq0)
logistic_func <- function(i0, a, b, c, d) (a - d) / (1 + (i0 / c)^b) + d

#' Inhibitor series ([I]0, total, ascending).
#' "default":  log-spaced from apparent IC50/1000 (at least 1e-4 nM) to
#'             apparent IC50*1000
#' "dilution": top_conc / dil_factor^k, k = n_conc-1, ..., 0
compute_i0_concs <- function(p0, ki, l0, kd, n_conc, mode = "default",
                             top_conc = NA_real_, dil_factor = NA_real_) {
  if (identical(mode, "dilution")) return(top_conc / dil_factor^((n_conc - 1L):0L))
  app <- inhibition_calc(p0, l0, kd, ki)$app_IC50
  10^seq(log10(max(1e-4, app / 1000)), log10(app * 1000), length.out = n_conc)
}

#' IC50 curve data for plotting: log-spaced [I]0 around the apparent IC50,
#' extended to cover the simulated series when one is given.
compute_ic50_curve <- function(p0, l0, kd, ki, i0_series = NULL, n_points = 200L) {
  app <- inhibition_calc(p0, l0, kd, ki)$app_IC50
  if (!is.finite(app) || app <= 0) return(NULL)
  lo <- max(1e-4, app / 1000); hi <- app * 1000
  if (!is.null(i0_series)) { lo <- min(lo, min(i0_series)); hi <- max(hi, max(i0_series)) }
  i0_total <- 10^seq(log10(lo), log10(hi), length.out = n_points)
  sp <- competitive_binding_species(p0, l0, i0_total, kd, ki)
  list(i0_total = i0_total, i_free = i0_total - sp$pi_eq, pl = sp$pl_eq,
       pl_eq_0 = pl_no_inhibitor(p0, l0, kd))
}

# Per-experiment nlsLM fitting wrappers. Each returns c(est, conv):
# conv = 1 converged, 0 not converged, NA error.

#' Fit the 4-parameter logistic -> apparent IC50 (parameter c).
#' Starting values: a from the data at the highest [I]0, d from the data at
#' the lowest [I]0, b = -1, c = app_ic50_init (0.5 x the true apparent IC50).
fit_logistic_exp <- function(pl_mat, i0_concs, repN, app_ic50_init) {
  df   <- data.frame(pl = as.vector(pl_mat), i0 = rep(i0_concs, each = repN))
  para <- list(a = mean(pl_mat[, ncol(pl_mat)]), b = -1, c = app_ic50_init,
               d = max(mean(pl_mat[, 1L]), 0))
  fit <- NULL; retries <- 0L
  while (is.null(fit) && retries < 12L) {
    fit <- tryCatch(suppressWarnings(
      nlsLM(pl ~ logistic_func(i0, a, b, c, d), data = df, start = para,
            lower = c(-Inf, -Inf, 0, 0))),        # c (IC50) and d must be >= 0
      error = function(e) NULL)
    if (is.null(fit)) { para$a <- para$a / 2; retries <- retries + 1L }
  }
  if (is.null(fit)) return(c(est = NA_real_, conv = NA_real_))
  c(est = unname(coef(fit)["c"]), conv = as.numeric(isTRUE(fit$convInfo$isConv)))
}

#' Fit the Wang equation (formula_competitive_pl) -> Ki; [L]0 and Kd fixed,
#' [P]0 and Ki estimated.
fit_ki_exp <- function(pl_mat, i0_concs, repN, l0_fixed, kd_fixed, p0_init, ki_init) {
  df <- data.frame(pl = as.vector(pl_mat), i0 = rep(i0_concs, each = repN))
  tryCatch({
    fit <- suppressWarnings(
      nlsLM(pl ~ formula_competitive_pl(p0, l0 = l0_fixed, i0, kd = kd_fixed, ki),
            data = df, start = list(p0 = p0_init, ki = ki_init),
            lower = c(0, 0)))                     # p0 and ki must be >= 0
    c(est = unname(coef(fit)["ki"]), conv = as.numeric(isTRUE(fit$convInfo$isConv)))
  }, error = function(e) c(est = NA_real_, conv = NA_real_))
}

#' Outcome of each estimate: "ok", "boundary" (at the lower bound 0),
#' "nonconv" (not converged) or "failed" (error).
fit_status <- function(est, conv, truth) {
  ifelse(is.na(est) | is.na(conv), "failed",
    ifelse(conv == 0, "nonconv",
      ifelse(est <= BOUNDARY_TOL * truth, "boundary", "ok")))
}

#' KS distance between x and the normal distribution fitted to x
#' (mean and SD estimated from x; the Lilliefors statistic).
#' Identical to ks.test(x, "pnorm", mean(x), sd(x))$statistic. Written out because
#' only D is needed: the p-value that ks.test reports assumes a fully specified
#' null and is not valid here, and ks.test warns on tied values.
ks_stat_normal <- function(x) {
  n <- length(x)
  s <- sd(x)
  if (n < 2L || !is.finite(s) || s <= 0) return(NA_real_)
  z  <- sort((x - mean(x)) / s)
  Fz <- pnorm(z)
  i  <- seq_len(n)
  max(Fz - (i - 1) / n, i / n - Fz)
}

#' Parametric bootstrap KS test of normality. The null distribution depends
#' only on n, so bootstrap samples are drawn from N(0, 1) with the same n as the
#' data tested; one null distribution per n is cached in `null_cache`.
#' p = (1 + #{D* >= D}) / (B + 1).
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

#' Format a bootstrap p-value as "<i>p</i> = x" or "<i>p</i> < 0.005" (B = 200)
fmt_p <- function(p, B = KS_B) {
  if (length(p) != 1L || is.na(p)) return("<i>p</i> = NA")
  if (p <= 1 / (B + 1) + 1e-12) return(sprintf("<i>p</i> &lt; %.3f", 1 / (B + 1)))
  paste0("<i>p</i> = ", formatC(p, format = "f", digits = 3))
}

# Section 3 — Plot functions ----------------------------------------------

#' IC50 inhibition curve: [PL] vs [I]total (apparent IC50) and [PL] vs [I]free
#' (IC50) on a log x-axis, with the current experiment's simulated points.
make_ic50_plot <- function(curve_data, app_ic50, ic50_val, scatter_df = NULL) {
  if (is.null(curve_data)) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, size = 4, colour = "grey50",
                      label = "Inhibition curve not available for these inputs") +
             theme_void())
  }
  pl_half <- curve_data$pl_eq_0 / 2
  df <- rbind(
    data.frame(x = curve_data$i0_total, y = curve_data$pl, curve = "[I]total"),
    data.frame(x = curve_data$i_free,   y = curve_data$pl, curve = "[I]free"))
  df <- df[is.finite(df$x) & df$x > 0, ]

  p <- ggplot(df, aes(x = x, y = y, colour = curve)) +
    geom_line(linewidth = 1.0) +
    geom_hline(yintercept = pl_half, linetype = "dashed", colour = "grey50", alpha = 0.8) +
    geom_vline(xintercept = app_ic50, linetype = "dotdash", colour = "#fc8d62", alpha = 0.85) +
    geom_vline(xintercept = ic50_val, linetype = "dotdash", colour = "#66c2a5", alpha = 0.85) +
    scale_x_log10() +
    scale_colour_manual(values = c("[I]total" = "#fc8d62", "[I]free" = "#66c2a5"), name = NULL) +
    labs(x = "Concentration of Free/Total Inhibitor (nM)",
         y = "Concentration of Protein\u2013Ligand Complex ([PL], nM)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.95, 0.95),
          legend.justification = c(1, 1),
          legend.background    = element_rect(fill = alpha("white", 0.7), colour = NA),
          legend.key           = element_blank(),
          legend.text          = element_text(size = 13),
          legend.key.width     = unit(1.5, "cm"))

  if (!is.null(scatter_df) && nrow(scatter_df) > 0) {
    sc <- scatter_df[scatter_df$i0 > 0, ]
    if (nrow(sc) > 0)
      p <- p + geom_point(data = sc, aes(x = i0, y = pl),
                          colour = "#e6550d", alpha = 0.70, size = 2.4, inherit.aes = FALSE)
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

#' Normal Q-Q plot. Red line: the fitted normal (mean and SD of the plotted
#' values), the same reference as the KS test. Grey band: pointwise 95% range
#' of the i-th order statistic of a normal sample of the same size, from
#' U(i) ~ Beta(i, n - i + 1); computed directly, no random numbers used.
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
    theme(aspect.ratio = 1)
}

# Section 4 — UI ----------------------------------------------------------

ui <- fluidPage(
  title = "Binding Curve Viewer \u2014 Competitive Binding",
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),

  tags$head(tags$style(HTML("
    body { padding: 0 14px 20px 14px; }
    .p3-header {
      padding: 10px 2px 8px 2px;
      border-bottom: 2px solid #2c3e50;
      margin-bottom: 16px;
    }
    .p3-header h1 { display:inline; font-size:1.42rem; font-weight:650; margin:0; }
    .p3-header .home-link {
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
      background:#ffffff; font-size:13px; min-height:80px; margin-bottom:10px;
      line-height:1.6;
    }
    .res-warn { color:#a04000; }
    .hist-label {
      font-weight:650; color:#495057; margin:12px 0 4px 2px;
      font-size:0.82rem; text-transform:uppercase; letter-spacing:0.06em;
    }
    .dist-bar {
      display:flex; align-items:center; flex-wrap:wrap; gap:6px 22px;
      margin:12px 0 4px 0;
    }
    .dist-bar .hist-label { margin:0 0 0 2px; }
    .dist-bar .form-group, .dist-bar .shiny-input-container { margin-bottom:0; }
    .dist-note { font-size:12.5px; color:#555; margin:2px 0 0 2px; }
    .param-note { font-size:11.5px; color:#666; font-weight:400; }
    .conc-table td, .conc-table th { padding:1px 8px !important; }
  "))),

  # Page header
  div(class = "p3-header",
    h1(HTML("Competitive Binding &mdash; Inhibition Curve")),
    tags$a(class = "home-link", href = "/binding-curve-viewer/", "[Back to Home]")
  ),

  # Main row: config | inhibition curve + nav | binding parameters card.
  # Curve and card refresh only when Start Simulation or Clear Simulation is clicked.
  fluidRow(
    column(2,
      numericInput("p0_val",
        HTML("Initial Conc. of Fixed Component (<b>[P]<sub>0</sub></b>, nM)
              <br><span class='param-note'>0.001 &ndash; 50 000 nM;
              [P]<sub>0</sub>/<i>K</i><sub>i</sub> &le; 1000</span>"),
        value = DEFAULTS$p0, min = CONC_MIN, max = CONC_MAX, step = 1),
      numericInput("l0_val",
        HTML("Ligand Conc. (<b>[L]<sub>0</sub></b>, nM)
              <br><span class='param-note'>0.001 &ndash; 50 000 nM</span>"),
        value = DEFAULTS$l0, min = CONC_MIN, max = CONC_MAX, step = 1),
      numericInput("kd_val",
        HTML("Ligand's Equilibrium Dissociation Constant (<b><i>K</i><sub>d</sub></b>, nM)
              <br><span class='param-note'>0.001 &ndash; 50 000 nM</span>"),
        value = DEFAULTS$kd, min = CONC_MIN, max = CONC_MAX, step = 1),
      numericInput("ki_val",
        HTML("Inhibitor's Equilibrium Dissociation Constant (<b><i>K</i><sub>i</sub></b>, nM)
              <br><span class='param-note'>0.001 &ndash; 50 000 nM</span>"),
        value = DEFAULTS$ki, min = CONC_MIN, max = CONC_MAX, step = 1)
    ),

    column(6,
      plotOutput("ic50_fig", height = "420px"),
      div(class = "nav-bar",
        actionButton("nav_left",  HTML("&#9664;"), class = "btn btn-primary btn-sm"),
        uiOutput("nav_text"),
        actionButton("nav_right", HTML("&#9654;"), class = "btn btn-primary btn-sm")
      )
    ),

    column(4,
      bslib::card(
        height = 500,
        bslib::card_header(HTML("<b>Binding Parameters</b>")),
        bslib::card_body(uiOutput("params_card"))
      )
    )
  ),

  # Monte Carlo Simulation Console
  fluidRow(column(12,
    div(class = "sim-console",
      tags$h5("Monte Carlo Simulation",
              style = "margin:0 0 12px 0; font-weight:700; color:#2b579a;"),
      # Row 1: ExpN | n_conc | replicates | [PL] noise | seed
      fluidRow(
        column(2,
          selectInput("expn", "Number of Simulated Experiments",
            choices = c("1000" = "1000", "5000" = "5000", "10,000" = "10000",
                        "100,000" = "100000"), selected = "1000")
        ),
        column(2,
          selectInput("n_conc", "Number of Inhibitor Concentrations per Experiment",
            choices = as.character(5:10), selected = "9")
        ),
        column(2,
          selectInput("repn", "Replicates per Inhibitor Concentration",
            choices = as.character(3:10), selected = "5")
        ),
        column(2,
          selectInput("ratio_sd",
            HTML("[PL] Noise SD as Fraction of [P]<sub>0</sub>"),
            choices = c("0.0125", "0.0250", "0.0500"), selected = "0.0250")
        ),
        column(2,
          numericInput("seed_val",
            HTML("Random Seed <small style='color:#666;font-weight:400'>(positive integer)</small>"),
            value = 1, min = 1, step = 1)
        )
      ),
      # Row 2: inhibitor series | top [I]0 | dilution factor | buttons
      fluidRow(
        column(4,
          radioButtons("design_mode", "Inhibitor Concentration Series",
            choiceNames = list(
              HTML("Default: log-spaced, apparent IC50/1000 to apparent IC50&times;1000"),
              HTML("Serial dilution: top concentration and dilution factor")),
            choiceValues = list("default", "dilution"),
            selected = "default")
        ),
        column(2,
          conditionalPanel("input.design_mode == 'dilution'",
            numericInput("top_conc", HTML("Top [I]<sub>0</sub> (nM)"),
                         value = 10000, min = TOP_MIN, max = TOP_MAX))
        ),
        column(2,
          conditionalPanel("input.design_mode == 'dilution'",
            selectInput("dil_factor", "Dilution Factor",
              choices  = c("2" = "2", "3" = "3",
                           "3.16 (half-log)" = as.character(sqrt(10)), "10" = "10"),
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
      uiOutput("progress_ui")
    )
  )),

  # Result panels (logistic + mechanistic)
  fluidRow(
    column(6, div(class = "res-panel", uiOutput("res_logistic"))),
    column(6, div(class = "res-panel", uiOutput("res_mechanistic")))
  ),

  # Distributions (histograms or normal Q-Q plots), drawn square
  div(class = "dist-bar",
    div(class = "hist-label", "Fitted Parameter Distributions"),
    radioButtons("dist_plot", NULL, inline = TRUE,
      choices  = c("Histogram" = "hist", "Normal Q-Q plot" = "qq"),
      selected = "hist")
  ),
  fluidRow(
    column(6,
      div(class = "hist-label", HTML(
        "Logistic Model &mdash; Fitted Apparent IC50 Distributions")),
      fluidRow(
        column(6, plotOutput("hist_app_ic50",  height = "auto")),
        column(6, plotOutput("hist_lapp_ic50", height = "auto"))
      )
    ),
    column(6,
      div(class = "hist-label", HTML(
        "Mechanistic Model (Wang Equation) &mdash; Fitted
         <span style='text-transform:none'>Ki</span> Distributions")),
      fluidRow(
        column(6, plotOutput("hist_ki",  height = "auto")),
        column(6, plotOutput("hist_lki", height = "auto"))
      )
    )
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

  .init_ic <- inhibition_calc(DEFAULTS$p0, DEFAULTS$l0, DEFAULTS$kd, DEFAULTS$ki)

  sim <- reactiveValues(
    # Committed binding parameters — only updated on Start or Clear.
    committed_params = DEFAULTS,
    binding_params   = .init_ic,
    ki_cp            = calc_ki_app1(DEFAULTS$l0, .init_ic$app_IC50, DEFAULTS$kd),
    ki_ex            = calc_ki_exact(DEFAULTS$p0, DEFAULTS$l0, .init_ic$app_IC50, DEFAULTS$kd),
    curve_data       = compute_ic50_curve(DEFAULTS$p0, DEFAULTS$l0, DEFAULTS$kd, DEFAULTS$ki),

    scatter_data = NULL,          # data.frame(exp_idx, i0, pl), first <= 10 experiments

    # Estimates (all experiments) and outcomes
    app_ic50_vec = NULL, st_app = NULL,      # logistic fit
    ki_cp_vec    = NULL,                     # Cheng-Prusoff from each fitted apparent IC50
    ki_wg_vec    = NULL,                     # Wang's group equation from each fitted apparent IC50
    ki_vec       = NULL, st_ki = NULL,       # mechanistic fit

    ks_results = NULL,
    i0_concs   = NULL,

    current_idx = 1L,
    total_expn  = 0L,
    status_msg  = "Simulation not started."
  )

  clear_results <- function() {
    sim$scatter_data <- NULL
    sim$app_ic50_vec <- NULL; sim$st_app <- NULL
    sim$ki_cp_vec <- NULL; sim$ki_wg_vec <- NULL
    sim$ki_vec <- NULL; sim$st_ki <- NULL
    sim$ks_results <- NULL; sim$i0_concs <- NULL
    sim$total_expn <- 0L; sim$current_idx <- 1L
  }

  commit_params <- function(p0, l0, kd, ki, i0_series = NULL) {
    ic <- inhibition_calc(p0, l0, kd, ki)
    sim$committed_params <- list(p0 = p0, l0 = l0, kd = kd, ki = ki)
    sim$binding_params   <- ic
    sim$ki_cp            <- calc_ki_app1(l0, ic$app_IC50, kd)
    sim$ki_ex            <- calc_ki_exact(p0, l0, ic$app_IC50, kd)
    sim$curve_data       <- compute_ic50_curve(p0, l0, kd, ki, i0_series)
  }

  current_scatter <- reactive({
    if (is.null(sim$scatter_data)) return(NULL)
    sim$scatter_data[sim$scatter_data$exp_idx == sim$current_idx, ]
  })

  output$ic50_fig <- renderPlot({
    make_ic50_plot(sim$curve_data, sim$binding_params$app_IC50,
                   sim$binding_params$IC50, current_scatter())
  })

  # Binding Parameters card
  output$params_card <- renderUI({
    bp <- sim$binding_params; cp <- sim$committed_params
    cond_html <- paste0(
      "<p style='font-size:13px; line-height:1.6; margin:0 0 4px 0;'>",
      "Under the Condition of ",
      "[P]<sub>0</sub> = <b>", fmt(cp$p0), "</b> nM, ",
      "[L]<sub>0</sub> = <b>", fmt(cp$l0), "</b> nM, ",
      "<i>K</i><sub>d</sub> = <b>", fmt(cp$kd), "</b> nM, ",
      "and <i>K</i><sub>i</sub> = <b>", fmt(cp$ki), "</b> nM:<br>",
      "IC50 Equals <b>", fmt(bp$IC50, 5), "</b> nM;<br>",
      "Apparent IC50 Equals <b>", fmt(bp$app_IC50, 5), "</b> nM.",
      "</p>",
      "<hr style='margin:3px 0;'>",
      "<p style='font-size:13px; line-height:1.6; margin:0 0 4px 0;'>",
      "<b><i>K</i><sub>i</sub></b> from Apparent IC50 via Cheng&ndash;Prusoff Equation<sup>1</sup>:<br>",
      fmt(sim$ki_cp, 5), " nM<br><br>",
      "<b><i>K</i><sub>i</sub></b> from Apparent IC50 via Wang's Group Equation<sup>2</sup>:<br>",
      fmt(sim$ki_ex, 5), " nM",
      "</p>",
      "<small style='color:#555;'>",
      "<sup>1</sup> Cheng &amp; Prusoff, <i>Biochem. Pharmacol.</i> 1973, 22(23), 3099.<br>",
      "<sup>2</sup> Nikolovska-Coleska, et al., <i>Anal. Biochem.</i> 2004, 332(2), 261.<br>",
      "<sup>3</sup> DeLean, et al., <i>Am J Physiol.</i> 1978, 235(2), E97.<br>",
      "<sup>4</sup> Wang, <i>FEBS Lett.</i> 1995, 360(2), 111.<br></small>"
    )
    # Inhibitor series with the theoretical inhibition — shown after simulation
    conc_html <- if (!is.null(sim$i0_concs)) {
      i0  <- sort(unique(sim$i0_concs))
      inh <- 1 - formula_competitive_pl(cp$p0, cp$l0, i0, cp$kd, cp$ki) / bp$pl_eq_0
      informative <- inh >= 0.05 & inh <= 0.95
      rows <- paste0(
        "<tr", ifelse(informative, " style='font-weight:700;'", ""), ">",
        "<td>", vapply(i0, fmt, character(1)), "</td>",
        "<td>", sprintf("%.1f%%", 100 * inh), "</td></tr>", collapse = "")
      paste0(
        "<hr style='margin:3px 0;'>",
        "<span style='display:block; font-size:13px; font-weight:bold; margin:4px 0;'>",
        "Inhibitor Concentrations Used</span>",
        "<table class='table table-sm mb-1 conc-table' style='font-size:12.5px;'>",
        "<thead><tr><th>[I]<sub>0</sub> (nM)</th><th>Inhibition</th></tr></thead>",
        "<tbody>", rows, "</tbody></table>",
        "<small style='color:#555'>", sum(informative), " of ", length(i0),
        " concentrations lie between 5% and 95% inhibition (bold).</small>")
    } else ""
    HTML(paste0(cond_html, conc_html))
  })

  # Hint for the serial-dilution design (live, not committed)
  output$design_hint <- renderUI({
    ok <- function(v) length(v) == 1L && is.numeric(v) && is.finite(v) && v > 0
    p0 <- input$p0_val; l0 <- input$l0_val; kd <- input$kd_val; ki <- input$ki_val
    if (!ok(p0) || !ok(l0) || !ok(kd) || !ok(ki)) return(NULL)
    app <- inhibition_calc(p0, l0, kd, ki)$app_IC50
    top <- input$top_conc; f <- suppressWarnings(as.numeric(input$dil_factor))
    n   <- as.integer(input$n_conc)
    series_txt <- if (ok(top) && ok(f))
      sprintf(" Series: %s &rarr; %s nM (%d points, %s-fold steps).",
              fmt(top / f^(n - 1L)), fmt(top), n, fmt(f, 3)) else ""
    tags$p(style = "font-size:12.5px; color:#555; margin:0 0 4px 0;",
      HTML(sprintf("For the current inputs, apparent IC50 = %s nM.%s", fmt(app), series_txt)))
  })

  # Navigation
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
    if (max_idx > 0L) sim$current_idx <- (sim$current_idx - 2L) %% max_idx + 1L
  })
  observeEvent(input$nav_right, {
    max_idx <- min(10L, sim$total_expn)
    if (max_idx > 0L) sim$current_idx <- sim$current_idx %% max_idx + 1L
  })

  output$progress_ui <- renderUI({
    tags$p(HTML(sim$status_msg),
           style = "font-weight:bold; color:#1a5276; margin:10px 0 0 0; font-size:13px;")
  })

  # Result panels -------------------------------------------------------------

  ks_line <- function(label, r) {
    if (is.null(r) || is.na(r$D))
      return(paste0(label, ": not available",
                    if (!is.null(r$msg)) paste0(" (", r$msg, ")") else ""))
    paste0(label, ": <i>D</i> = ", fmt(r$D, 3), ", ", fmt_p(r$p), " (<i>n</i> = ", r$n, ")")
  }

  outcome_line <- function(status) {
    n_all <- length(status)
    cnt <- table(factor(status, levels = c("ok", "boundary", "nonconv", "failed")))
    if (cnt[["ok"]] == n_all)
      return(sprintf("Fit outcomes: all %d fits converged normally.<br>", n_all))
    sprintf(paste0("Fit outcomes (of %d): %d normal &middot; %d at the lower limit (0 nM) ",
                   "&middot; %d not converged &middot; %d errors<br>"),
            n_all, cnt[["ok"]], cnt[["boundary"]], cnt[["nonconv"]], cnt[["failed"]])
  }

  # Linear-scale summary: "ok" + "boundary" estimates
  stats_line <- function(label_html, vec, status) {
    v <- vec[status %in% c("ok", "boundary")]
    v <- v[is.finite(v)]
    if (length(v) < 2L) return(paste0(label_html, ": too few usable estimates<br>"))
    q <- quantile(v, c(0.025, 0.5, 0.975), names = FALSE)
    paste0(label_html, " (nM), <i>n</i> = ", length(v), ": ",
           "mean &plusmn; SD = ", fmt(mean(v)), " &plusmn; ", fmt(sd(v)),
           "; median = ", fmt(q[2]),
           " (2.5&ndash;97.5%: ", fmt(q[1]), "&ndash;", fmt(q[3]), ")<br>")
  }

  warn_line <- function(status, extra = character(0)) {
    w <- extra
    if (any(status == "boundary"))
      w <- c(w, sprintf(paste0("%s%% of fits are at the lower limit (0 nM): kept on the ",
                               "linear scale, excluded on log10"),
                        fmt(100 * mean(status == "boundary"), 3)))
    if (any(status %in% c("nonconv", "failed")))
      w <- c(w, sprintf(paste0("%d fits that did not converge and %d fits that ended in ",
                               "an error were excluded from all statistics"),
                        sum(status == "nonconv"), sum(status == "failed")))
    if (length(w)) paste0("<span class='res-warn'>&#9888; ", paste(w, collapse = "; "), ".</span>") else ""
  }

  ks_header <- paste0(
    "Normality test (Kolmogorov&ndash;Smirnov, mean and SD estimated from the data; ",
    "parametric bootstrap, <i>B</i> = ", KS_B, "):<br>")

  output$res_logistic <- renderUI({
    if (is.null(sim$app_ic50_vec) || is.null(sim$ks_results))
      return(HTML("<b>Logistic Model<sup>3</sup> Results:</b><br>Simulation not started."))
    st <- sim$st_app; ks <- sim$ks_results
    use <- st %in% c("ok", "boundary")
    n_neg_wg <- sum(use & is.finite(sim$ki_wg_vec) & sim$ki_wg_vec <= 0)
    extra <- if (n_neg_wg > 0)
      sprintf("%d <i>K</i><sub>i</sub> values from Wang's group equation are &le; 0", n_neg_wg) else character(0)
    HTML(paste0(
      "<b>Logistic Model<sup>3</sup> Results:</b><br>",
      outcome_line(st),
      stats_line("Fitted Apparent IC50", sim$app_ic50_vec, st),
      stats_line("<i>K</i><sub>i</sub> by Cheng&ndash;Prusoff Equation<sup>1</sup>", sim$ki_cp_vec, st),
      stats_line("<i>K</i><sub>i</sub> by Wang's Group Equation<sup>2</sup>", sim$ki_wg_vec, st),
      ks_header,
      "&nbsp;&nbsp;", ks_line("Apparent IC50", ks$app_lin), "<br>",
      "&nbsp;&nbsp;", ks_line("log10(Apparent IC50)", ks$app_log), "<br>",
      warn_line(st, extra)))
  })

  output$res_mechanistic <- renderUI({
    if (is.null(sim$ki_vec) || is.null(sim$ks_results))
      return(HTML("<b>Mechanistic Model<sup>4</sup> (Wang Equation) Results:</b><br>Simulation not started."))
    st <- sim$st_ki; ks <- sim$ks_results
    HTML(paste0(
      "<b>Mechanistic Model<sup>4</sup> (Wang Equation) Results:</b><br>",
      outcome_line(st),
      stats_line("Fitted <i>K</i><sub>i</sub>", sim$ki_vec, st),
      ks_header,
      "&nbsp;&nbsp;", ks_line("<i>K</i><sub>i</sub>", ks$ki_lin), "<br>",
      "&nbsp;&nbsp;", ks_line("log10(<i>K</i><sub>i</sub>)", ks$ki_log), "<br>",
      warn_line(st)))
  })

  # Distribution plots: linear scale "ok" + "boundary", log10 scale "ok" only
  make_dist_output <- function(out_id, vec_name, st_name, title, unit_lab, log_scale) {
    output[[out_id]] <- renderPlot({
      vec <- sim[[vec_name]]; st <- sim[[st_name]]
      if (is.null(vec)) return(NULL)
      if (log_scale) {
        v <- log10(vec[st == "ok"])
        ttl <- sprintf("%s: log10  (n = %d)", title, length(v))
        lab <- paste0("log10(", unit_lab, ")")
      } else {
        v <- vec[st %in% c("ok", "boundary")]
        ttl <- sprintf("%s  (n = %d)", title, length(v))
        lab <- unit_lab
      }
      if (identical(input$dist_plot, "qq")) make_qq_plot(v, ttl, lab)
      else make_histogram_plot(v, ttl, lab)
    },
    height = function() {
      w <- session$clientData[[paste0("output_", out_id, "_width")]]
      if (is.null(w) || !is.finite(w) || w <= 0) 300 else w
    })
  }
  make_dist_output("hist_app_ic50",  "app_ic50_vec", "st_app", "Logistic: Apparent IC50", "Apparent IC50 (nM)", FALSE)
  make_dist_output("hist_lapp_ic50", "app_ic50_vec", "st_app", "Logistic: Apparent IC50", "Apparent IC50 (nM)", TRUE)
  make_dist_output("hist_ki",        "ki_vec",       "st_ki",  "Mechanistic: Ki",         "Ki (nM)",            FALSE)
  make_dist_output("hist_lki",       "ki_vec",       "st_ki",  "Mechanistic: Ki",         "Ki (nM)",            TRUE)

  # Input validation; returns NULL when all inputs are valid --------------------
  check_inputs <- function(need_seed = TRUE) {
    is_num <- function(v) length(v) == 1L && is.numeric(v) && is.finite(v)
    in_rng <- function(v, lo, hi) is_num(v) && v >= lo && v <= hi
    msgs <- character(0)
    labs <- c(p0_val = "[P]<sub>0</sub>", l0_val = "[L]<sub>0</sub>",
              kd_val = "<i>K</i><sub>d</sub>", ki_val = "<i>K</i><sub>i</sub>")
    for (id in names(labs))
      if (!in_rng(input[[id]], CONC_MIN, CONC_MAX))
        msgs <- c(msgs, paste0(labs[[id]], " must be between 0.001 and 50 000 nM."))
    if (need_seed) {
      s <- input$seed_val
      if (!is_num(s) || s < 1 || s > .Machine$integer.max || s != round(s))
        msgs <- c(msgs, "Random seed must be a positive integer.")
    }
    if (identical(input$design_mode, "dilution") &&
        !in_rng(input$top_conc, TOP_MIN, TOP_MAX))
      msgs <- c(msgs, "Top [I]<sub>0</sub> must be between 0.001 and 1e9 nM.")
    if (length(msgs)) paste(msgs, collapse = "<br>") else NULL
  }

  # Start Simulation ----------------------------------------------------------
  observeEvent(input$start_sim, {
    expN <- as.integer(input$expn)

    if (!IS_LOCAL && expN > 1000L) {
      showModal(modalDialog(
        title = tags$span(style = "color:#c0392b; font-weight:700;",
                          HTML("&#9888;&#65039; Large Simulation &mdash; Server Load Warning")),
        HTML(paste0(
          "Running <b>Number of Simulated Experiments &gt; 1000</b> on the web server ",
          "may cause heavy load and slow the application for all users.<br><br>",
          "Please <b>download and run the app locally</b>:<br>",
          "<a href='https://github.com/ydu-sci/Binding_Curve_Viewer' target='_blank' style='font-weight:600;'>",
          "github.com/ydu-sci/Binding_Curve_Viewer</a><br>(app file: <code>src/shiny-apps/comp-sim/app.R</code>)<br><br>",
          "Once downloaded, launch it from that folder in <b>local mode</b>:<br>",
          "<code style='background:#f1f1f1;padding:3px 6px;border-radius:3px;'>",
          "BINDING_APP_LOCAL=true Rscript -e \"shiny::runApp('app.R')\"</code>",
          "<br><small style='color:#555;'>R console: ",
          "<code>Sys.setenv(BINDING_APP_LOCAL=\"true\"); shiny::runApp(\"app.R\")</code></small>"
        )),
        footer = modalButton("OK"), easyClose = TRUE
      ))
      return()
    }

    err <- check_inputs()
    if (!is.null(err)) {
      showNotification(HTML(err), type = "error", duration = 8)
      return()
    }

    p0 <- as.numeric(input$p0_val); l0 <- as.numeric(input$l0_val)
    kd <- as.numeric(input$kd_val); ki <- as.numeric(input$ki_val)
    repN   <- as.integer(input$repn)
    n_conc <- as.integer(input$n_conc)
    ratio  <- as.numeric(input$ratio_sd)
    seed   <- as.integer(input$seed_val)
    mode   <- input$design_mode
    top    <- if (identical(mode, "dilution")) as.numeric(input$top_conc) else NA_real_
    dil    <- if (identical(mode, "dilution")) as.numeric(input$dil_factor) else NA_real_

    # [P]0/Ki limits: > RATIO_MAX refused, > RATIO_WARN warned (tight binding)
    ratio_pk <- p0 / ki
    if (ratio_pk > RATIO_MAX) {
      showModal(modalDialog(
        title = tags$span(style = "color:#c0392b; font-weight:700;",
                          HTML("[P]<sub>0</sub> / <i>K</i><sub>i</sub> too high")),
        HTML(sprintf(paste0(
          "[P]<sub>0</sub> / <i>K</i><sub>i</sub> = %s exceeds %d.<br><br>",
          "Under such tight-binding conditions the inhibition curve reflects ",
          "[P]<sub>0</sub> rather than <i>K</i><sub>i</sub>, and <i>K</i><sub>i</sub> ",
          "cannot be estimated. Please lower [P]<sub>0</sub> or raise <i>K</i><sub>i</sub>."),
          fmt(ratio_pk, 3), RATIO_MAX)),
        footer = modalButton("OK"), easyClose = TRUE))
      return()
    }
    ratio_note <- NULL
    if (ratio_pk > RATIO_WARN) {
      ratio_note <- sprintf(paste0(
        "[P]<sub>0</sub> / <i>K</i><sub>i</sub> = %s (&gt; %d): tight binding. ",
        "Fitted <i>K</i><sub>i</sub> values are imprecise and some fits may reach ",
        "<i>K</i><sub>i</sub> = 0; consider lowering [P]<sub>0</sub> below <i>K</i><sub>i</sub>."),
        fmt(ratio_pk, 3), RATIO_WARN)
      showNotification(HTML(ratio_note), type = "warning", duration = 12)
    }

    ic <- inhibition_calc(p0, l0, kd, ki)
    i0_concs <- compute_i0_concs(p0, ki, l0, kd, n_conc, mode, top, dil)
    commit_params(p0, l0, kd, ki, i0_concs)
    clear_results()
    sim$status_msg <- sprintf("Running (ExpN=%d, seed=%d) ...", expN, seed)

    tryCatch(
      withProgress(
        message = sprintf("Monte Carlo Simulation (ExpN = %d)", expN),
        value = 0,
        {
          set.seed(seed)

          # Noise-free [PL] at each [I]0 (identical in every experiment) plus
          # additive [PL] readout noise — the only random draws, in the same
          # order as v4 with zero concentration errors.
          n_total  <- repN * n_conc * expN
          pl_theo  <- formula_competitive_pl(p0, l0, rep(i0_concs, each = repN), kd, ki)
          noisy_array <- array(rep(pl_theo, expN) + rnorm(n_total, sd = ratio * p0),
                               dim = c(repN, n_conc, expN))

          app_est <- app_conv <- numeric(expN)
          ki_est  <- ki_conv  <- numeric(expN)
          store_idx <- integer(0); store_i0 <- numeric(0); store_pl <- numeric(0)
          app_init <- ic$app_IC50 * 0.5
          p0_init  <- p0 * 0.5
          ki_init  <- ki * 0.5
          upd_every <- max(1L, expN %/% 100L)

          for (i in seq_len(expN)) {
            pl_mat_i <- noisy_array[, , i]
            if (i <= 10L) {
              store_idx <- c(store_idx, rep(i, n_conc * repN))
              store_i0  <- c(store_i0,  rep(i0_concs, each = repN))
              store_pl  <- c(store_pl,  as.vector(pl_mat_i))
            }
            r1 <- fit_logistic_exp(pl_mat_i, i0_concs, repN, app_init)
            app_est[i] <- r1[["est"]]; app_conv[i] <- r1[["conv"]]
            r2 <- fit_ki_exp(pl_mat_i, i0_concs, repN, l0_fixed = l0, kd_fixed = kd,
                             p0_init = p0_init, ki_init = ki_init)
            ki_est[i] <- r2[["est"]]; ki_conv[i] <- r2[["conv"]]

            if (i %% upd_every == 0L || i == expN)
              incProgress(0.80 * upd_every / expN,
                          detail = sprintf("Fitting experiment %d / %d  (%.0f%%)",
                                           i, expN, 100 * i / expN))
          }

          st_app <- fit_status(app_est, app_conv, ic$app_IC50)
          st_ki  <- fit_status(ki_est,  ki_conv,  ki)

          # Ki from each experiment's fitted apparent IC50
          ki_cp_vec <- calc_ki_app1(l0, app_est, kd)
          ki_wg_vec <- calc_ki_exact(p0, l0, app_est, kd)

          # Bootstrap KS tests (progress 80% -> 100%)
          null_cache <- new.env(parent = emptyenv())
          ks_specs <- list(
            app_lin = list(app_est[st_app %in% c("ok", "boundary")], "Apparent IC50"),
            app_log = list(log10(app_est[st_app == "ok"]),          "log10(Apparent IC50)"),
            ki_lin  = list(ki_est[st_ki %in% c("ok", "boundary")],  "Ki"),
            ki_log  = list(log10(ki_est[st_ki == "ok"]),            "log10(Ki)"))
          ks_res <- list()
          for (k in names(ks_specs)) {
            incProgress(0.05, detail = paste0("KS test: ", ks_specs[[k]][[2]], " ..."))
            x <- ks_specs[[k]][[1]]
            ks_res[[k]] <- tryCatch(
              ks_boot_test(x, null_cache),
              error = function(e) list(n = length(x), D = NA_real_, p = NA_real_,
                                       msg = conditionMessage(e)))
          }

          sim$scatter_data <- data.frame(exp_idx = store_idx, i0 = store_i0, pl = store_pl)
          sim$app_ic50_vec <- app_est; sim$st_app <- st_app
          sim$ki_cp_vec    <- ki_cp_vec
          sim$ki_wg_vec    <- ki_wg_vec
          sim$ki_vec       <- ki_est;  sim$st_ki  <- st_ki
          sim$ks_results   <- ks_res
          sim$i0_concs     <- i0_concs
          sim$total_expn   <- expN
          sim$current_idx  <- 1L
          sim$status_msg   <- paste0(
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

  # Clear Simulation ----------------------------------------------------------
  # Clears all results and redraws the theoretical curve from the current inputs.
  observeEvent(input$clear_sim, {
    err <- check_inputs(need_seed = FALSE)
    clear_results()
    if (!is.null(err)) {
      showNotification(HTML(err), type = "error", duration = 8)
      sim$status_msg <- "Simulation cleared (curve not updated: invalid inputs)."
      return()
    }
    commit_params(as.numeric(input$p0_val), as.numeric(input$l0_val),
                  as.numeric(input$kd_val), as.numeric(input$ki_val))
    sim$status_msg <- "Simulation cleared."
  })
}

# Section 6 — Launch ------------------------------------------------------

shinyApp(ui = ui, server = server)
