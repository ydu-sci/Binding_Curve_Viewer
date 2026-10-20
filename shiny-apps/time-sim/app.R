# Time Dependency of Dissociation Constant (Kd)
#
# This app demonstrates how the Kd obtained from a saturation binding
# experiment depends on the incubation time allowed before measurement.
# No experimental noise is considered — all [PL] values are theoretical.
#
# Three curves are shown:
#   1. Theoretical equilibrium [PL]        (formula1, t -> infinity)
#   2. [PL] reached at incubation time t1  (exact second-order kinetics)
#   3. [PL] reached at incubation time t2  (exact second-order kinetics)
#
# The t1 and t2 data are then fitted with formula1 (the exact equilibrium
# binding equation) by Levenberg-Marquardt least squares (nlsLM), with
# [P]total and Kd both free and constrained to be >= 0. The resulting
# apparent Kd values are reported together with the outcome of each fit.
# When the incubation time is too short, the apparent Kd is overestimated.
#
# LIGAND CONCENTRATION DESIGN — each time point gets its OWN [L]total series.
# A single shared series is not used. Instead, for each incubation time the
# series is solved so that [PL] AT THAT TIME spans 1% to 99% of [P]total
# (20 log-spaced points), so both saturation curves are sampled over the same
# occupancy range and are equally well determined. The choice hardly changes
# the apparent Kd: with the defaults, fitting the 2 h data on the equilibrium
# 1-99% series gives 13.37x the true Kd instead of 13.36x. On that shared
# series the 2 h curve still reaches 98.9% occupancy at the top (high [L]total
# equilibrates fast); the shortfall is at the bottom, where the slow
# low-[L]total points have reached only 0.1% occupancy instead of 1%.
#
# INCUBATION-TIME GUIDANCE — the card reports the incubation time after which
# every [L]total has reached 95% and 99% of its equilibrium [PL]. For each
# [L]total the time is obtained from the closed-form inverse of the exact
# second-order solution; the reported value is its maximum over [L]total.
#
# UNITS: kon is entered in M^-1 s^-1 and converted internally to nM^-1 s^-1
# (multiply by 1e-9), so that Kd = koff / kon is expressed in nM. With the
# defaults below, Kd = 1e-5 / (6.1e5 * 1e-9) = 0.0163934 nM.
#
# Requirements:
#   install.packages(c("shiny", "bslib", "ggplot2", "minpack.lm"))
#   bslib >= 0.4.0, ggplot2 >= 3.5.0
#
# Run:
#   Rscript -e "shiny::runApp('app.R')"

library(shiny)
library(bslib)
library(ggplot2)
library(minpack.lm)

# Fractional occupancy limits spanned by each concentration series
OCC_LO <- 0.01     # 1%  of [P]total
OCC_HI <- 0.99     # 99% of [P]total
N_CONC <- 20L      # concentrations per series

# Allowed input ranges (checked on the server before every run)
P0_MIN   <- 0.001;  P0_MAX   <- 50000    # nM
KOFF_MIN <- 1e-6;   KOFF_MAX <- 1        # s^-1
KON_MIN  <- 1e3;    KON_MAX  <- 1e8      # M^-1 s^-1
T_MIN    <- 0.001;  T_MAX    <- 1e5      # hours

RATIO_WARN   <- 10      # [P]total/Kd above this: warning (titration regime)
RATIO_MAX    <- 1000    # [P]total/Kd above this: run is refused
BOUNDARY_TOL <- 1e-6    # fitted Kd <= BOUNDARY_TOL * true Kd counts as the Kd = 0 bound
START_RATIO  <- 0.5     # nlsLM starting values = START_RATIO * true [P]total and Kd

# Default parameters (the default state is computed from these at app load)
DEFAULTS <- list(
  p0    = 0.008,    # nM
  koff  = 1e-5,     # s^-1
  kon_M = 6.1e5,    # M^-1 s^-1 (user-facing unit)
  t1_h  = 2,        # hours
  t2_h  = 20        # hours
)

# Section 1 — Mathematical helper functions

#' Exact equilibrium binding equation (Formula 1).
#' Smaller root of [PL]^2 - S*[PL] + [P]t*[L]t = 0, S = [P]t + [L]t + Kd.
#' Written as 2*[P]t*[L]t / (S + sqrt(D)) rather than (S - sqrt(D)) / 2: the two
#' are algebraically identical, but the second subtracts two nearly equal
#' numbers when [P]t << [L]t + Kd and loses precision (catastrophic
#' cancellation). The discriminant D = S^2 - 4[P]t[L]t is evaluated as
#' ([P]t - [L]t)^2 + Kd*(2([P]t + [L]t) + Kd), which is the same quantity
#' without the subtraction of two large numbers.
#' @param x   Total ligand concentration [L]total (nM), vectorised
#' @param p0  Total fixed component concentration [P]total (nM)
#' @param kd  Dissociation constant (nM)
formula1 <- function(x, p0, kd) {
  S <- p0 + x + kd
  D <- (p0 - x)^2 + kd * (2 * (p0 + x) + kd)
  2 * p0 * x / (S + sqrt(pmax(D, 0)))
}

#' Roots of the binding quadratic, evaluated without cancellation.
#' pl_2 is the larger root, pl_1 = [PL]eq the smaller one (pl_1 * pl_2 = [P]t*[L]t),
#' and rD = pl_2 - pl_1 = sqrt(D) > 0 for any Kd > 0.
quad_roots <- function(l0, p0, kd) {
  S    <- p0 + l0 + kd
  rD   <- sqrt(pmax((p0 - l0)^2 + kd * (2 * (p0 + l0) + kd), 0))
  pl_2 <- (S + rD) / 2
  list(rD = rD, pl_2 = pl_2, pl_1 = p0 * l0 / pl_2)
}

#' Exact second-order association kinetics: [PL] at time t.
#' Same solution as formula_second_order_kin() in the reference R script,
#'     [PL](t) = pl_1*pl_2*(1 - E) / (pl_2 - pl_1*E),  E = exp(-kon*rD*t),
#' rearranged with pl_2 - pl_1*E = rD + pl_1*(1 - E) and 1 - E = -expm1(-kon*rD*t)
#' so that no step subtracts nearly equal numbers. As t -> infinity this
#' converges to formula1(l0, p0, koff/kon).
#' @param t     Incubation time (s), scalar
#' @param l0    Total ligand concentration (nM), vectorised
#' @param p0    Total fixed component concentration (nM)
#' @param koff  Off-rate constant (s^-1)
#' @param kon   On-rate constant (nM^-1 s^-1, already converted from M^-1 s^-1)
formula_second_order_kin <- function(t, l0, p0, koff, kon) {
  q   <- quad_roots(l0, p0, koff / kon)
  omE <- -expm1(-kon * q$rD * t)                # 1 - E, in [0, 1]
  p0 * l0 * omE / (q$rD + q$pl_1 * omE)
}

#' Time (s) at which [PL] reaches the fraction frac of its equilibrium value.
#' Closed-form inverse of the fraction of equilibrium reached at time t,
#'     [PL](t) / [PL]eq = pl_2 * (1 - E) / (rD + pl_1 * (1 - E)),  E = exp(-kon*rD*t):
#'     t = log(1 + frac*rD / ((1 - frac)*pl_2)) / (kon*rD)
time_to_fraction <- function(frac, l0, p0, koff, kon) {
  q <- quad_roots(l0, p0, koff / kon)
  log1p(frac * q$rD / ((1 - frac) * q$pl_2)) / (kon * q$rD)
}

#' Incubation time after which every [L]total has reached 95% and 99% of its
#' equilibrium [PL]. The maximum of time_to_fraction() over [L]total is located
#' on a log grid spanning 8 decades beyond [P]total and Kd on either side, then
#' refined with optimize() between the neighbouring grid points.
equilibration_summary <- function(p0, koff, kon) {
  kd <- koff / kon
  lg <- seq(log10(min(p0, kd)) - 8, log10(max(p0, kd)) + 8, length.out = 2001L)
  longest <- function(frac) {
    f <- function(x) time_to_fraction(frac, 10^x, p0, koff, kon)
    v <- f(lg)
    i <- which.max(v)
    o <- optimize(f, lg[c(max(1L, i - 1L), min(length(lg), i + 1L))],
                  maximum = TRUE, tol = 1e-10)
    max(v[i], o$objective)
  }
  list(t95 = longest(0.95), t99 = longest(0.99))
}

#' Total ligand concentration giving a target fractional occupancy AT EQUILIBRIUM.
#'
#' Inverting formula1 at [PL] = f * [P]total gives
#'     [L]total = ([P]total + Kd - [PL]) * [PL] / ([P]total - [PL])
#'              = f * [P]total + Kd * f / (1 - f)
#'
#' Because [PL] at any finite time is at most its equilibrium value, this
#' closed-form result is also a LOWER BOUND for the time-resolved solution
#' below, which makes it a safe left bracket for the root search.
l0_at_occupancy_eq <- function(f, p0, kd) {
  f * p0 + kd * f / (1 - f)
}

#' Total ligand concentration giving a target fractional occupancy AT TIME t.
#'
#' Solves  [PL](t, [L]total) = f * [P]total  for [L]total.
#' [PL](t, .) is strictly increasing in [L]total and approaches [P]total as
#' [L]total -> infinity for any t > 0, so a unique root always exists and
#' bisection is safe. The search runs on a log10 scale.
#'
#' @param f     Target fractional occupancy, strictly between 0 and 1
#' @param p0    Total fixed component concentration (nM)
#' @param koff  Off-rate constant (s^-1)
#' @param kon   On-rate constant (nM^-1 s^-1)
#' @param t     Incubation time (s)
#' @return      [L]total in nM, or NA if the target cannot be bracketed
l0_at_occupancy_time <- function(f, p0, koff, kon, t) {
  target <- f * p0
  kd     <- koff / kon
  fn <- function(lg) formula_second_order_kin(t, 10^lg, p0, koff, kon) - target

  # Left bracket: the equilibrium solution (smallest possible [L]total)
  lg_lo <- log10(l0_at_occupancy_eq(f, p0, kd))

  # If the reaction is already equilibrated at this concentration, the
  # equilibrium bracket IS the answer. This case is common at high occupancy
  # and long t. A plain fn(lg_lo) >= 0 test is not sufficient: when fully
  # equilibrated, fn(lg_lo) is zero only to within floating-point noise and
  # can come back very slightly negative, which would send the expansion loop
  # below on a fruitless 60-step search ending in NA. Hence the tolerance.
  tol <- target * 1e-9
  if (fn(lg_lo) >= -tol) return(10^lg_lo)

  # Expand the right bracket upward until the target is exceeded
  lg_hi <- lg_lo
  k <- 0L
  while (fn(lg_hi) < 0 && k < 60L) {
    lg_hi <- lg_hi + 0.5
    k <- k + 1L
  }
  if (fn(lg_hi) < 0) return(NA_real_)          # target unreachable

  10^uniroot(fn, c(lg_lo, lg_hi), tol = 1e-12)$root
}

#' Adaptive ligand concentration series for one incubation time.
#'
#' Spans OCC_LO to OCC_HI fractional occupancy AT TIME t, log-spaced. The top
#' concentration therefore always drives [PL] to 99% of [P]total at that
#' specific time, whatever [P]total, Kd and t the user enters.
#'
#' @param n  Number of concentration points
make_l0_series_at_time <- function(p0, koff, kon, t, n = N_CONC) {
  lo <- l0_at_occupancy_time(OCC_LO, p0, koff, kon, t)
  hi <- l0_at_occupancy_time(OCC_HI, p0, koff, kon, t)
  if (is.na(lo) || is.na(hi) || hi <= lo) return(rep(NA_real_, n))
  10^seq(log10(lo), log10(hi), length.out = n)
}

#' Fit formula1 to a set of ([L]total, [PL]) points to recover the apparent Kd.
#' [P]total is also a free parameter of the fit (as in a real experiment), but
#' only the Kd estimate is reported. Levenberg-Marquardt (nlsLM) with both
#' parameters constrained to be >= 0, started at START_RATIO x the true values.
#' (Base nls() is not used: its convergence test fails on zero-residual data,
#' which is exactly what a noise-free fit at long incubation times produces.)
#' Returns c(p0, kd, conv): conv = 1 converged, 0 not converged, NA failed.
fit_kd <- function(pl_obs, l0_list, p0_true, kd_true) {
  na <- c(p0 = NA_real_, kd = NA_real_, conv = NA_real_)
  if (any(!is.finite(l0_list)) || any(!is.finite(pl_obs))) return(na)
  df <- data.frame(pl = pl_obs, x = l0_list)
  tryCatch({
    fit <- suppressWarnings(nlsLM(
      pl ~ formula1(x, p0, kd),
      data    = df,
      start   = list(p0 = p0_true * START_RATIO, kd = kd_true * START_RATIO),
      lower   = c(p0 = 0, kd = 0),
      control = nls.lm.control(maxiter = 200)))
    cf <- coef(fit)
    c(p0 = unname(cf[["p0"]]), kd = unname(cf[["kd"]]),
      conv = as.numeric(isTRUE(fit$convInfo$isConv)))
  }, error = function(e) na)
}

#' Classify a fit: "ok", "boundary" (Kd at the lower bound 0),
#' "nonconv" (iteration limit reached) or "failed" (no fit returned).
fit_status <- function(fit, kd_true) {
  if (is.na(fit[["kd"]]) || is.na(fit[["conv"]])) return("failed")
  if (fit[["conv"]] == 0) return("nonconv")
  if (fit[["kd"]] <= BOUNDARY_TOL * kd_true) return("boundary")
  "ok"
}

#' Compute the full simulation state from a set of input parameters.
#' Called once at app load for the defaults, and on every "Start Simulation".
#' @param kon_M  On-rate constant as entered by the user, in M^-1 s^-1
compute_state <- function(p0, koff, kon_M, t1_h, t2_h) {
  kon <- kon_M * 1e-9                          # M^-1 s^-1 -> nM^-1 s^-1
  kd  <- koff / kon
  t1  <- t1_h * 3600                           # hours -> seconds
  t2  <- t2_h * 3600

  # Independent concentration series, one per incubation time
  l0_t1 <- make_l0_series_at_time(p0, koff, kon, t1)
  l0_t2 <- make_l0_series_at_time(p0, koff, kon, t2)

  pl_t1 <- formula_second_order_kin(t1, l0_t1, p0, koff, kon)
  pl_t2 <- formula_second_order_kin(t2, l0_t2, p0, koff, kon)

  fit_t1 <- fit_kd(pl_t1, l0_t1, p0, kd)
  fit_t2 <- fit_kd(pl_t2, l0_t2, p0, kd)

  list(
    p0 = p0, koff = koff, kon_M = kon_M, kon = kon, kd = kd,
    t1_h = t1_h, t2_h = t2_h,
    l0_t1 = l0_t1, pl_t1 = pl_t1,
    l0_t2 = l0_t2, pl_t2 = pl_t2,
    fit_t1 = fit_t1, st_t1 = fit_status(fit_t1, kd),
    fit_t2 = fit_t2, st_t2 = fit_status(fit_t2, kd),
    eq = equilibration_summary(p0, koff, kon)
  )
}

# Section 2 — Default state
#
# Computed once when the app is loaded (about 0.05 s), so the app shows a
# complete result immediately and Reset restores exactly this state.
# With the defaults: Kd = 0.0163934 nM; series 0.002413 -> 1.718 nM (2 h) and
# 0.0003731 -> 1.631 nM (20 h); apparent Kd 13.36x (2 h) and 1.46x (20 h) the
# true value; every [L]total reaches 95% / 99% of equilibrium after 55.9 / 86.0 h.

DEFAULT_STATE <- compute_state(DEFAULTS$p0, DEFAULTS$koff, DEFAULTS$kon_M,
                               DEFAULTS$t1_h, DEFAULTS$t2_h)

# Section 3 — Formatting helpers and plot function

#' Format a number with 6 significant digits (Kd values, parameters)
fmt6 <- function(v) {
  if (length(v) != 1L || !is.finite(v)) return("\u2014")
  formatC(v, format = "g", digits = 6)
}

#' Format a positive number with 3 significant digits, keeping trailing zeros
#' (86.0 rather than 86)
sig3 <- function(x) {
  if (!is.finite(x)) return("\u2014")
  if (x >= 100) return(formatC(round(x), format = "d", big.mark = ""))
  formatC(x, format = "f", digits = max(0L, 2L - floor(log10(x))))
}

#' Format a duration given in seconds
fmt_time <- function(sec) {
  h <- sec / 3600
  if (!is.finite(h)) return("\u2014")
  if (sec < 60) return(paste(sig3(sec), "s"))
  if (h < 1) return(paste(sig3(h * 60), "min"))
  if (h < 48) return(paste(sig3(h), "h"))
  sprintf("%s h (%s d)", sig3(h), sig3(h / 24))
}

#' Build the three-curve saturation binding plot.
#'
#' The two time curves use different concentration series, so the x range spans
#' the union of both. Curves are evaluated densely from the stored parameters
#' (a pure analytical evaluation, not a re-fit); points mark the 20
#' concentrations actually used for each non-linear fit.
make_binding_plot <- function(st) {

  # Series keys are plain ASCII. The legend shows t1 and t2 with plotmath
  # subscripts, which work with any font; Unicode subscript digits are missing
  # from some server fonts and were drawn as empty boxes.
  lab_eq <- "Theoretical equilibrium"
  lab_t1 <- sprintf("t1 = %g h", st$t1_h)
  lab_t2 <- sprintf("t2 = %g h", st$t2_h)
  lvls   <- c(lab_eq, lab_t1, lab_t2)
  leg_labels <- as.expression(list(
    lab_eq,
    bquote(t[1] == .(sprintf("%g h", st$t1_h))),
    bquote(t[2] == .(sprintf("%g h", st$t2_h)))
  ))

  # x range covers the union of both concentration series
  all_l0 <- c(st$l0_t1, st$l0_t2)
  all_l0 <- all_l0[is.finite(all_l0) & all_l0 > 0]
  if (length(all_l0) == 0L) return(ggplot() + theme_void())

  l0_dense <- 10^seq(log10(min(all_l0)), log10(max(all_l0)), length.out = 400L)
  t1_s <- st$t1_h * 3600
  t2_s <- st$t2_h * 3600

  df_curve <- rbind(
    data.frame(l0 = l0_dense,
               pl = formula1(l0_dense, st$p0, st$kd),
               series = lab_eq),
    data.frame(l0 = l0_dense,
               pl = formula_second_order_kin(t1_s, l0_dense, st$p0, st$koff, st$kon),
               series = lab_t1),
    data.frame(l0 = l0_dense,
               pl = formula_second_order_kin(t2_s, l0_dense, st$p0, st$koff, st$kon),
               series = lab_t2)
  )
  df_curve$series <- factor(df_curve$series, levels = lvls)

  # Each time point contributes its own 20 measured concentrations
  df_pts <- rbind(
    data.frame(l0 = st$l0_t1, pl = st$pl_t1, series = lab_t1),
    data.frame(l0 = st$l0_t2, pl = st$pl_t2, series = lab_t2)
  )
  df_pts$series <- factor(df_pts$series, levels = lvls)

  pal <- setNames(c("black", "#fc8d62", "#66c2a5"), lvls)
  lty <- setNames(c("dashed", "solid", "solid"),    lvls)
  shp <- setNames(c(NA, 16, 17),                    lvls)

  ggplot() +
    geom_line(data = df_curve,
              aes(x = l0, y = pl, colour = series, linetype = series),
              linewidth = 0.9, na.rm = TRUE) +
    geom_point(data = df_pts,
               aes(x = l0, y = pl, colour = series, shape = series),
               size = 2.3, na.rm = TRUE) +
    scale_x_log10() +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE,
                        breaks = lvls, labels = leg_labels) +
    scale_linetype_manual(values = lty, name = NULL, drop = FALSE,
                          breaks = lvls, labels = leg_labels) +
    scale_shape_manual(values = shp, name = NULL, drop = FALSE,
                       breaks = lvls, labels = leg_labels,
                       na.translate = FALSE) +
    labs(x = "Concentration of Total Ligand (nM)",
         y = "Concentration of Bound Protein ([PL], nM)") +
    theme_bw(base_size = 13) +
    theme(legend.position        = "inside",
          legend.position.inside = c(0.02, 0.98),
          legend.justification   = c(0, 1),
          legend.background      = element_rect(fill = alpha("white", 0.75),
                                                colour = NA),
          legend.key             = element_blank(),
          legend.text            = element_text(size = 12, hjust = 0),
          legend.key.width       = unit(1.4, "cm"),
          axis.title             = element_text(size = 16),
          axis.title.x           = element_text(margin = margin(t = 8)),
          axis.title.y           = element_text(margin = margin(r = 8)),
          axis.text              = element_text(size = 12))
}

# Section 4 — UI

ui <- fluidPage(
  title = "Binding Curve Viewer \u2014 Time Dependency of Kd",
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),

  tags$head(tags$style(HTML("
    body { padding: 0 14px 20px 14px; }
    .p4-header {
      padding: 10px 2px 8px 2px;
      border-bottom: 2px solid #2c3e50;
      margin-bottom: 16px;
    }
    .p4-header h1 { display:inline; font-size:1.42rem; font-weight:650; margin:0; }
    .p4-header .home-link {
      margin-left:14px; font-size:0.95rem; font-weight:500; white-space:nowrap;
    }
    .param-note { font-size:11.5px; color:#666; font-weight:400; }
    .live-kd {
      font-size:12px; color:#2b579a; font-weight:600;
      margin:-6px 0 12px 2px;
    }
    .res-block { font-size:13.5px; line-height:1.75; }
    .res-block p { margin:0 0 6px 0; }
    .res-warn { color:#a04000; font-size:12px; font-weight:600; }
    .conc-list {
      font-size:11px; color:#333; display:block;
      margin:1px 0 0 0; line-height:1.35; word-break:break-all;
    }
  "))),

  # Page header
  div(class = "p4-header",
    h1(HTML("Time Dependency of Dissociation Constant (<i>K</i><sub>d</sub>)")),
    tags$a(class = "home-link", href = "/binding-curve-viewer/", "[Back to Home]")
  ),

  fluidRow(

    # config column
    column(3,
      numericInput("p0_val",
        HTML("Total Fixed Component Conc. (<b>[P]<sub>total</sub></b>, nM)
              <br><span class='param-note'>0.001 \u2013 50 000 nM</span>"),
        value = DEFAULTS$p0, min = P0_MIN, max = P0_MAX, step = 0.001),

      numericInput("koff_val",
        HTML("Off-rate Constant (<b><i>k</i><sub>off</sub></b>, s<sup>-1</sup>)
              <br><span class='param-note'>1e-6 \u2013 1 s<sup>-1</sup></span>"),
        value = DEFAULTS$koff, min = KOFF_MIN, max = KOFF_MAX, step = 1e-6),

      numericInput("kon_val",
        HTML("On-rate Constant (<b><i>k</i><sub>on</sub></b>,
              M<sup>-1</sup>s<sup>-1</sup>)
              <br><span class='param-note'>1e3 \u2013 1e8
              M<sup>-1</sup>s<sup>-1</sup></span>"),
        value = DEFAULTS$kon_M, min = KON_MIN, max = KON_MAX, step = 1e3),

      # Live feedback: Kd and [P]total/Kd implied by the current inputs.
      # Updates while typing; the committed values used for the plot and the
      # fits are the ones shown in the results card on the right.
      uiOutput("live_kd"),

      numericInput("t1_val",
        HTML("Incubation Time <b>t<sub>1</sub></b> (hours)
              <br><span class='param-note'>0.001 \u2013 100 000 h</span>"),
        value = DEFAULTS$t1_h, min = T_MIN, max = T_MAX, step = 0.5),

      numericInput("t2_val",
        HTML("Incubation Time <b>t<sub>2</sub></b> (hours)
              <br><span class='param-note'>0.001 \u2013 100 000 h</span>"),
        value = DEFAULTS$t2_h, min = T_MIN, max = T_MAX, step = 1),

      tags$div(style = "margin-top:14px;",
        actionButton("start_sim", "Start Simulation",
          class = "btn btn-success btn-sm", style = "margin-right:6px;"),
        actionButton("reset_sim", "Reset Simulation",
          class = "btn btn-danger btn-sm")
      ),

      uiOutput("status_ui")
    ),

    # plot
    column(5,
      plotOutput("binding_plot", height = "500px")
    ),

    # results card
    column(4,
      bslib::card(
        height = 560,
        bslib::card_header(HTML("<b>Binding Parameters</b>")),
        bslib::card_body(
          style = "overflow-y:auto;",
          uiOutput("results_card")
        )
      )
    )
  ),

  tags$br()
)

# Section 5 — Server

#' Validate the inputs; returns NULL if all are valid, otherwise an HTML
#' string listing every problem.
check_inputs <- function(p0, koff, kon_M, t1_h, t2_h) {
  is_num <- function(v) length(v) == 1L && is.numeric(v) && is.finite(v)
  in_rng <- function(v, lo, hi) is_num(v) && v >= lo && v <= hi
  msgs <- character(0)
  if (!in_rng(p0, P0_MIN, P0_MAX))
    msgs <- c(msgs, "[P]<sub>total</sub> must be between 0.001 and 50 000 nM.")
  if (!in_rng(koff, KOFF_MIN, KOFF_MAX))
    msgs <- c(msgs, "<i>k</i><sub>off</sub> must be between 1e-6 and 1 s<sup>-1</sup>.")
  if (!in_rng(kon_M, KON_MIN, KON_MAX))
    msgs <- c(msgs, "<i>k</i><sub>on</sub> must be between 1e3 and 1e8 M<sup>-1</sup>s<sup>-1</sup>.")
  if (!in_rng(t1_h, T_MIN, T_MAX))
    msgs <- c(msgs, "Incubation time t<sub>1</sub> must be between 0.001 and 100 000 h.")
  if (!in_rng(t2_h, T_MIN, T_MAX))
    msgs <- c(msgs, "Incubation time t<sub>2</sub> must be between 0.001 and 100 000 h.")
  if (length(msgs)) paste(msgs, collapse = "<br>") else NULL
}

server <- function(input, output, session) {

  # Committed state: the default state computed at app load.
  # Only "Start Simulation" replaces it with a freshly computed state.
  sim <- reactiveValues(
    st         = DEFAULT_STATE,
    status_msg = "Showing default parameters."
  )

  # ── Live Kd hint (updates while typing; not the committed value) ─────────
  output$live_kd <- renderUI({
    ok <- is.numeric(input$koff_val) && is.numeric(input$kon_val) &&
          is.finite(input$koff_val) && is.finite(input$kon_val) &&
          input$kon_val > 0
    if (!ok) return(div(class = "live-kd", "Current input: \u2014"))
    # Convert kon from M^-1 s^-1 to nM^-1 s^-1 before dividing
    kd_live <- input$koff_val / (input$kon_val * 1e-9)
    ratio_txt <- ""
    if (is.numeric(input$p0_val) && is.finite(input$p0_val) && kd_live > 0)
      ratio_txt <- sprintf("<br>[P]<sub>total</sub> / <i>K</i><sub>d</sub> = %s",
                           formatC(input$p0_val / kd_live, format = "g", digits = 4))
    div(class = "live-kd",
        HTML(sprintf("Current input: <i>K</i><sub>d</sub> = %s nM%s",
                     formatC(kd_live, format = "g", digits = 6), ratio_txt)))
  })

  # ── Status message ───────────────────────────────────────────────────────
  output$status_ui <- renderUI({
    tags$p(HTML(sim$status_msg),
           style = "font-weight:bold; color:#1a5276; margin:10px 0 0 0;
                    font-size:12.5px;")
  })

  # ── Plot ─────────────────────────────────────────────────────────────────
  output$binding_plot <- renderPlot({
    make_binding_plot(sim$st)
  })

  # ── Results card ─────────────────────────────────────────────────────────
  output$results_card <- renderUI({
    st <- sim$st
    eq <- st$eq

    # Compact rendering of one concentration series
    conc_str <- function(v) {
      if (all(is.na(v))) return("\u2014")
      paste(formatC(sort(v), format = "g", digits = 4), collapse = ", ")
    }

    # One fitted-Kd paragraph, with the fit outcome
    fit_block <- function(fit, status, lab, t_h) {
      kd_hat <- fit[["kd"]]
      head <- paste0("<p><b>Fitted <i>K</i><sub>d</sub> (", lab, " = ", t_h,
                     " h)</b> = <b>",
                     if (status == "failed") "\u2014" else fmt6(kd_hat), "</b> nM<br>")
      note <- switch(status,
        ok = paste0("<span class='param-note'>", sprintf("%.2f", kd_hat / st$kd),
                    "\u00d7 the true <i>K</i><sub>d</sub></span>"),
        boundary = paste0("<span class='res-warn'>Fit reached the lower limit ",
                          "<i>K</i><sub>d</sub> = 0 nM; not a usable estimate.</span>"),
        nonconv = paste0("<span class='param-note'>", sprintf("%.2f", kd_hat / st$kd),
                         "\u00d7 the true <i>K</i><sub>d</sub></span><br>",
                         "<span class='res-warn'>Fit did not converge ",
                         "(iteration limit reached); value not reliable.</span>"),
        failed = if (all(is.na(if (lab == "t<sub>1</sub>") st$l0_t1 else st$l0_t2)))
                   "<span class='res-warn'>Concentration series could not be computed; no fit.</span>"
                 else
                   "<span class='res-warn'>Fit failed (no solution returned).</span>")
      paste0(head, note, "</p>")
    }

    HTML(paste0(
      "<div class='res-block'>",

      # True kinetic Kd
      "<p><b>True <i>K</i><sub>d</sub></b> ",
      "(<i>k</i><sub>off</sub> / <i>k</i><sub>on</sub>) = <b>",
      fmt6(st$kd), "</b> nM</p>",

      "<hr style='margin:8px 0;'>",

      # Fits at t1 and t2
      fit_block(st$fit_t1, st$st_t1, "t<sub>1</sub>", st$t1_h),
      fit_block(st$fit_t2, st$st_t2, "t<sub>2</sub>", st$t2_h),

      "<hr style='margin:8px 0;'>",

      # Incubation-time guidance
      "<p><b>Incubation time for every [L]<sub>total</sub> to reach 95% / 99% ",
      "of its equilibrium [PL]</b>: <b>", fmt_time(eq$t95), "</b> / <b>",
      fmt_time(eq$t99), "</b></p>",

      "<hr style='margin:8px 0;'>",

      # Committed input parameters
      "<p><span style='font-weight:bold; font-size:12px;'>",
      "Parameters used:</span><br>",
      "<span class='param-note'>",
      "[P]<sub>total</sub> = ", fmt6(st$p0), " nM &nbsp;|&nbsp; ",
      "<i>k</i><sub>off</sub> = ", fmt6(st$koff), " s<sup>-1</sup><br>",
      "<i>k</i><sub>on</sub> = ", fmt6(st$kon_M), " M<sup>-1</sup>s<sup>-1</sup>",
      " &nbsp;|&nbsp; [P]<sub>total</sub> / <i>K</i><sub>d</sub> = ",
      formatC(st$p0 / st$kd, format = "g", digits = 4),
      "</span></p>",

      "<hr style='margin:8px 0;'>",

      # Two independent concentration series, one per incubation time
      "<p><span style='font-weight:bold; font-size:12px;'>",
      "[L]<sub>total</sub> series for t<sub>1</sub> = ", st$t1_h,
      " h (nM), 1\u201399% of [P]<sub>total</sub> at t<sub>1</sub>:</span>",
      "<span class='conc-list'>", conc_str(st$l0_t1), "</span></p>",

      "<p><span style='font-weight:bold; font-size:12px;'>",
      "[L]<sub>total</sub> series for t<sub>2</sub> = ", st$t2_h,
      " h (nM), 1\u201399% of [P]<sub>total</sub> at t<sub>2</sub>:</span>",
      "<span class='conc-list'>", conc_str(st$l0_t2), "</span></p>",

      "</div>"
    ))
  })

  # ── Start Simulation: recompute everything from the current inputs ───────
  observeEvent(input$start_sim, {

    # Validate inputs
    err <- check_inputs(input$p0_val, input$koff_val, input$kon_val,
                        input$t1_val, input$t2_val)
    if (!is.null(err)) {
      showNotification(HTML(err), type = "error", duration = 8)
      sim$status_msg <- "Inputs out of range; simulation not run."
      return()
    }

    p0    <- as.numeric(input$p0_val)
    koff  <- as.numeric(input$koff_val)
    kon_M <- as.numeric(input$kon_val)
    kd    <- koff / (kon_M * 1e-9)

    # [P]total/Kd limits: > RATIO_MAX refused, > RATIO_WARN warned
    ratio_pk <- p0 / kd
    if (ratio_pk > RATIO_MAX) {
      showModal(modalDialog(
        title = tags$span(style = "color:#c0392b; font-weight:700;",
                          HTML("[P]<sub>total</sub> / <i>K</i><sub>d</sub> too high")),
        HTML(sprintf(paste0(
          "[P]<sub>total</sub> / <i>K</i><sub>d</sub> = %s exceeds %d.<br><br>",
          "Under such strong titration conditions the binding curve reflects ",
          "[P]<sub>total</sub> rather than <i>K</i><sub>d</sub>, and <i>K</i><sub>d</sub> ",
          "cannot be estimated from a saturation experiment. Please lower ",
          "[P]<sub>total</sub>, or change <i>k</i><sub>off</sub> / <i>k</i><sub>on</sub> ",
          "to raise <i>K</i><sub>d</sub>."),
          formatC(ratio_pk, format = "g", digits = 4), RATIO_MAX)),
        footer = modalButton("OK"),
        easyClose = TRUE
      ))
      sim$status_msg <- "[P]<sub>total</sub> / <i>K</i><sub>d</sub> too high; simulation not run."
      return()
    }
    ratio_note <- NULL
    if (ratio_pk > RATIO_WARN) {
      ratio_note <- sprintf(paste0(
        "[P]<sub>total</sub> / <i>K</i><sub>d</sub> = %s (&gt; %d): titration regime. ",
        "The fitted <i>K</i><sub>d</sub> is poorly determined; consider lowering ",
        "[P]<sub>total</sub> below <i>K</i><sub>d</sub>."),
        formatC(ratio_pk, format = "g", digits = 4), RATIO_WARN)
      showNotification(HTML(ratio_note), type = "warning", duration = 12)
    }

    res <- tryCatch(
      compute_state(
        p0    = p0,
        koff  = koff,
        kon_M = kon_M,                         # M^-1 s^-1; converted inside
        t1_h  = as.numeric(input$t1_val),
        t2_h  = as.numeric(input$t2_val)
      ),
      error = function(e) e)

    if (inherits(res, "error")) {
      showNotification(paste("Simulation error:", conditionMessage(res)),
                       type = "error", duration = 10)
      sim$status_msg <- "Simulation failed; previous result kept."
      return()
    }

    sim$st <- res
    msg <- "Simulation updated."
    if (!all(c(res$st_t1, res$st_t2) == "ok"))
      msg <- paste(msg, "See the fit warnings in the results card.")
    if (!is.null(ratio_note)) msg <- paste0(msg, "<br><span style='color:#a04000;'>",
                                            ratio_note, "</span>")
    sim$status_msg <- msg
  })

  # ── Reset Simulation: restore the default state and inputs ───────────────
  observeEvent(input$reset_sim, {
    updateNumericInput(session, "p0_val",   value = DEFAULTS$p0)
    updateNumericInput(session, "koff_val", value = DEFAULTS$koff)
    updateNumericInput(session, "kon_val",  value = DEFAULTS$kon_M)
    updateNumericInput(session, "t1_val",   value = DEFAULTS$t1_h)
    updateNumericInput(session, "t2_val",   value = DEFAULTS$t2_h)

    sim$st         <- DEFAULT_STATE
    sim$status_msg <- "Reset to default parameters."
  })
}

# Section 6 — Launch
shinyApp(ui = ui, server = server)
