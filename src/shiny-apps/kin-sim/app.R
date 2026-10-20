#
# Kinetics of Association and Dissociation — v3
#
# Pure R Shiny rewrite of kinetics_association_dissociation() from
# binding_curve_viewer.py. No Plotly / client-side JS interactivity is used;
# all computation happens server-side in R.
#
# Mathematical core (association/dissociation theoretical curves) is ported
# from binding_curve_viewer.py.
#
# Monte Carlo simulation design (concentration-point generation, time-point
# generation, and the three Gauss-Newton fitting models) is ported from
# Binding_Para_MC_Simulation_kin_upgrade.R, "Part 2. Binding Kinetics
# Experiment" (l0_calc, second_order_asso_points, gen_plNoisy_array_kin,
# formula_binding_kinetics_disso_group, formula_second_order_kin_group,
# formula_binding_kinetics_group).
#
# Requirements:
#   install.packages(c("shiny", "bslib", "ggplot2"))
#   bslib >= 0.4.0, ggplot2 >= 3.5.0
#
# Run modes:
#   Server mode (default) — Number of Simulated Experiments > 1000 shows a
#   "run locally" warning modal instead of starting the simulation:
#     Rscript -e "shiny::runApp('app.R')"
#
#   Local mode — Number of Simulated Experiments of 5000, 10000 or 100000 allowed:
#     Shell:     BINDING_APP_LOCAL=true Rscript -e "shiny::runApp('app.R')"
#     R console: Sys.setenv(BINDING_APP_LOCAL = "true"); shiny::runApp("app.R")
#
# Changes in v3 (relative to v2):
#   1. Normality test: KSgeneral::cont_ks_test against a normal whose mean and
#      SD were estimated from the same data (invalid p-values) is replaced by a
#      parametric bootstrap KS test with B = 200. The bootstrap sample size
#      equals the number of values actually tested. D and n are reported with p.
#      KSgeneral is no longer required.
#   2. The exact 2nd-order solution uses a cancellation-free form
#      (pl_2 = (S + sqrt(D))/2, pl_1 = P0*L0/pl_2, pl_1 - pl_2 = -sqrt(D)),
#      algebraically identical to v2. In v2 the 2nd-order fit failed when
#      Kd >> [P]0 (e.g. [P]0 = 0.001 nM, Kd = 1e6 nM).
#   3. Fit outcomes are recorded (normal, estimate <= 0, not converged, error).
#      With nls(warnOnly = TRUE), v2 silently kept fits that had not converged.
#      Not-converged and error fits are excluded from all statistics; estimates
#      <= 0 are kept on the linear scale and excluded on log10. kon is fitted
#      only when the koff fit of the same experiment is usable.
#   4. Every test is wrapped in tryCatch and the whole simulation is guarded,
#      so an error is reported in the app instead of disconnecting the session.
#   5. The single-[L]0 input and its theoretical curves are removed. The plots
#      show the true curves of the [L]0 series used in the simulation, and the
#      card lists per-[L]0 quantities, including [L]0/[P]0. A note appears when
#      the pseudo-first-order condition [L]0 >= 10*[P]0 is not met.
#   6. Default [P]0 = 10 nM. Time points keep the v2 design (equal steps of
#      [PL]: association 0-99.9%, dissociation 100-0.1% of [PL]eq, per [L]0),
#      which is assumed to sample every curve sufficiently.
#   7. Inputs are validated on the server; Number of Simulated Experiments
#      1000, 5000, 10000 or 100000; values above 1000 run in local mode only.
#      All noise is drawn before fitting (fixed RNG order), so 100000 experiments
#      with the largest design (10 [L]0 x 15 time points x 5 replicates) need about
#      1.2 GB of memory; the default design needs about 130 MB.
#   8. Results: fit outcomes, median and 2.5-97.5% range, significant-digit
#      formatting; distributions as histograms or normal Q-Q plots, drawn square.
#   9. [L]0 shown with an ordered (light to dark) colour scale and labelled by
#      concentration; "Back to Home" link in the header.
#  10. The exact 2nd-order solution evaluates 1 - E as -expm1(-kon*rD*t). With
#      1 - exp(...) the relative error reached 3.5e-8 at very short times; it is
#      now below 1e-14 over the whole input range (checked against 60-digit
#      reference values for 3000 random parameter sets).
#

library(shiny)
library(bslib)
library(ggplot2)

# Run-mode flag
IS_LOCAL <- identical(Sys.getenv("BINDING_APP_LOCAL"), "true")

# Limits and settings
P0_MIN   <- 0.001;  P0_MAX   <- 50000   # nM
KOFF_MIN <- 1e-6;   KOFF_MAX <- 1       # s^-1
KON_MIN  <- 1e3;    KON_MAX  <- 1e8     # M^-1 s^-1
PSEUDO_RATIO <- 10   # pseudo-first-order condition: [L]0 >= 10 * [P]0
KS_B     <- 200L     # bootstrap replicates for the KS test

# Defaults used on page load (before any simulation)
DEFAULTS <- list(p0 = 10, koff = 0.01, kon_M = 1e5, n_groups = 8L, time_pts = 10L)

# Section 1 — Mathematical helper functions -------------------------------
#
# Units: concentrations in nM, koff in s^-1, kon in nM^-1 s^-1 internally
# (the user enters kon in M^-1 s^-1; kon_nM = kon_M * 1e-9).

#' Roots of [PL]^2 - S*[PL] + P0*L0 = 0, S = P0 + L0 + Kd, without cancellation.
#' pl_1 = [PL]eq (smaller root), pl_2 = larger root, rD = sqrt(S^2 - 4*P0*L0)
#' = pl_2 - pl_1. pl_1 is computed as P0*L0/pl_2 instead of (S - rD)/2, which
#' loses precision when P0*L0 << S^2 (e.g. Kd >> [P]0).
kin_roots <- function(p0, l0, kd) {
  S    <- p0 + l0 + kd
  rD   <- sqrt(pmax(S^2 - 4 * p0 * l0, 0))
  pl_2 <- (S + rD) / 2
  list(pl_1 = p0 * l0 / pl_2, pl_2 = pl_2, rD = rD)
}

#' Exact 2nd-order association model, [PL](t) for [PL](0) = 0.
#' Algebraically identical to pl_1*pl_2*(1 - E)/(pl_2 - pl_1*E) with
#' E = exp(kon*(pl_1 - pl_2)*t), written with the cancellation-free roots and
#' rearranged with pl_2 - pl_1*E = rD + pl_1*(1 - E) and 1 - E = -expm1(-kon*rD*t)
#' (as in the time-dependency app), so that no step subtracts nearly equal
#' numbers when kon*rD*t is small.
formula_second_order_kin_group <- function(t, l0, p0, koff, kon) {
  S    <- p0 + l0 + koff / kon
  rD   <- sqrt(pmax(S^2 - 4 * p0 * l0, 0))
  pl_1 <- p0 * l0 / ((S + rD) / 2)
  omE  <- -expm1(-kon * rD * t)                 # 1 - E, in [0, 1]
  p0 * l0 * omE / (rD + pl_1 * omE)
}

#' Time at which the exact 2nd-order association reaches a fraction `frac`
#' of [PL]eq: t = ln((pl_2 - frac*pl_1) / (pl_2*(1 - frac))) / (kon*rD)
assoc_time_exact <- function(frac, p0, l0, koff, kon) {
  r <- kin_roots(p0, l0, koff / kon)
  log((r$pl_2 - frac * r$pl_1) / (r$pl_2 * (1 - frac))) / (kon * r$rD)
}

#' Dissociation model: [PL](t) = pl0[g] * exp(-koff * t)
#' (dissociation without rebinding, starting from the equilibrium complex)
formula_binding_kinetics_disso_group <- function(t, pl0, g, koff) {
  pl0[g] * exp(-koff * t)
}

#' Pseudo-1st-order association model: [PL](t) = pl_eq[g]*(1-exp(-(kon*l0+koff)*t))
formula_binding_kinetics_group <- function(t, l0, pl_eq, g, koff, kon) {
  pl_eq[g] * (1 - exp(-(kon * l0 + koff) * t))
}

#' Generalised l0_calc(): n ligand concentrations ([L]0, total) giving 20%-90%
#' fractional occupancy of [P]0 at equilibrium. n = 8 reproduces the original
#' l0_calc() exactly (seq(0.2, 0.9, 0.1) has 8 elements).
l0_calc_general <- function(p0, kd, n) {
  points <- seq(0.2, 0.9, length.out = n)
  pl <- points * p0
  (p0 + kd - pl) * pl / (p0 - pl)
}

#' Experimental design: [L]0 series, time points and noise-free [PL] for every
#' [L]0, per-[L]0 summary quantities, and smooth true curves for the plots.
#' Time points (as in v2): association times at equal steps of [PL] (0 to 99.9%
#' of [PL]eq) and dissociation times at equal steps of [PL] (100% to 0.1%),
#' both set per [L]0 from the true kinetics. This design is assumed to sample
#' every curve sufficiently.
build_kin_design <- function(p0, koff, kon_M, n_groups, time_pts) {
  kon <- kon_M * 1e-9
  kd  <- koff / kon
  l0  <- l0_calc_general(p0, kd, n_groups)
  pl1 <- kin_roots(p0, l0, kd)$pl_1

  asso <- vector("list", n_groups); disso <- vector("list", n_groups)
  for (j in seq_len(n_groups)) {
    fa <- seq(0, 0.999, length.out = time_pts)
    ta <- assoc_time_exact(fa, p0, l0[j], koff, kon)
    pa <- pl1[j] * fa
    fd <- rev(seq(0.001, 1, length.out = time_pts))
    td <- -log(fd) / koff
    pd <- pl1[j] * fd
    asso[[j]]  <- data.frame(g = j, t = ta, pl = pa)
    disso[[j]] <- data.frame(g = j, t = td, pl = pd)
  }
  asso <- do.call(rbind, asso); disso <- do.call(rbind, disso)

  # Smooth true curves for the plots
  grid_a <- seq(0, max(asso$t),  length.out = 300L)
  grid_d <- seq(0, max(disso$t), length.out = 300L)
  curve_a <- do.call(rbind, lapply(seq_len(n_groups), function(j) {
    pseudo_eq <- p0 * l0[j] / (l0[j] + kd)
    kobs <- kon * l0[j] + koff
    rbind(
      data.frame(g = j, t = grid_a, model = "Exact 2nd-order",
                 pl = formula_second_order_kin_group(grid_a, l0[j], p0, koff, kon)),
      data.frame(g = j, t = grid_a, model = "Pseudo-1st-order",
                 pl = pseudo_eq * (1 - exp(-kobs * grid_a))))
  }))
  curve_d <- do.call(rbind, lapply(seq_len(n_groups), function(j)
    data.frame(g = j, t = grid_d, pl = pl1[j] * exp(-koff * grid_d))))

  t99        <- assoc_time_exact(0.99, p0, l0, koff, kon)
  pseudo_t99 <- log(100) / (kon * l0 + koff)

  list(p0 = p0, koff = koff, kon = kon, kon_M = kon_M, kd = kd,
       n_groups = n_groups, time_pts = time_pts,
       l0 = l0, pl1 = pl1, asso = asso, disso = disso,
       curve_a = curve_a, curve_d = curve_d,
       t99 = t99, pseudo_t99 = pseudo_t99, pt99_ratio = pseudo_t99 / t99,
       disso_t50 = log(2) / koff, disso_t99 = log(100) / koff)
}

# Per-experiment nls() wrappers (Gauss-Newton). nls() handles the per-group
# parameters (pl0[g], pl_eq[g]). Each returns c(est, conv):
# conv = 1 converged, 0 not converged (nls warnOnly), NA error.

nls_est <- function(fit_fun, coef_name) {
  tryCatch({
    fit <- suppressWarnings(fit_fun())
    c(est = unname(coef(fit)[coef_name]), conv = as.numeric(isTRUE(fit$convInfo$isConv)))
  }, error = function(e) c(est = NA_real_, conv = NA_real_))
}

#' koff from one experiment's dissociation data (all [L]0 groups together)
fit_disso_exp <- function(t, pl, g, pl0_init, koff_init) {
  df <- data.frame(t = t, pl = pl, g = g)
  nls_est(function()
    nls(pl ~ formula_binding_kinetics_disso_group(t, pl0, g, koff),
        data = df, start = list(pl0 = pl0_init, koff = koff_init),
        control = nls.control(maxiter = 50L, warnOnly = TRUE)), "koff")
}

#' kon (nM^-1 s^-1) from association data with the EXACT 2nd-order model,
#' koff fixed to the same experiment's dissociation estimate; p0 is fitted as
#' one shared scalar.
fit_sec_order_exp <- function(t, l0, g, pl, p0_init, koff_fixed, kon_init) {
  df <- data.frame(t = t, l0 = l0, pl = pl, g = g)
  nls_est(function()
    nls(pl ~ formula_second_order_kin_group(t, l0, p0, koff = koff_fixed, kon),
        data = df, start = list(p0 = p0_init, kon = kon_init),
        control = nls.control(maxiter = 50L, warnOnly = TRUE)), "kon")
}

#' kon (nM^-1 s^-1) from the SAME association data with the PSEUDO-1ST-ORDER
#' model, koff fixed as above; one plateau pl_eq per [L]0 group.
fit_pseudo_exp <- function(t, l0, g, pl, pleq_init, koff_fixed, kon_init) {
  df <- data.frame(t = t, l0 = l0, pl = pl, g = g)
  nls_est(function()
    nls(pl ~ formula_binding_kinetics_group(t, l0, pl_eq, g, koff = koff_fixed, kon),
        data = df, start = list(pl_eq = pleq_init, kon = kon_init),
        control = nls.control(maxiter = 50L, warnOnly = TRUE)), "kon")
}

#' Outcome of each estimate: "ok", "nonpos" (estimate <= 0), "nonconv" (not
#' converged) or "failed" (error, or not fitted).
est_status <- function(est, conv) {
  ifelse(is.na(est) | is.na(conv), "failed",
    ifelse(conv == 0, "nonconv",
      ifelse(est <= 0, "nonpos", "ok")))
}

#' Outcome of Kd = koff/kon from the outcomes of its two fits
kd_status <- function(st_off, st_on, kd) {
  ifelse(st_off == "failed" | st_on == "failed" | is.na(kd), "failed",
    ifelse(st_off == "nonconv" | st_on == "nonconv", "nonconv",
      ifelse(kd <= 0, "nonpos", "ok")))
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

#' Parametric bootstrap KS test of normality.
#' The null distribution of the statistic depends only on n (normal is a
#' location-scale family), so bootstrap samples are drawn from N(0, 1) with the
#' same n as the data tested; one null distribution per n is cached in
#' `null_cache` (an environment). p = (1 + #{D* >= D}) / (B + 1).
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

#' Ordered colour scale for the [L]0 series: light (lowest [L]0) to dark (highest)
get_group_colors <- function(n) {
  # the 3 lightest steps are dropped so the lowest [L]0 stays visible on white
  rev(grDevices::hcl.colors(n + 3L, palette = "Blues 3")[seq_len(n)])
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

l0_labels <- function(d) vapply(d$l0, fmt, character(1), digits = 3)

#' Association plot: true curves of every [L]0 (solid: exact 2nd-order;
#' dashed: pseudo-1st-order), plus the current experiment's simulated points.
make_kin_on_plot <- function(d, scatter_df = NULL) {
  cols <- setNames(get_group_colors(d$n_groups), seq_len(d$n_groups))
  p <- ggplot() +
    geom_line(data = d$curve_a,
              aes(x = t, y = pl, colour = factor(g), linetype = model),
              linewidth = 0.8) +
    scale_linetype_manual(values = c("Exact 2nd-order" = "solid",
                                     "Pseudo-1st-order" = "22"), name = NULL)
  if (!is.null(scatter_df) && nrow(scatter_df) > 0)
    p <- p + geom_point(data = scatter_df, aes(x = t, y = pl, colour = factor(g)),
                        size = 2.2, alpha = 0.85)
  p + scale_colour_manual(values = cols, labels = l0_labels(d),
                          name = "[L]0 (nM)") +
    labs(title = "Association", x = "Time (s)",
         y = "Concentration of Bound Protein ([PL], nM)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "right", legend.text = element_text(size = 11),
          legend.key.width = unit(1.1, "cm"),
          plot.title = element_text(size = 13, face = "bold"))
}

#' Dissociation plot: true exponential decay of every [L]0's equilibrium
#' complex, plus the current experiment's simulated points.
make_kin_disso_plot <- function(d, scatter_df = NULL) {
  cols <- setNames(get_group_colors(d$n_groups), seq_len(d$n_groups))
  p <- ggplot() +
    geom_line(data = d$curve_d, aes(x = t, y = pl, colour = factor(g)),
              linewidth = 0.8)
  if (!is.null(scatter_df) && nrow(scatter_df) > 0)
    p <- p + geom_point(data = scatter_df, aes(x = t, y = pl, colour = factor(g)),
                        size = 2.2, alpha = 0.85)
  p + scale_colour_manual(values = cols, labels = l0_labels(d),
                          name = "[L]0 (nM)") +
    labs(title = "Dissociation", x = "Time (s)",
         y = "Concentration of Bound Protein ([PL], nM)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "right", legend.text = element_text(size = 11),
          plot.title = element_text(size = 13, face = "bold"))
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
  title = "Binding Curve Viewer — Kinetics of Association and Dissociation",
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),

  tags$head(tags$style(HTML("
    body { padding: 0 14px 20px 14px; }
    .p2-header {
      padding: 10px 2px 8px 2px;
      border-bottom: 2px solid #2c3e50;
      margin-bottom: 16px;
    }
    .p2-header h1 { display:inline; font-size:1.42rem; font-weight:650; margin:0; }
    .p2-header .home-link {
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
    .dist-bar {
      display:flex; align-items:center; flex-wrap:wrap; gap:6px 22px;
      margin:12px 0 4px 0;
    }
    .dist-bar .hist-label { margin:0 0 0 2px; }
    .dist-bar .form-group, .dist-bar .shiny-input-container { margin-bottom:0; }
    .dist-note { font-size:12.5px; color:#555; margin:2px 0 0 2px; }
    .param-card { margin-top:14px; }
    .group-table { font-size:12.5px; margin-bottom:4px; }
    .group-table th, .group-table td { padding:2px 8px !important; white-space:nowrap;
                                       text-align:right; }
    .group-table th:first-child, .group-table td:first-child { text-align:left; }
    .group-table .unmet { color:#a04000; font-weight:700; }
  "))),

  # p2_header
  div(class = "p2-header",
    h1("Kinetics of Association and Dissociation"),
    tags$a(class = "home-link", href = "/binding-curve-viewer/", "[Back to Home]")
  ),

  # Main row: kinetic parameters | association plot | dissociation plot
  fluidRow(
    column(2,
      numericInput("p0_val",
        HTML("Initial Conc. of Fixed Component (<b>[P]<sub>0</sub></b>, nM)
              <br><small style='color:#666;font-weight:400'>0.001 &ndash; 50 000 nM</small>"),
        value = DEFAULTS$p0, min = P0_MIN, max = P0_MAX, step = 1),
      numericInput("koff_val",
        HTML("Off-rate Constant (<b><i>k</i><sub>off</sub></b>, s<sup>-1</sup>)
              <br><small style='color:#666;font-weight:400'>1e-6 &ndash; 1 s<sup>-1</sup></small>"),
        value = DEFAULTS$koff, min = KOFF_MIN, max = KOFF_MAX, step = 0.000001),
      numericInput("kon_val",
        HTML("On-rate Constant (<b><i>k</i><sub>on</sub></b>, M<sup>-1</sup>s<sup>-1</sup>)
              <br><small style='color:#666;font-weight:400'>1e3 &ndash; 1e8 M<sup>-1</sup>s<sup>-1</sup></small>"),
        value = DEFAULTS$kon_M, min = KON_MIN, max = KON_MAX, step = 1000),
      tags$p(style = "font-size:12px; color:#666;", HTML(
        "The ligand concentrations [L]<sub>0</sub> are set by the simulation design
         (20%&ndash;90% of [P]<sub>0</sub> bound at equilibrium) and listed below.
         Plots and table update when a simulation is started or cleared."))
    ),
    column(5, plotOutput("plot_kinetics_on",    height = "420px")),
    column(5, plotOutput("plot_kinetics_disso", height = "420px"))
  ),

  # Navigation arrows — only the first 10 simulated experiments are browsable
  fluidRow(
    column(10, offset = 2,
      div(class = "nav-bar",
        actionButton("nav_left",  HTML("&#9664;"), class = "btn btn-primary btn-sm"),
        uiOutput("nav_text"),
        actionButton("nav_right", HTML("&#9654;"), class = "btn btn-primary btn-sm")
      )
    )
  ),

  # Binding Parameters card: Kd, dissociation times, per-[L]0 table
  div(class = "param-card",
    bslib::card(
      bslib::card_header(HTML("<b>Binding Parameters</b>")),
      bslib::card_body(uiOutput("kin_card_body"))
    )
  ),

  # Monte Carlo Simulation Console
  fluidRow(column(12,
    div(class = "sim-console",
      tags$h5("Monte Carlo Simulation",
              style = "margin:0 0 12px 0; font-weight:700; color:#2b579a;"),

      fluidRow(
        column(2,
          selectInput("expn", "Number of Simulated Experiments",
            choices = c("1000" = "1000", "5000" = "5000", "10,000" = "10000",
                        "100,000" = "100000"), selected = "1000")
        ),
        column(2,
          selectInput("conc_pts", "Number of Ligand Concentrations per Experiment",
            choices = as.character(3:10), selected = as.character(DEFAULTS$n_groups))
        ),
        column(2,
          selectInput("time_pts", "Number of Time Points per Experiment",
            choices = as.character(10:15), selected = as.character(DEFAULTS$time_pts))
        ),
        column(2,
          selectInput("repn", "Replicates per Time Point",
            choices = as.character(1:5), selected = "1")
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

      # Time-point design note | Start + Clear
      fluidRow(
        column(6,
          tags$p(style = "font-size:12.5px; color:#555; margin:8px 0 0 0;", HTML(
            "Time points: equal steps of [PL] for every [L]<sub>0</sub> (association
             0&ndash;99.9%, dissociation 100&ndash;0.1% of [PL]<sub>eq</sub>)."))
        ),
        column(4,
          tags$div(style = "margin-top:4px;",
            actionButton("start_sim", "Start Simulation",
              class = "btn btn-success btn-sm", style = "margin-right:6px;"),
            actionButton("clear_sim", "Clear Simulation",
              class = "btn btn-danger btn-sm")
          )
        )
      ),

      uiOutput("progress_ui")
    )
  )),

  # Fitting result panels
  fluidRow(
    column(4, div(class = "res-panel", uiOutput("res_disso"))),
    column(4, div(class = "res-panel", uiOutput("res_sec_order"))),
    column(4, div(class = "res-panel", uiOutput("res_pseudo")))
  ),

  # Distribution plots (histograms or normal Q-Q plots), one row per model.
  # height = "auto": each plot is drawn as tall as it is wide.
  div(class = "dist-bar",
    div(class = "hist-label", "Fitted Parameter Distributions"),
    radioButtons("dist_plot", NULL, inline = TRUE,
      choices  = c("Histogram" = "hist", "Normal Q-Q plot" = "qq"),
      selected = "hist")
  ),
  div(class = "hist-label", HTML(
    "Dissociation Exp. &mdash; <span style='text-transform:none'>fitted <i>k</i><sub>off</sub></span>")),
  fluidRow(
    column(3, plotOutput("hist_koff",  height = "auto")),
    column(3, plotOutput("hist_lkoff", height = "auto"))
  ),
  div(class = "hist-label", HTML(
    "Second Order Binding Equation &mdash; <span style='text-transform:none'>fitted
     <i>k</i><sub>on</sub> and <i>K</i><sub>d</sub> = <i>k</i><sub>off</sub>/<i>k</i><sub>on</sub></span>")),
  fluidRow(
    column(3, plotOutput("hist_kon_sec",  height = "auto")),
    column(3, plotOutput("hist_lkon_sec", height = "auto")),
    column(3, plotOutput("hist_kd_sec",   height = "auto")),
    column(3, plotOutput("hist_lkd_sec",  height = "auto"))
  ),
  div(class = "hist-label", HTML(
    "Pseudo-First-Order Binding Equation &mdash; <span style='text-transform:none'>fitted
     <i>k</i><sub>on</sub> and <i>K</i><sub>d</sub> = <i>k</i><sub>off</sub>/<i>k</i><sub>on</sub></span>")),
  fluidRow(
    column(3, plotOutput("hist_kon_pseudo",  height = "auto")),
    column(3, plotOutput("hist_lkon_pseudo", height = "auto")),
    column(3, plotOutput("hist_kd_pseudo",   height = "auto")),
    column(3, plotOutput("hist_lkd_pseudo",  height = "auto"))
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
    # Committed design (parameters, [L]0 series, time points, true curves).
    # Updated only by Start Simulation or Clear Simulation.
    design = build_kin_design(DEFAULTS$p0, DEFAULTS$koff, DEFAULTS$kon_M,
                              DEFAULTS$n_groups, DEFAULTS$time_pts),

    # Raw data for the first <= 10 experiments (navigation / scatter overlay)
    all_exps_asso  = NULL,   # data.frame(exp_idx, t, g, pl)
    all_exps_disso = NULL,   # data.frame(exp_idx, t, g, pl)

    # Estimates across all experiments (kon in M^-1 s^-1, Kd in nM) and outcomes
    koff_vec = NULL, kon_sec_vec = NULL, kon_pseudo_vec = NULL,
    kd_sec_vec = NULL, kd_pseudo_vec = NULL,
    st_koff = NULL, st_kon_sec = NULL, st_kon_pseudo = NULL,
    st_kd_sec = NULL, st_kd_pseudo = NULL,
    n_skipped = 0L,          # kon fits not run because the koff fit was unusable

    ks_results = NULL,       # bootstrap KS results (NULL: not run yet)

    current_idx = 1L,
    total_expn  = 0L,
    status_msg  = "Simulation not started."
  )

  clear_results <- function() {
    sim$all_exps_asso <- NULL; sim$all_exps_disso <- NULL
    sim$koff_vec <- NULL; sim$kon_sec_vec <- NULL; sim$kon_pseudo_vec <- NULL
    sim$kd_sec_vec <- NULL; sim$kd_pseudo_vec <- NULL
    sim$st_koff <- NULL; sim$st_kon_sec <- NULL; sim$st_kon_pseudo <- NULL
    sim$st_kd_sec <- NULL; sim$st_kd_pseudo <- NULL
    sim$n_skipped <- 0L; sim$ks_results <- NULL
    sim$total_expn <- 0L; sim$current_idx <- 1L
  }

  # Current experiment's scatter subsets (reactive on navigation only)
  current_sc_asso <- reactive({
    if (is.null(sim$all_exps_asso)) return(NULL)
    sim$all_exps_asso[sim$all_exps_asso$exp_idx == sim$current_idx, ]
  })
  current_sc_disso <- reactive({
    if (is.null(sim$all_exps_disso)) return(NULL)
    sim$all_exps_disso[sim$all_exps_disso$exp_idx == sim$current_idx, ]
  })

  # Plots — depend only on the committed design and navigation
  output$plot_kinetics_on <- renderPlot({
    make_kin_on_plot(sim$design, current_sc_asso())
  })
  output$plot_kinetics_disso <- renderPlot({
    make_kin_disso_plot(sim$design, current_sc_disso())
  })

  # Binding Parameters card
  output$kin_card_body <- renderUI({
    d <- sim$design
    unmet <- d$l0 / d$p0 < PSEUDO_RATIO
    cell <- function(x, flag = rep(FALSE, length(x)))
      paste0("<td", ifelse(flag, " class='unmet'", ""), ">", x, "</td>", collapse = "")
    row <- function(label, cells) paste0("<tr><td>", label, "</td>", cells, "</tr>")
    table_html <- paste0(
      "<div style='overflow-x:auto;'>",
      "<table class='table table-sm group-table'>",
      row("<b>[L]<sub>0</sub> (nM)</b>",
          cell(paste0("<b>", vapply(d$l0, fmt, character(1), digits = 4), "</b>"))),
      row("[L]<sub>0</sub>/[P]<sub>0</sub>",
          cell(vapply(d$l0 / d$p0, fmt, character(1), digits = 3), unmet)),
      row("Fraction bound at equilibrium",
          cell(sprintf("%.1f%%", 100 * d$pl1 / d$p0))),
      row("<i>t</i><sub>0.99</sub>, 2nd-order (s)",
          cell(vapply(d$t99, fmt, character(1), digits = 3))),
      row("<i>t</i><sub>pseudo-0.99</sub> / <i>t</i><sub>0.99</sub>",
          cell(sprintf("%.2f", d$pt99_ratio))),
      "</table></div>")
    note_html <- if (any(unmet)) {
      sprintf(paste0(
        "<span class='res-warn'>&#9888; The pseudo-first-order condition ",
        "[L]<sub>0</sub> &ge; %d&times;[P]<sub>0</sub> is not met for %d of %d ",
        "concentrations (orange); pseudo-first-order <i>k</i><sub>on</sub> and ",
        "<i>K</i><sub>d</sub> are expected to be biased.</span>"),
        PSEUDO_RATIO, sum(unmet), length(unmet))
    } else {
      sprintf(paste0("The pseudo-first-order condition [L]<sub>0</sub> &ge; ",
                     "%d&times;[P]<sub>0</sub> is met for all concentrations."), PSEUDO_RATIO)
    }
    design_html <- sprintf(
      "Time points: equal steps of [PL] per [L]<sub>0</sub> (%d points each).", d$time_pts)
    HTML(paste0(
      "<div class='row'>",
      "<div class='col-md-3' style='font-size:13px; line-height:1.6;'>",
        "<b><i>K</i><sub>d</sub></b> = <i>k</i><sub>off</sub>/<i>k</i><sub>on</sub> = ",
        fmt(d$kd), " nM<br>",
        "[P]<sub>0</sub> = ", fmt(d$p0), " nM<br>",
        "Half-life of dissociation (ln 2/<i>k</i><sub>off</sub>) = ", fmt(d$disso_t50, 3), " s<br>",
        "Time to 99% dissociation = ", fmt(d$disso_t99, 3), " s",
      "</div>",
      "<div class='col-md-9'>", table_html,
        "<div style='font-size:12.5px;'>", note_html, "<br>", design_html, "</div>",
      "</div></div>"))
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

  # Progress / status message
  output$progress_ui <- renderUI({
    tags$p(HTML(sim$status_msg),
           style = "font-weight:bold; color:#1a5276; margin:10px 0 0 0; font-size:13px;")
  })

  # Result panels -----------------------------------------------------------

  ks_line <- function(label, r) {
    if (is.null(r) || is.na(r$D))
      return(paste0(label, ": not available",
                    if (!is.null(r$msg)) paste0(" (", r$msg, ")") else ""))
    paste0(label, ": <i>D</i> = ", fmt(r$D, 3), ", ", fmt_p(r$p),
           " (<i>n</i> = ", r$n, ")")
  }

  # "Fit outcomes" line for one fitted parameter
  outcome_line <- function(what, status, n_skipped = 0L) {
    n_all <- length(status)
    cnt <- table(factor(status, levels = c("ok", "nonpos", "nonconv", "failed")))
    if (cnt[["ok"]] == n_all)
      return(sprintf("Fit outcomes (%s): all %d fits converged normally.<br>", what, n_all))
    skipped <- if (n_skipped > 0L)
      sprintf(" (incl. %d not fitted because the <i>k</i><sub>off</sub> fit was unusable)", n_skipped) else ""
    sprintf(paste0("Fit outcomes (%s, of %d): %d normal &middot; %d with estimate &le; 0 ",
                   "&middot; %d not converged &middot; %d errors%s<br>"),
            what, n_all, cnt[["ok"]], cnt[["nonpos"]], cnt[["nonconv"]], cnt[["failed"]], skipped)
  }

  # Summary line for one parameter: linear scale uses "ok" + "nonpos"
  stats_line <- function(label_html, unit_html, vec, status) {
    v <- vec[status %in% c("ok", "nonpos")]
    if (length(v) < 2L)
      return(paste0(label_html, ": too few usable estimates<br>"))
    q <- quantile(v, c(0.025, 0.5, 0.975), names = FALSE)
    paste0(label_html, " (", unit_html, "), <i>n</i> = ", length(v), ": ",
           "mean &plusmn; SD = ", fmt(mean(v)), " &plusmn; ", fmt(sd(v)),
           "; median = ", fmt(q[2]),
           " (2.5&ndash;97.5%: ", fmt(q[1]), "&ndash;", fmt(q[3]), ")<br>")
  }

  # Warning line about excluded / non-positive estimates
  warn_line <- function(status_list) {
    nonpos  <- sum(vapply(status_list, function(s) sum(s == "nonpos"), numeric(1)))
    exclude <- sum(vapply(status_list, function(s) sum(s %in% c("nonconv", "failed")), numeric(1)))
    w <- character(0)
    if (nonpos > 0)
      w <- c(w, "estimates &le; 0 are kept on the linear scale and excluded on log10")
    if (exclude > 0)
      w <- c(w, "fits that did not converge or ended in an error are excluded from all statistics")
    if (length(w)) paste0("<span class='res-warn'>&#9888; ", paste(w, collapse = "; "), ".</span><br>") else ""
  }

  ks_header <- paste0(
    "Normality test (Kolmogorov&ndash;Smirnov, mean and SD estimated from the data; ",
    "parametric bootstrap, <i>B</i> = ", KS_B, "):<br>")

  kon_lab <- "Fitted <i>k</i><sub>on</sub>"
  kon_unit <- "M<sup>-1</sup>s<sup>-1</sup>"
  kd_lab  <- "<i>K</i><sub>d</sub> = <i>k</i><sub>off</sub>/<i>k</i><sub>on</sub>"

  output$res_disso <- renderUI({
    if (is.null(sim$koff_vec) || is.null(sim$ks_results))
      return(HTML("<b>Dissociation Exp. Results:</b><br>Simulation not started."))
    ks <- sim$ks_results
    HTML(paste0(
      "<b>Dissociation Exp. Results:</b><br>",
      outcome_line("<i>k</i><sub>off</sub>", sim$st_koff),
      stats_line("Fitted <i>k</i><sub>off</sub>", "s<sup>-1</sup>", sim$koff_vec, sim$st_koff),
      ks_header,
      "&nbsp;&nbsp;", ks_line("<i>k</i><sub>off</sub>", ks$koff_lin), "<br>",
      "&nbsp;&nbsp;", ks_line("log10(<i>k</i><sub>off</sub>)", ks$koff_log), "<br>",
      warn_line(list(sim$st_koff))))
  })

  make_assoc_panel <- function(title, kon, st_kon, kd, st_kd, ks_kon_lin, ks_kon_log,
                               ks_kd_lin, ks_kd_log, extra = "") {
    HTML(paste0(
      "<b>", title, " Results:</b><br>",
      outcome_line("<i>k</i><sub>on</sub>", st_kon, sim$n_skipped),
      stats_line(kon_lab, kon_unit, kon, st_kon),
      stats_line(kd_lab, "nM", kd, st_kd),
      ks_header,
      "&nbsp;&nbsp;", ks_line("<i>k</i><sub>on</sub>", ks_kon_lin), "<br>",
      "&nbsp;&nbsp;", ks_line("log10(<i>k</i><sub>on</sub>)", ks_kon_log), "<br>",
      "&nbsp;&nbsp;", ks_line("<i>K</i><sub>d</sub>", ks_kd_lin), "<br>",
      "&nbsp;&nbsp;", ks_line("log10(<i>K</i><sub>d</sub>)", ks_kd_log), "<br>",
      warn_line(list(st_kon, st_kd)), extra))
  }

  output$res_sec_order <- renderUI({
    if (is.null(sim$kon_sec_vec) || is.null(sim$ks_results))
      return(HTML("<b>Second Order Binding Equation Results:</b><br>Simulation not started."))
    ks <- sim$ks_results
    make_assoc_panel("Second Order Binding Equation",
                     sim$kon_sec_vec, sim$st_kon_sec, sim$kd_sec_vec, sim$st_kd_sec,
                     ks$kon_sec_lin, ks$kon_sec_log, ks$kd_sec_lin, ks$kd_sec_log)
  })

  output$res_pseudo <- renderUI({
    if (is.null(sim$kon_pseudo_vec) || is.null(sim$ks_results))
      return(HTML("<b>Pseudo-First-Order Binding Equation Results:</b><br>Simulation not started."))
    ks <- sim$ks_results
    d <- sim$design
    n_unmet <- sum(d$l0 / d$p0 < PSEUDO_RATIO)
    extra <- if (n_unmet > 0)
      sprintf(paste0("<span class='res-warn'>&#9888; [L]<sub>0</sub> &lt; %d&times;[P]<sub>0</sub> ",
                     "for %d of %d concentrations: the pseudo-first-order assumption ",
                     "is not met and bias is expected.</span>"),
              PSEUDO_RATIO, n_unmet, length(d$l0)) else ""
    make_assoc_panel("Pseudo-First-Order Binding Equation",
                     sim$kon_pseudo_vec, sim$st_kon_pseudo, sim$kd_pseudo_vec, sim$st_kd_pseudo,
                     ks$kon_pseudo_lin, ks$kon_pseudo_log, ks$kd_pseudo_lin, ks$kd_pseudo_log,
                     extra)
  })

  # Distribution plots ------------------------------------------------------
  # Linear scale: "ok" + "nonpos" estimates; log10 scale: "ok" only.
  make_dist_output <- function(out_id, vec_name, st_name, title, unit_lab, log_scale) {
    output[[out_id]] <- renderPlot({
      vec <- sim[[vec_name]]; st <- sim[[st_name]]
      if (is.null(vec)) return(NULL)
      if (log_scale) {
        v <- log10(vec[st == "ok"])
        ttl <- sprintf("%s: log10  (n = %d)", title, length(v))
        lab <- paste0("log10(", unit_lab, ")")
      } else {
        v <- vec[st %in% c("ok", "nonpos")]
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
  make_dist_output("hist_koff",        "koff_vec",       "st_koff",       "koff",          "koff (s-1)",    FALSE)
  make_dist_output("hist_lkoff",       "koff_vec",       "st_koff",       "koff",          "koff (s-1)",    TRUE)
  make_dist_output("hist_kon_sec",     "kon_sec_vec",    "st_kon_sec",    "2nd-order kon", "kon (M-1 s-1)", FALSE)
  make_dist_output("hist_lkon_sec",    "kon_sec_vec",    "st_kon_sec",    "2nd-order kon", "kon (M-1 s-1)", TRUE)
  make_dist_output("hist_kd_sec",      "kd_sec_vec",     "st_kd_sec",     "2nd-order Kd",  "Kd (nM)",       FALSE)
  make_dist_output("hist_lkd_sec",     "kd_sec_vec",     "st_kd_sec",     "2nd-order Kd",  "Kd (nM)",       TRUE)
  make_dist_output("hist_kon_pseudo",  "kon_pseudo_vec", "st_kon_pseudo", "Pseudo kon",    "kon (M-1 s-1)", FALSE)
  make_dist_output("hist_lkon_pseudo", "kon_pseudo_vec", "st_kon_pseudo", "Pseudo kon",    "kon (M-1 s-1)", TRUE)
  make_dist_output("hist_kd_pseudo",   "kd_pseudo_vec",  "st_kd_pseudo",  "Pseudo Kd",     "Kd (nM)",       FALSE)
  make_dist_output("hist_lkd_pseudo",  "kd_pseudo_vec",  "st_kd_pseudo",  "Pseudo Kd",     "Kd (nM)",       TRUE)

  # Input validation; returns NULL when all inputs are valid --------------------
  check_inputs <- function(need_seed = TRUE) {
    is_num <- function(v) length(v) == 1L && is.numeric(v) && is.finite(v)
    in_rng <- function(v, lo, hi) is_num(v) && v >= lo && v <= hi
    msgs <- character(0)
    if (!in_rng(input$p0_val, P0_MIN, P0_MAX))
      msgs <- c(msgs, "[P]<sub>0</sub> must be between 0.001 and 50 000 nM.")
    if (!in_rng(input$koff_val, KOFF_MIN, KOFF_MAX))
      msgs <- c(msgs, "<i>k</i><sub>off</sub> must be between 1e-6 and 1 s<sup>-1</sup>.")
    if (!in_rng(input$kon_val, KON_MIN, KON_MAX))
      msgs <- c(msgs, "<i>k</i><sub>on</sub> must be between 1e3 and 1e8 M<sup>-1</sup>s<sup>-1</sup>.")
    if (need_seed) {
      s <- input$seed_val
      if (!is_num(s) || s < 1 || s > .Machine$integer.max || s != round(s))
        msgs <- c(msgs, "Random seed must be a positive integer.")
    }
    if (length(msgs)) paste(msgs, collapse = "<br>") else NULL
  }

  # Build the design from the current inputs (after validation)
  design_from_inputs <- function() {
    build_kin_design(as.numeric(input$p0_val), as.numeric(input$koff_val),
                     as.numeric(input$kon_val), as.integer(input$conc_pts),
                     as.integer(input$time_pts))
  }

  # Start Simulation ----------------------------------------------------------
  observeEvent(input$start_sim, {

    expN <- as.integer(input$expn)

    # Server-mode guard: ExpN > 1000 triggers a "run locally" warning
    if (!IS_LOCAL && expN > 1000L) {
      showModal(modalDialog(
        title = tags$span(style = "color:#c0392b; font-weight:700;",
                          "\u26a0\ufe0f  Large Simulation \u2014 Server Load Warning"),
        HTML(paste0(
          "Running <b>Number of Simulated Experiments &gt; 1000</b> on the web server ",
          "may cause heavy load and slow the application for all users.<br><br>",
          "Please <b>download and run the app locally</b>:<br>",
          "<a href='https://github.com/ydu-sci/Binding_Curve_Viewer' target='_blank' style='font-weight:600;'>",
          "github.com/ydu-sci/Binding_Curve_Viewer</a><br>(app file: <code>src/shiny-apps/kin-sim/app.R</code>)<br><br>",
          "Once downloaded, launch it from that folder in <b>local mode</b>:<br>",
          "<code style='background:#f1f1f1; padding:3px 6px; border-radius:3px;'>",
          "BINDING_APP_LOCAL=true Rscript -e \"shiny::runApp('app.R')\"</code><br>",
          "<small style='color:#555;'>Or from R console: ",
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

    repN  <- as.integer(input$repn)
    ratio <- as.numeric(input$ratio_sd)
    seed  <- as.integer(input$seed_val)

    # Commit the design and reset previous results
    d <- design_from_inputs()
    sim$design <- d
    clear_results()
    sim$status_msg <- sprintf("Running (ExpN=%d, seed=%d) ...", expN, seed)

    tryCatch(
      withProgress(
        message = sprintf("Monte Carlo Simulation (ExpN = %d)", expN),
        value = 0,
        {
          set.seed(seed)

          # Full data vectors: every time point repeated repN times, groups stacked
          t_asso_full  <- rep(d$asso$t,  each = repN)
          pl_asso_theo <- rep(d$asso$pl, each = repN)
          g_asso_full  <- rep(d$asso$g,  each = repN)
          l0_asso_full <- d$l0[g_asso_full]
          t_disso_full  <- rep(d$disso$t,  each = repN)
          pl_disso_theo <- rep(d$disso$pl, each = repN)
          g_disso_full  <- rep(d$disso$g,  each = repN)
          n_asso  <- length(t_asso_full)
          n_disso <- length(t_disso_full)

          # Noise (drawn in the same order as v2) and starting values
          noise_sd <- ratio * d$p0
          noise_asso_mat  <- matrix(rnorm(n_asso  * expN, sd = noise_sd), nrow = n_asso,  ncol = expN)
          noise_disso_mat <- matrix(rnorm(n_disso * expN, sd = noise_sd), nrow = n_disso, ncol = expN)
          pl0_init_vec  <- d$pl1 * 0.5
          p0_init       <- d$p0 * 0.5
          pleq_init_vec <- d$pl1 * 0.5
          koff_init     <- d$koff * 0.5
          kon_init      <- d$kon * 0.5

          koff_est <- conv_off <- numeric(expN)
          kon2_est <- conv_2   <- rep(NA_real_, expN)
          konp_est <- conv_p   <- rep(NA_real_, expN)
          store_a <- vector("list", min(10L, expN)); store_d <- vector("list", min(10L, expN))
          upd_every <- max(1L, expN %/% 100L)

          # Main loop: (1) dissociation fit -> koff; (2) both association fits on
          # the SAME noisy association data with koff fixed, only when the koff
          # fit is usable (converged, koff > 0).
          for (i in seq_len(expN)) {
            y_a <- pl_asso_theo  + noise_asso_mat[, i]
            y_d <- pl_disso_theo + noise_disso_mat[, i]
            if (i <= 10L) {
              store_a[[i]] <- data.frame(exp_idx = i, t = t_asso_full,  g = g_asso_full,  pl = y_a)
              store_d[[i]] <- data.frame(exp_idx = i, t = t_disso_full, g = g_disso_full, pl = y_d)
            }

            r_off <- fit_disso_exp(t_disso_full, y_d, g_disso_full, pl0_init_vec, koff_init)
            koff_est[i] <- r_off[["est"]]; conv_off[i] <- r_off[["conv"]]

            if (identical(est_status(r_off[["est"]], r_off[["conv"]]), "ok")) {
              r2 <- fit_sec_order_exp(t_asso_full, l0_asso_full, g_asso_full, y_a,
                                      p0_init, r_off[["est"]], kon_init)
              rp <- fit_pseudo_exp(t_asso_full, l0_asso_full, g_asso_full, y_a,
                                   pleq_init_vec, r_off[["est"]], kon_init)
              kon2_est[i] <- r2[["est"]]; conv_2[i] <- r2[["conv"]]
              konp_est[i] <- rp[["est"]]; conv_p[i] <- rp[["conv"]]
            }

            if (i %% upd_every == 0L || i == expN)
              incProgress(0.80 * upd_every / expN,
                          detail = sprintf("Fitting experiment %d / %d  (%.0f%%)",
                                           i, expN, 100 * i / expN))
          }

          st_koff <- est_status(koff_est, conv_off)
          st_kon2 <- est_status(kon2_est, conv_2)
          st_konp <- est_status(konp_est, conv_p)
          n_skipped <- sum(st_koff != "ok")

          # Kinetic Kd = koff / kon in nM (kon still in nM^-1 s^-1 here)
          kd2 <- koff_est / kon2_est
          kdp <- koff_est / konp_est
          st_kd2 <- kd_status(st_koff, st_kon2, kd2)
          st_kdp <- kd_status(st_koff, st_konp, kdp)

          # kon in M^-1 s^-1 for display
          kon2_M <- kon2_est * 1e9
          konp_M <- konp_est * 1e9

          # Bootstrap KS tests (progress 80% -> 100%)
          null_cache <- new.env(parent = emptyenv())
          lin <- function(v, s) v[s %in% c("ok", "nonpos")]
          lg  <- function(v, s) log10(v[s == "ok"])
          ks_specs <- list(
            koff_lin       = list(lin(koff_est, st_koff), "koff"),
            koff_log       = list(lg(koff_est,  st_koff), "log10(koff)"),
            kon_sec_lin    = list(lin(kon2_M, st_kon2),   "2nd-order kon"),
            kon_sec_log    = list(lg(kon2_M,  st_kon2),   "2nd-order log10(kon)"),
            kon_pseudo_lin = list(lin(konp_M, st_konp),   "pseudo kon"),
            kon_pseudo_log = list(lg(konp_M,  st_konp),   "pseudo log10(kon)"),
            kd_sec_lin     = list(lin(kd2, st_kd2),       "2nd-order Kd"),
            kd_sec_log     = list(lg(kd2,  st_kd2),       "2nd-order log10(Kd)"),
            kd_pseudo_lin  = list(lin(kdp, st_kdp),       "pseudo Kd"),
            kd_pseudo_log  = list(lg(kdp,  st_kdp),       "pseudo log10(Kd)"))
          ks_res <- list()
          for (k in names(ks_specs)) {
            incProgress(0.20 / length(ks_specs), detail = paste0("KS test: ", ks_specs[[k]][[2]], " ..."))
            x <- ks_specs[[k]][[1]]
            ks_res[[k]] <- tryCatch(
              ks_boot_test(x, null_cache),
              error = function(e) list(n = length(x), D = NA_real_, p = NA_real_,
                                       msg = conditionMessage(e)))
          }

          # Commit results
          sim$all_exps_asso  <- do.call(rbind, store_a)
          sim$all_exps_disso <- do.call(rbind, store_d)
          sim$koff_vec <- koff_est; sim$st_koff <- st_koff
          sim$kon_sec_vec <- kon2_M; sim$st_kon_sec <- st_kon2
          sim$kon_pseudo_vec <- konp_M; sim$st_kon_pseudo <- st_konp
          sim$kd_sec_vec <- kd2; sim$st_kd_sec <- st_kd2
          sim$kd_pseudo_vec <- kdp; sim$st_kd_pseudo <- st_kdp
          sim$n_skipped <- n_skipped
          sim$ks_results <- ks_res
          sim$total_expn  <- expN
          sim$current_idx <- 1L
          sim$status_msg  <- sprintf("Completed: %d runs, seed = %d.", expN, seed)
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
  # Clears all results and recomputes the design (plots, table) from the
  # current inputs, so the design can be previewed before a simulation.
  observeEvent(input$clear_sim, {
    err <- check_inputs(need_seed = FALSE)
    clear_results()
    if (!is.null(err)) {
      showNotification(HTML(err), type = "error", duration = 8)
      sim$status_msg <- "Simulation cleared (design not updated: invalid inputs)."
      return()
    }
    sim$design <- design_from_inputs()
    sim$status_msg <- "Simulation cleared."
  })
}

# Section 6 — Launch ------------------------------------------------------

shinyApp(ui = ui, server = server)
