Release Notes
=============

Version 2.0 (October 2026)
--------------------------

Binding Curve Viewer v2.0 introduces exact IC50 calculations and four
R Shiny applications: three for Monte Carlo simulation of binding
experiments and one for the effect of incubation time on the apparent
:math:`K_\mathrm{d}`.

**Exact Calculation of IC50 and Apparent IC50**

Previous versions of Binding Curve Viewer approximated the IC50 and
apparent IC50. Version 2.0 replaces this approximation with exact 
analytical calculations based on the equilibrium mass-balance 
equations of the competitive binding system (see :doc:`appendix3`).

The calculation is implemented in two functions:

- ``inhibition_calc(p0, l0, kd, ki)`` — a Python function used in the
  Binding Curve Viewer web application.
- ``inhibitionCalc(p0, l0, kd, ki)`` — the equivalent JavaScript
  function for client-side computation.

Both functions derive :math:`[\mathrm{PL}]_{\mathrm{eq\text{-}0}}` from
the quadratic binding equation, determine the concentrations of all
species at 50% inhibition from the equilibrium mass-balance equations,
and return the exact IC50 (free inhibitor concentration) and the exact
apparent IC50 (total inhibitor concentration) without approximation.

**Four New R Shiny Applications**

Version 2.0 adds four standalone R Shiny applications. Three of them
simulate binding experiments by Monte Carlo: they add Gaussian noise to
the exact [PL] values, fit every simulated experiment with an exact
mechanistic model and with a common approximation, and show the
distributions of the fitted parameters. The fourth computes, without
noise, how a saturation experiment read before equilibrium overestimates
:math:`K_\mathrm{d}`. For details, see :doc:`simulations`.

- **Determination of Dissociation Constant** (``src/shiny-apps/kd-sim/app.R``):
  Simulates saturation binding experiments and fits :math:`K_\mathrm{d}`
  with the exact quadratic binding equation (Formula 1) and with the
  hyperbolic approximation (Formula 2), which neglects ligand depletion.

- **Kinetics of Association and Dissociation** (``src/shiny-apps/kin-sim/app.R``):
  Simulates association and dissociation time courses at several ligand
  concentrations, estimates :math:`k_\mathrm{off}` from the dissociation
  data and :math:`k_\mathrm{on}` from the association data with the exact
  second-order equation and with the pseudo-first-order approximation,
  and reports the kinetic
  :math:`K_\mathrm{d} = k_\mathrm{off}/k_\mathrm{on}` for both.

- **Competitive Binding — Inhibition Curve** (``src/shiny-apps/comp-sim/app.R``):
  Simulates inhibition curves, fits the apparent IC50 with a
  four-parameter logistic model, converts it to :math:`K_\mathrm{i}`
  with the Cheng–Prusoff equation and with Wang's group equation, and
  also fits :math:`K_\mathrm{i}` directly with Wang's exact competitive
  binding equation.

- **Time Dependency of Dissociation Constant** (``src/shiny-apps/time-sim/app.R``):
  Computes the saturation curves reached after two incubation times with
  the exact second-order kinetic solution, fits them with the equilibrium
  equation to give the apparent :math:`K_\mathrm{d}`, and reports the
  incubation time needed for every ligand concentration to reach 95% and
  99% of equilibrium.

The three Monte Carlo applications share the following features:

- Seeded random numbers drawn in a fixed order, so every run can be
  reproduced exactly.
- 1000 simulated experiments per run on the web server; 5000, 10,000 or
  100,000 in local mode.
- Classification of every fit (normal, at the lower bound, not
  converged, error); fits that did not converge or failed are excluded
  from the statistics.
- Summary statistics (mean, SD, median and 2.5–97.5% range) on the linear
  and log10 scales, shown as histograms or as normal Q–Q plots with
  pointwise 95% bands.
- Normality tests by a parametric bootstrap Kolmogorov–Smirnov test
  (:math:`B` = 200), which replaces the standard KS p-value; the
  KSgeneral package is no longer required.
- Warnings for designs outside the validity of the approximations, such as
  the titration regime or ligand concentrations below
  :math:`10[\mathrm{P}]_0` in the kinetic simulation.

Simulation of errors in the concentrations themselves (pipetting errors in
:math:`[\mathrm{P}]_\mathrm{total}`, :math:`[\mathrm{L}]_\mathrm{total}`
and :math:`[\mathrm{I}]_\mathrm{total}`), available in earlier versions
of the applications, has been removed; only random error in the [PL]
readout is simulated.

**Numerical Accuracy**

All binding equations in the R applications are evaluated in forms that
avoid subtracting nearly equal numbers: the smaller root of the quadratic
binding equation is computed as
:math:`2[\mathrm{P}][\mathrm{L}]/(S + \sqrt{S^2 - 4[\mathrm{P}][\mathrm{L}]})`,
the exact second-order kinetic solution uses ``expm1``, and the free
protein concentration from Wang's cubic equation is refined by bounded
Newton steps on the mass balance. Against 60-digit reference values the
maximum relative error is below :math:`10^{-13}`.
