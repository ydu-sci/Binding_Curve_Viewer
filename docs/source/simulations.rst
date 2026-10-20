Simulations in R
================

Four R Shiny applications complement the visualization tools of Binding Curve
Viewer. Three of them simulate binding experiments by Monte Carlo and show how
readout noise, the concentration design and the choice of fitting model affect
the estimated binding constants: for each simulated experiment, noise is added
to the exact [PL] values and the data are fitted with an exact mechanistic model
and with a common approximation. The fourth application shows, without noise,
how the incubation time of a saturation experiment affects the apparent
*K*\ :sub:`d`. All four can be used on the
`website <https://www.eplatton.net/binding-curve-viewer/>`_ or run locally.

Requirements and Run Modes
--------------------------

The following R packages are required (see :doc:`overview` for the minimum
versions)::

   install.packages(c("shiny", "bslib", "ggplot2", "minpack.lm"))

The Monte Carlo applications can be launched in two modes. In **server mode**
(the default), a run is limited to 1000 simulated experiments; selecting a larger
number shows a message explaining how to run the application locally. In
**local mode**, 5000, 10,000 or 100,000 simulated experiments per run are
available:

.. code-block:: bash

   # Server mode
   Rscript -e "shiny::runApp('app.R')"

   # Local mode, from the shell
   BINDING_APP_LOCAL=true Rscript -e "shiny::runApp('app.R')"

.. code-block:: r

   # Local mode, from the R console
   Sys.setenv(BINDING_APP_LOCAL = "true")
   shiny::runApp("app.R")

Run time grows in proportion to the number of simulated experiments. As a rough
guide, the *K*\ :sub:`d` application fits 10,000 experiments in about half a
minute on one CPU core, and the kinetics and competition applications take
about 2.5 times as long. The kinetics application draws all noise before
fitting: 100,000 experiments with the largest design (10 ligand concentrations
× 15 time points × 5 replicates) need about 1.2 GB of memory, and the default
design about 130 MB.

The Time Dependency of Dissociation Constant application contains no random
component and is launched in the same way, without a local mode.

Common Features of the Monte Carlo Applications
-----------------------------------------------

- **Committed inputs.** Parameters are read when *Start Simulation* is pressed,
  so plots and results always refer to the same parameter set. All inputs are
  checked on the server before a run starts.
- **Readout noise.** Noise-free [PL] values are computed with exact equations,
  and additive Gaussian noise with an SD of 1.25%, 2.5% (default) or 5% of the
  protein concentration is added to every replicate. Errors in the
  concentrations themselves are not simulated.
- **Reproducibility.** The random seed is set and all noise is drawn in a fixed
  order before fitting, so the same inputs and seed reproduce a run exactly.
- **Individual experiments.** The first ten simulated experiments can be browsed
  with the ◀ and ▶ buttons, overlaid on the true curves.
- **Fitting.** Every simulated experiment is fitted by non-linear least squares,
  starting at half the true parameter values: with the Levenberg–Marquardt
  algorithm (``minpack.lm::nlsLM``) in the *K*\ :sub:`d` and competition
  applications, and with the Gauss–Newton algorithm (``nls``) in the kinetics
  application.
- **Fit outcomes.** Every fit is classified as normal, at the lower bound
  (estimate of 0 or below), not converged, or error, and the number of fits in
  each class is reported. Fits that did not converge or failed are excluded
  from all statistics; fits at the lower bound are kept on the linear scale and
  excluded on the log10 scale.
- **Summary statistics.** Mean, SD, median and 2.5–97.5% range of each
  estimated parameter on the linear and log10 scales.
- **Distributions.** Histograms with a fitted normal density (default) or
  normal Q–Q plots with a pointwise 95% band derived from the beta distribution
  of order statistics.
- **Normality tests.** Normality on the linear and log10 scales is tested with
  the Kolmogorov–Smirnov statistic *D* against a normal distribution whose mean
  and SD are estimated from the data. Because estimating these parameters
  invalidates the standard KS p-value, *p* is obtained by parametric bootstrap
  with *B* = 200: bootstrap samples of the same size are drawn from a standard
  normal distribution and their *D* is computed in the same way. This null
  distribution depends only on the sample size, so it is generated once per size
  and *p* = (1 + #{*D*\ :sup:`*` ≥ *D*})/(*B* + 1); the smallest reportable
  value is *p* < 0.005. *D*, *n* and *p* are reported; samples with fewer than
  five values are not tested.
- **Error handling.** Errors in a fit, a test or the whole simulation are
  reported in the application instead of ending the session.


Determination of Dissociation Constant (*K*\ :sub:`d`)
------------------------------------------------------

**File:** ``shiny-apps/kd-sim/app.R``

This application simulates saturation binding experiments. Noisy [PL] values are
generated at a series of ligand concentrations, and each simulated experiment is
fitted with two models:

- **Formula 1 (Exact):** the quadratic binding equation, which accounts for
  ligand depletion and is valid at any [P]\ :sub:`total`-to-*K*\ :sub:`d` ratio.
  The smaller root is computed in a cancellation-free form,
  2[P]\ :sub:`total`\ [L]\ :sub:`total` / (*S* + √(*S*\ :sup:`2` −
  4[P]\ :sub:`total`\ [L]\ :sub:`total`)), with
  *S* = [P]\ :sub:`total` + [L]\ :sub:`total` + *K*\ :sub:`d`.
- **Formula 2 (Approximate):** the hyperbolic model, which assumes
  [L]\ :sub:`free` ≈ [L]\ :sub:`total` and is valid only when
  [P]\ :sub:`total` ≪ *K*\ :sub:`d`. Otherwise it estimates approximately the
  apparent *K*\ :sub:`d`, the total ligand concentration at which half of the
  protein is bound.

Both models estimate *K*\ :sub:`d` and [P]\ :sub:`total`; the results focus on
*K*\ :sub:`d`. The ligand series is log-spaced from 0.001 × *K*\ :sub:`d` to the
concentration at which 99% of the protein is bound (default), or defined as a
serial dilution from a top concentration with a dilution factor of 2, 3, 3.16
(half-log) or 10. The theoretical fraction bound is shown for every
concentration. A [P]\ :sub:`total`/*K*\ :sub:`d` ratio above 10 gives a
titration-regime warning, and a ratio above 1000 stops the run.

**Key parameters:** *K*\ :sub:`d` and [P]\ :sub:`total` (default 100 nM each),
number of simulated experiments, number of ligand concentrations per experiment
(5–10, default 9), replicates per concentration (3–10, default 5), [PL] noise SD
as a fraction of [P]\ :sub:`total`, random seed, and the ligand concentration
series.


Kinetics of Association and Dissociation
----------------------------------------

**File:** ``shiny-apps/kin-sim/app.R``

This application simulates association and dissociation time courses to estimate
the dissociation rate constant *k*\ :sub:`off` and the association rate constant
*k*\ :sub:`on`. The ligand concentrations [L]\ :sub:`0` give 20–90% occupancy of
[P]\ :sub:`0` at equilibrium. For each [L]\ :sub:`0`, the time points are placed
at equal steps of [PL], from 0 to 99.9% of [PL]\ :sub:`eq` for association and
from 100% to 0.1% for dissociation. Noise-free [PL] values are computed with the
exact second-order equations, in which 1 − e\ :sup:`−x` is evaluated as
−expm1(−\ *x*) to keep full precision at short times.

Each simulated experiment is fitted in three steps:

- **Dissociation:** *k*\ :sub:`off` is estimated by fitting an exponential
  decay to the dissociation data of all ligand concentrations together.
- **Second-order binding equation:** *k*\ :sub:`on` is estimated from the
  association data with the exact second-order kinetic model, with
  *k*\ :sub:`off` fixed to the value from the dissociation fit.
- **Pseudo-first-order binding equation:** *k*\ :sub:`on` is estimated from the
  same association data with the pseudo-first-order approximation,
  [PL](*t*) = [PL]\ :sub:`eq`\ (1 − exp(−(*k*\ :sub:`on`\ [L]\ :sub:`0` +
  *k*\ :sub:`off`)\ *t*)), again with *k*\ :sub:`off` fixed.

*k*\ :sub:`on` is fitted only when the *k*\ :sub:`off` fit of the same experiment
is usable. The kinetic *K*\ :sub:`d` = *k*\ :sub:`off`/*k*\ :sub:`on` is reported
for both association models. For every ligand concentration, the application
lists [L]\ :sub:`0`/[P]\ :sub:`0`, the time to reach 99% of equilibrium and the
ratio of the pseudo-first-order to the exact time, and it flags concentrations
that do not satisfy the pseudo-first-order condition
[L]\ :sub:`0` ≥ 10[P]\ :sub:`0`. Both models are kept because they correspond to
different experiments: in surface plasmon resonance the analyte is continuously
replenished and the pseudo-first-order equation applies, whereas in radioligand
and other solution assays the binding partners are often present at comparable
concentrations and the second-order equation applies.

**Key parameters:** [P]\ :sub:`0` (default 10 nM), *k*\ :sub:`off` (default
0.01 s\ :sup:`−1`), *k*\ :sub:`on` (default 10\ :sup:`5` M\ :sup:`−1` s\ :sup:`−1`),
number of simulated experiments, number of ligand concentrations per experiment
(3–10, default 8), number of time points per experiment (10–15, default 10),
replicates per time point (1–5, default 1), [PL] noise SD as a fraction of
[P]\ :sub:`0`, and random seed.


Competitive Binding — Inhibition Curve
--------------------------------------

**File:** ``shiny-apps/comp-sim/app.R``

This application simulates competitive binding experiments in which an inhibitor
displaces a ligand from the protein. Equilibrium [PL] values at a series of
inhibitor concentrations are computed with Wang's exact cubic equation
(FEBS Lett. 1995, 360, 111–114): the trigonometric solution for the free protein
concentration is refined by bounded Newton steps on the mass balance, which keeps
full precision when [I]\ :sub:`0` or [L]\ :sub:`0` greatly exceeds
[P]\ :sub:`0`. The exact IC50 (free inhibitor) and apparent IC50 (total
inhibitor) are calculated in closed form (see :doc:`appendix3`), and the
theoretical inhibition is shown for every inhibitor concentration.

The inhibitor series is log-spaced from apparent IC50/1000 to apparent
IC50 × 1000 (default), or defined as a serial dilution from a top concentration
with a dilution factor of 2, 3, 3.16 (half-log) or 10. Each simulated experiment
is analysed in three ways:

- **Logistic model (4-parameter):** the apparent IC50 is obtained from a
  four-parameter logistic fit.
- **Conversion of the apparent IC50:** each fitted apparent IC50 is converted to
  *K*\ :sub:`i` with the Cheng–Prusoff equation (Biochem. Pharmacol. 1973, 22,
  3099–3108) and with Wang's group equation (Anal. Biochem. 2004, 332, 261–273).
- **Mechanistic model (Wang equation):** *K*\ :sub:`i` is fitted directly to the
  [PL]-versus-[I]\ :sub:`total` data with the exact competitive binding
  equation, with [L]\ :sub:`0` and *K*\ :sub:`d` fixed.

A [P]\ :sub:`0`/*K*\ :sub:`i` ratio above 10 gives a tight-binding warning, and a
ratio above 1000 stops the run.

**Key parameters:** [P]\ :sub:`0` (default 100 nM), [L]\ :sub:`0` (default
75 nM), *K*\ :sub:`d` and *K*\ :sub:`i` (default 100 nM each), number of
simulated experiments, number of inhibitor concentrations per experiment (5–10,
default 9), replicates per concentration (3–10, default 5), [PL] noise SD as a
fraction of [P]\ :sub:`0`, random seed, and the inhibitor concentration series.


Time Dependency of Dissociation Constant
----------------------------------------

**File:** ``shiny-apps/time-sim/app.R``

This application shows how the *K*\ :sub:`d` obtained from a saturation binding
experiment depends on the incubation time before measurement. It was inspired by
the study of Bedini et al. (J. Pineal Res. 2024, 76, e12941), who measured the
binding of the slowly dissociating radioligand 2-[\ :sup:`125`\ I]iodomelatonin
to melatonin receptors after 2 h and 20 h of incubation. No noise is added; all
[PL] values are theoretical.

Three curves are shown: the equilibrium saturation curve and the curves reached
after incubation times *t*\ :sub:`1` and *t*\ :sub:`2`, computed with the exact
second-order kinetic solution. For each incubation time, a series of 20
log-spaced ligand concentrations is chosen so that [PL] at that time spans 1–99%
of [P]\ :sub:`total`. Each series is fitted with the exact equilibrium equation
(Formula 1) by Levenberg–Marquardt least squares, with [P]\ :sub:`total` and
*K*\ :sub:`d` free and constrained to be non-negative. The apparent
*K*\ :sub:`d` is reported as a value and as a multiple of the true
*K*\ :sub:`d` = *k*\ :sub:`off`/*k*\ :sub:`on`, together with the outcome of the
fit. When the incubation time is too short, the apparent *K*\ :sub:`d` is
overestimated.

The application also reports the incubation time after which every ligand
concentration has reached 95% and 99% of its equilibrium [PL]. These times are
obtained from the closed-form inverse of the exact kinetic solution and
maximised over all ligand concentrations. A [P]\ :sub:`total`/*K*\ :sub:`d` ratio
above 10 gives a titration-regime warning, and a ratio above 1000 stops the run.

**Key parameters:** [P]\ :sub:`total` (0.001–50,000 nM, default 0.008 nM),
*k*\ :sub:`off` (10\ :sup:`−6`–1 s\ :sup:`−1`, default
1 × 10\ :sup:`−5` s\ :sup:`−1`), *k*\ :sub:`on` (10\ :sup:`3`–10\ :sup:`8`
M\ :sup:`−1` s\ :sup:`−1`, default 6.1 × 10\ :sup:`5` M\ :sup:`−1` s\ :sup:`−1`),
and the incubation times *t*\ :sub:`1` and *t*\ :sub:`2` (0.001–100,000 h,
default 2 h and 20 h).
