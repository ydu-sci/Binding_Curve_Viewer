# Binding Curve Viewer

**Binding Curve Viewer** is a set of tools for visualizing and simulating the
equilibrium and kinetics of binding reactions with 1:1 stoichiometry and no
cooperativity, for planning and validating protein–ligand binding experiments.

- Online tools: <https://www.eplatton.net/binding-curve-viewer/>
- Documentation: <https://www.eplatton.net/binding-curve-viewer/doc/>

## Purpose

Binding constants (K<sub>d</sub>, k<sub>on</sub>, k<sub>off</sub>, IC<sub>50</sub>
and K<sub>i</sub>) are estimated by fitting a model to assay data, so their
accuracy depends on the assay design, the readout noise, the fitting model and
the incubation time. Binding Curve Viewer shows the exact theoretical curves and,
in version 2.0, simulates complete experiments, so that a design and an analysis
method can be tested before the experiment is run.

## Contents

**Visualization tools (Python/Bokeh).** `src/binding_curve_viewer_v2.py`
generates six interactive web pages: determination of the dissociation constant,
kinetics of association and dissociation, iso-affinity plot, competitive binding
kinetics (association and dissociation), and simulation of high-throughput
screening. Version 2.0 calculates the exact IC<sub>50</sub> and apparent
IC<sub>50</sub> in closed form.

**Interactive simulation tools (R Shiny).** Four standalone applications in
`src/shiny-apps/`:

| Application | Description |
| --- | --- |
| `kd-sim` | Monte Carlo simulation of saturation binding; K<sub>d</sub> fitted with the exact quadratic equation and with the hyperbola |
| `kin-sim` | Monte Carlo simulation of association and dissociation; k<sub>on</sub> fitted with the exact second-order and the pseudo-first-order equations |
| `comp-sim` | Monte Carlo simulation of competitive binding; apparent IC<sub>50</sub>, K<sub>i</sub> by the Cheng–Prusoff and Wang's group equations, and K<sub>i</sub> fitted directly with Wang's exact equation |
| `time-sim` | Apparent K<sub>d</sub> of a saturation experiment read before equilibrium, and the incubation time needed |

The three Monte Carlo applications report the distributions of the fitted
parameters (histograms or normal Q–Q plots), fit-outcome counts and parametric
bootstrap Kolmogorov–Smirnov tests. Runs are reproducible from the random seed.

## Installation

Requirements:

- Python ≥ 3.7.6 with Bokeh ≥ 3.4.0 and NumPy ≥ 1.26.0
- R ≥ 4.3 with shiny ≥ 1.13.0, bslib ≥ 0.4.0, ggplot2 ≥ 3.5.0 and minpack.lm

```bash
git clone https://github.com/ydu-sci/Binding_Curve_Viewer.git
cd Binding_Curve_Viewer
pip install -r requirements.txt
```

```r
install.packages(c("shiny", "bslib", "ggplot2", "minpack.lm"))
```

## Usage

Visualization tools (the six HTML pages are written to the current directory):

```bash
python src/binding_curve_viewer_v2.py
```

Simulation tools, from the R console in the repository folder:

```r
shiny::runApp("src/shiny-apps/kd-sim")
```

On the web server, a Monte Carlo run is limited to 1000 simulated experiments.
In local mode, 5000, 10,000 or 100,000 are available:

```r
Sys.setenv(BINDING_APP_LOCAL = "true")
shiny::runApp("src/shiny-apps/kd-sim")
```

Documentation:

```bash
cd docs
pip install -r requirements.txt
make html
```

## Repository structure

```text
src/binding_curve_viewer_v2.py   visualization tools (Python/Bokeh)
src/shiny-apps/                  interactive simulation tools (R Shiny)
src/scripts/                     scripts for the figures of the 2024 article
docs/                            Sphinx documentation
```

## Version 1.0

The version described in the 2024 article in *J. Chem. Inf. Model.* is available
under the tag [v1.0.0](https://github.com/ydu-sci/Binding_Curve_Viewer/tree/v1.0.0).

## How to cite

Du, Y. Binding Curve Viewer: Visualizing the Equilibrium and Kinetics of
Protein–Ligand Binding and Competitive Binding. *J. Chem. Inf. Model.* 2024,
64 (10), 4180–4192. <https://doi.org/10.1021/acs.jcim.4c00130>

## License

The code is released under the MIT License (see [LICENSE.txt](LICENSE.txt)).
The web pages are made available under the CC BY 4.0 license.
