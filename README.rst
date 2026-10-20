Binding Curve Viewer
====================

**Binding Curve Viewer** is a set of tools for visualizing and simulating the
equilibrium and kinetics of binding reactions with 1:1 stoichiometry and no
cooperativity, for experiment planning and validation.

- Online tools: https://www.eplatton.net/binding-curve-viewer/
- Documentation: https://www.eplatton.net/binding-curve-viewer/doc/

Version 2.0 contains two parts.

**Visualization tools (Python/Bokeh).** ``binding_curve_viewer_v2.py`` generates
six interactive web pages: determination of the dissociation constant, kinetics
of association and dissociation, iso-affinity plot, competitive binding kinetics
(association and dissociation), and simulation of high-throughput screening.
Version 2.0 calculates the exact IC50 and apparent IC50 in closed form.

**Interactive simulation tools (R Shiny).** Four standalone applications in
``shiny-apps/``:

==========================  =================================================================
Application                 Description
==========================  =================================================================
``kd-sim``                  Monte Carlo simulation of saturation binding; *K*\ d fitted with
                            the exact quadratic equation and with the hyperbola
``kin-sim``                 Monte Carlo simulation of association and dissociation; *k*\ on
                            fitted with the exact second-order and the pseudo-first-order
                            equations
``comp-sim``                Monte Carlo simulation of competitive binding; apparent IC50,
                            *K*\ i by Cheng–Prusoff and Wang's group equation, and *K*\ i
                            fitted directly with Wang's exact equation
``time-sim``                Apparent *K*\ d of a saturation experiment read before
                            equilibrium, and the incubation time needed
==========================  =================================================================

Requirements
------------

- Python >= 3.7.6, Bokeh >= 3.4.0, NumPy >= 1.26.0
- R >= 4.3 with shiny >= 1.13.0, bslib >= 0.4.0, ggplot2 >= 3.5.0 and minpack.lm

Usage
-----

Visualization tools (the HTML pages are written to the current directory):

.. code-block:: bash

   pip install -r requirements.txt
   python binding_curve_viewer_v2.py

Simulation tools, from the R console:

.. code-block:: r

   install.packages(c("shiny", "bslib", "ggplot2", "minpack.lm"))
   shiny::runApp("shiny-apps/kd-sim")

On the web server, a Monte Carlo run is limited to 1000 simulated experiments.
In local mode, 5000, 10,000 or 100,000 are available:

.. code-block:: r

   Sys.setenv(BINDING_APP_LOCAL = "true")
   shiny::runApp("shiny-apps/kd-sim")

Documentation:

.. code-block:: bash

   cd docs
   pip install -r requirements.txt
   make html

Repository structure
--------------------

.. code-block:: text

   binding_curve_viewer_v2.py   visualization tools (Python/Bokeh)
   shiny-apps/                  interactive simulation tools (R Shiny)
   scripts/                     supporting scripts for the 2024 article
   docs/                        Sphinx documentation

Version 1.0
-----------

The version described in the 2024 article in *J. Chem. Inf. Model.* is
available under the tag
`v1.0.0 <https://github.com/ydu-sci/Binding_Curve_Viewer/tree/v1.0.0>`_.

How to cite
-----------

Du, Y. Binding Curve Viewer: Visualizing the Equilibrium and Kinetics of
Protein–Ligand Binding and Competitive Binding. *J. Chem. Inf. Model.* 2024,
64 (10), 4180–4192. https://doi.org/10.1021/acs.jcim.4c00130

License
-------

The code is released under the MIT License (see ``LICENSE``). The web pages are
made available under the CC BY 4.0 license.
