========
Overview
========

Binding Curve Viewer is a set of tools to visualize and simulate the equilibrium and kinetics 
of binding reactions with 1:1 stoichiometry and no cooperativity for experiment planning and 
validation. To use the tools online, please visit the `website <https://www.eplatton.net/binding-curve-viewer/>`_. Source code can be found at `GitHub <https://github.com/ydu-sci/Binding_Curve_Viewer>`_.

Requirements
============

1. Visualization Tools (Python Environment)
--------------------------------------------

To generate the visualization tools of Binding Curve Viewer, the following dependencies are required:

#. Python >= 3.7.6
#. Bokeh >= 3.4.0
#. Numpy >= 1.26.0

2. Interactive Simulation Tools (R Environment)
-----------------------------------------------

To run the interactive simulation tools, the following dependencies are required:

#. R >= 4.3
#. shiny >= 1.13.0
#. bslib >= 0.4.0
#. ggplot2 >= 3.5.0
#. minpack.lm (used by all apps except Kinetics of Association and Dissociation)

The packages can be installed in R with::

   install.packages(c("shiny", "bslib", "ggplot2", "minpack.lm"))

Usage
=====

Visit the `website <https://www.eplatton.net/binding-curve-viewer/>`_ to use the online tools.
Or download the source code from `GitHub Page <https://github.com/ydu-sci/Binding_Curve_Viewer/>`_
and run binding_curve_viewer_v2.py to get all the visualization tools of Binding Curve Viewer. 
To use the interactive simulation components, run the corresponding scripts in R environment. 

.. _bcv_paper:

How to Cite
===========

Du, Y. 2024. Binding Curve Viewer: Visualizing the Equilibrium and Kinetics of Protein–Ligand Binding and Competitive Binding. J. Chem. Inf. Model. 64(10):4180-4192. [`Free Full Text from Publisher <https://pubs.acs.org/doi/10.1021/acs.jcim.4c00130>`_] [`PubMed <https://pubmed.ncbi.nlm.nih.gov/38720179/>`_]

License of the code
===================

    MIT License

    Copyright (c) 2023-2026 Yu Du

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

License of the webpages
=======================

`CC BY 4.0 <https://creativecommons.org/licenses/by/4.0/>`_