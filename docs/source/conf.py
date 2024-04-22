# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))


# -- Project information -----------------------------------------------------

project = 'binding_curve_viewer'
copyright = '2023, Yu Du'
author = 'Yu Du'
html_title = "Binding Curve Viewer's Documentation"

# The full version, including alpha/beta/rc tags
release = '1.0.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'guzzle_sphinx_theme',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
import guzzle_sphinx_theme
html_static_path = ['_static']
html_translator_class = 'guzzle_sphinx_theme.HTMLTranslator'
html_theme_path = guzzle_sphinx_theme.html_theme_path()
html_theme = 'guzzle_sphinx_theme'
# Custom sidebar templates, maps document names to template names.
html_sidebars = {
    '**': ['logo-text.html', 'globaltoc.html', 'searchbox.html']
}
# Guzzle theme options (see theme.conf for more information)
html_theme_options = {
    # Set the name of the project to appear in the sidebar
    'project_nav_name': 'Binding Curve Viewer',
    #'globaltoc_maxdepth': 1,

}

latex_elements = {'pointsize': '10pt'}