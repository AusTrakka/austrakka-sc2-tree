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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import os
import sys
from pathlib import Path

sys.path.append(os.path.abspath("../../src/austrakka_sc2_tree"))
from austrakka_sc2_tree import __about__ # noqa: E402

# -- Project information -----------------------------------------------------

project = "Austrakka Covid Tree"
copyright = "2023-present, Austrakka"
author = "Austrakka"

# The full version, including alpha/beta/rc tags
release = __about__.__version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["snakedoc", "sphinx_immaterial"]

smk_linkcode_mapping = (str(Path(__file__).absolute().parents[2]), "https://github.com/AusTrakka/austrakka-sc2-tree/blob/main")

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_immaterial"
# html_logo = "img/snakedoc-logo.svg"
html_theme_options = {
    "toc_title_is_page_title": True,
    "site_url": "https://austrakka.github.io/austrakka-sc2-tree/",
    "repo_url": "https://github.com/AusTrakka/austrakka-sc2-tree",
    "repo_name": "austrakka-sc2-tree",
    "features": ["content.code.annotate", "navigation.instant", "toc.follow"],
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
