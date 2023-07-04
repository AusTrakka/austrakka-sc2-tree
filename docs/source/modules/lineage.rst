:tocdepth: 3

.. _lineage:

Lineage
==========

This module contains the rules to assign covid lineages.

Rules
-----
.. smk:autodoc:: ../../src/austrakka_covid_tree/workflow/Snakefile download_nextclade_db nextclade collapse_lineages mask_lineages rename_columns
  :configfile: config.yaml