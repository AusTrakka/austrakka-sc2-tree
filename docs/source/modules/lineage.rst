:tocdepth: 3

.. _lineage:

Lineage
==========

This module contains the rules to assign covid lineages.

Rules
-----
.. smk:autodoc:: ../../src/austrakka_sc2_tree/workflow/Snakefile download_nextclade_db nextclade collapse_lineages mask_lineages format_nextclade
  :configfile: config.yaml