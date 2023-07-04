.. _usage:

Usage
============

To run the Snakemake pipeline, use the following command:

.. code-block:: bash

    austrakka-covid-tree run [OPTIONS] [TARGET]

All unrecognized arguments are passed onto Snakemake. If the TARGET file is not specified, the pipeline will run the 'all' rule by default.

Options
-------

- ``target [TARGET]``: Specify the file to generate. If None will run the pipeline 'all' rule.

- ``--fasta [PATH]``: (Required) Path to the fasta file to build tree from.

- ``--outdir [PATH]``: Output directory for tree and metadata.

- ``--dated / --no-dated``: Option to include timestamp in the output.

- ``--starting-tree [PATH]``: Starting tree to use for tree building.

- ``--filter [TEXT]``: Filters to apply to the sequences before tree building. Default is 'Coverage >= 90'.

- ``--tree-threads [INTEGER]``: Number of threads to use for tree building.

- ``--lineage-pango-collapse-url [TEXT]``: URL to pango collapse file. Default is `https://raw.githubusercontent.com/MDU-PHL/pango-collapse/main/pango_collapse/collapse.txt`.

- ``--lineage-threads [INTEGER]``: Number of threads to use for lineage calling.

- ``--lineage-unreleased / --no-lineage-unreleased``: Option to use unreleased nextclade version.

- ``--config [FILE]``: Path to snakemake config file. Overrides existing config and defaults.

- ``--resource -r [PATH]``: Additional resources to copy to workdir at run time (relative to pipeline directory).

- ``--profile -p [TEXT]``: Name of profile to use for configuring Snakemake.

- ``--force -f``: Force the execution of pipeline regardless of already created output.

- ``--lock -l``: Lock the working directory.

- ``--keep-resources -R``: Keep resources after pipeline completes.

- ``--keep-snakemake -S``: Keep .snakemake folder after pipeline completes.

- ``--dag -d [PATH]``: Save directed acyclic graph to file. Must end in .pdf, .png or .svg.

- ``--cores -c [INTEGER]``: Set the number of cores to use. If None will use all cores.

- ``--verbose -v``: Run pipeline in verbose mode.

- ``--help-snakemake -hs``: Print the snakemake help and exit.

- ``--help -h``: Show this message and exit.

For more information, refer to the austrakka-covid-tree CLI help or the Snakemake documentation.
