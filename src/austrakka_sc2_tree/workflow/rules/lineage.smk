rule download_nextclade_db:
    """
    Downloads the latest Nextclade database for SARS-CoV-2.

    :output:                   Temporary directory to store the downloaded Nextclade database.

    :conda:                    Path to the Conda environment file (nextclade.yaml) in the ENVS directory.

    :param download_unreleased_tree: Flag to determine whether to download the unreleased Nextclade tree. Set to 1 if 
                                     'use_unreleased_nextclade' is set to True in the configuration, otherwise 0.

    .. note::
        This rule uses the 'nextclade dataset get' command to download the latest Nextclade database for SARS-CoV-2. 
        If 'use_unreleased_nextclade' is set to True in the configuration, it also downloads the unreleased Nextclade 
        tree from the Nextstrain staging area.
    """
    output:
        nextclade_data_dir=temp(directory("{outdir}/nextclade_data_dir"))
    conda:
        ENVS / "nextclade.yaml"
    params:
        download_unreleased_tree=1 if config.get('use_unreleased_nextclade', False) else 0
    shell:
        """
        nextclade dataset get \
            --verbose \
            --name 'nextstrain/sars-cov-2/wuhan-hu-1/orfs' \
            --output-dir {output.nextclade_data_dir}
        if [ {params.download_unreleased_tree} = 1 ]
        then
            echo "Using unreleased Nextclade tree!"
            wget https://nextstrain.org/charon/getDataset?prefix=staging/nextclade/sars-cov-2 -O {output}/tree.json
        fi
        """

rule nextclade:
    """
    Runs Nextclade to get the lineage calls for the input genome sequences.

    :input fasta:              Path to the input FASTA file containing genome sequences.
    :input nextclade_data_dir: Directory containing the downloaded Nextclade database.

    :output nextclade_tsv:     Temporary TSV file to store the lineage call results from Nextclade.
    :output alignment:         Temporary FASTA file to store the aligned sequences generated by Nextclade.

    :threads:                  Number of threads to use for this rule, which is fetched from the "lineage" section of the 
                               configuration. If not specified, the total number of available cores in the workflow 
                               is used.

    :conda:                    Path to the Conda environment file (nextclade.yaml) in the ENVS directory.

    .. note::
        This rule uses the 'nextclade run' command to call lineages for the input genome sequences. It uses the downloaded 
        Nextclade database for this purpose. The lineage calls are stored in a TSV file, and the aligned sequences are 
        stored in a FASTA file. The standard error and output from the command are logged in the specified log file.
    """
    input:
        fasta = "{outdir}/{name}.filtered.fasta" if config.get("data", None) else config["fasta"],
        nextclade_data_dir = rules.download_nextclade_db.output.nextclade_data_dir
    output:
        nextclade_tsv=temp("{outdir}/{name}.nextclade.raw.tsv"),
        alignment=temp("{outdir}/{name}.nextclade.afa")
    params:
        reference_sequence = RESOURCES / "MN908947.3.fna",
    threads:
        config["lineage"].get("threads") if config["lineage"].get("threads") else workflow.cores
    conda:
        ENVS / "nextclade.yaml"
    shell:
        """
        nextclade run \
            -j {threads} \
            -D {input.nextclade_data_dir} \
            --input-ref {params.reference_sequence} \
            --output-tsv {output.nextclade_tsv} \
            --output-fasta {output.alignment} \
            {input.fasta}

        # add nextclade version to the output file
        nextclade_data=$(grep '"tag":' {input.nextclade_data_dir}/pathogen.json | awk -F'"' '{{ print $4 }}')
        nextclade_version=$(nextclade -V)
        awk -v OFS='\t' -v val="$nextclade_version;nextclade_data $nextclade_data" \
            '{{if(NR==1) print $0, "Lineage_note"; else print $0, val}}' {output.nextclade_tsv} > {output.nextclade_tsv}.tmp
        mv {output.nextclade_tsv}.tmp {output.nextclade_tsv}
        """


rule mask_lineages:
    """
    Masks 'Unassigned' and 'LowCoverage' lineage calls in the input lineage call results.

    :input nextclade_tsv:      The TSV file produced by the :smk:ref:`collapse_lineages` rule, containing collapsed lineage call results.

    :output masked_nextclade_tsv: Temporary TSV file to store the masked lineage call results.

    :conda:                    Path to the Conda environment file (python.yaml) in the ENVS directory.

    :script:                   Python script (mask_lineages.py) in the SCRIPTS directory that is run to mask 'Unassigned' 
                               and 'LowCoverage' lineage calls.

    .. note::
        This rule runs a Python script to mask 'Unassigned' and 'LowCoverage' lineage calls in the input lineage call results. 
        The masked lineage calls are stored in a TSV file. The standard error and output from the command are logged in the 
        specified log file.
    """
    input:
        nextclade_tsv=rules.nextclade.output.nextclade_tsv
    output:
        masked_nextclade_tsv=temp("{outdir}/{name}.nextclade.masked.tsv")
    params:
        mask_low_coverage=config["lineage"]['mask']['low_coverage'] == "True",
        mask_unassigned=config["lineage"]['mask']['unassigned'] == "True",
    conda:
        ENVS / "python.yaml"
    script:
        SCRIPTS / "mask_lineages.py"


rule collapse_lineages:
    """
    Runs pango-collapse to collapse Nextclade lineages.

    :input nextclade_tsv:      The TSV file produced by the :smk:ref:`nextclade` rule, containing lineage call results.

    :output nextclade_collapsed_tsv: Temporary TSV file to store the collapsed lineage results.

    :param url:                The URL of the collapse.txt file used for lineage collapsing. This is fetched from the 
                               "lineage" section of the configuration, with a default value pointing to the collapse.txt 
                               file in the MDU-PHL/pango-collapse repository on GitHub.

    :conda:                    Path to the Conda environment file (pango_collapse.yaml) in the ENVS directory.


                               structured as "collapse_lineages.{name}.log".

    .. note::
        This rule uses the 'pango-collapse' command to collapse lineages in the input lineage call results. The collapsed 
        lineage calls are stored in a TSV file. The standard error and output from the command are logged in the 
        specified log file.
    """
    input:
        nextclade_tsv=rules.mask_lineages.output.masked_nextclade_tsv
    output:
        nextclade_collapsed_tsv=temp("{outdir}/{name}.nextclade.collapsed.tsv")
    params:
        file=config["lineage"].get("pango_collapse_file")
    conda:
        ENVS / "pango_collapse.yaml"
    shell:
        """
        if [[ "{params.file}" == http* ]]; then
            # If the file is a URL, run pango-collapse with the --url option
            pango-collapse -l Nextclade_pango --latest --url "{params.file}" -o "{output}" "{input}"
        else
            # If the file is a local path, run pango-collapse with the --collapse-file option
            pango-collapse -l Nextclade_pango --collapse-file "{params.file}" -o "{output}" "{input}"
        fi
        """

rule extract_upload_metadata:
    """
    Nextclade output to AT upload format.

    :input nextclade_tsv:       The TSV file produced by the :smk:ref:`mask_lineages` rule, containing masked lineage call results.

    :output at_matadata_tsv:    The output TSV file with column names in AT format.

    :conda:                     Path to the Conda environment file (python.yaml) in the ENVS directory.

    :script:                    Python script (rename_columns.py) in the SCRIPTS directory that is run to rename Nextclade 
                                columns to AT format.

    .. note::
        This rule runs a Python script to rename Nextclade columns to the AT format in the input lineage call results. 
        The renamed lineage calls are stored in a TSV file. The standard error and output from the command are logged in the 
        specified log file.
    """
    input:
        nextclade_tsv=rules.collapse_lineages.output.nextclade_collapsed_tsv
    output:
        metadata_csv="{outdir}/{name}.metadata.csv",
    conda:
        ENVS / "python.yaml"
    script:
        SCRIPTS / "extract_upload_metadata.py"
