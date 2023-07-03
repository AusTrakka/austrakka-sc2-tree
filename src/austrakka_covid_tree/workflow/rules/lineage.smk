rule download_nextclade_db:
    """
    Download latest nextclade covid db
    """
    output:
        temp(directory("nextclade_data_dir"))
    conda:
        ENVS / "nextclade.yaml"
    log:
        LOGS / "lineage" / "download_nextclade_db.log"
    params:
        download_unreleased_tree=1 if config.get('use_unreleased_nextclade', False) else 0
    shell:
        """
        nextclade dataset get \
            --verbose \
            --name 'sars-cov-2' \
            --output-dir {output} > {log}
        if [ {params.download_unreleased_tree} = 1 ]
        then
            echo "Using unreleased Nextclade tree!" >> {log}
            wget https://nextstrain.org/charon/getDataset?prefix=staging/nextclade/sars-cov-2 -O {output}/tree.json
        fi
        """

rule nextclade:
    """
    Run nextclade to get the lineage calls
    """
    input:
        fasta = config["fasta"],
        nextclade_data_dir = "nextclade_data_dir"
    output:
        nextclade_tsv=temp("{group}.nextclade.tsv"),
        alignment=temp("{group}.nextclade.afa")
    threads:
        config["lineage"].get("threads") if config["lineage"].get("threads") else workflow.cores
    conda:
        ENVS / "nextclade.yaml"
    log:
        LOGS / "lineage" / "nextclade.{group}.log"
    shell:
        """
        nextclade run \
            -j {threads} \
            -D {input.nextclade_data_dir} \
            --output-tsv {output.nextclade_tsv} \
            --output-fasta {output.alignment} \
            {input.fasta}  2>&1 | tee {log}
        """

rule collapse_lineages:
    """
    Use pango-collapse to collapse nextclade lineages 
    """
    input:
        nextclade_tsv=rules.nextclade.output.nextclade_tsv
    output:
        nextclade_collapsed_tsv=temp("{group}.nextclade.collapsed.tsv")
    params:
        url=config["lineage"].get("pango_collapse_url", "https://raw.githubusercontent.com/MDU-PHL/pango-collapse/main/pango_collapse/collapse.txt")
    conda:
        ENVS / "pango_collapse.yaml"
    log:
        LOGS / "lineage" / "collapse_lineages.{group}.log"
    shell:
        """
        pango-collapse -l Nextclade_pango --latest --url {params.url} -o {output} {input} 2>&1 | tee {log}
        """

rule mask_lineages:
    """
    Set Unassigned and LowCoverage
    """
    input:
        nextclade_tsv=rules.collapse_lineages.output.nextclade_collapsed_tsv
    output:
        masked_nextclade_tsv=temp("{group}.nextclade.masked.tsv")
    conda:
        ENVS / "python.yaml"
    log:
        LOGS / "lineage" / "mask_lineages.{group}.log"
    script:
        SCRIPTS / "mask_lineages.py"

rule rename_columns:
    """
    Rename Nextclade columns to AT format
    """
    input:
        nextclade_tsv=rules.mask_lineages.output.masked_nextclade_tsv
    output:
        at_matadata_tsv="{group}.metadata.tsv",
    conda:
        ENVS / "python.yaml"
    script:
        SCRIPTS / "rename_columns.py"
