rule phytest:
    """
    Run phytest on an upload data and tree
    """
    message: "Running phytest"
    input:
        data = rules.extract_upload_metadata.output.metadata_csv,
        tree = rules.ladderize_tree.output.newick,
    output:
        report = report("{outdir}/{name}.phytest-report.html")
    conda:
        ENVS / "phytest.yaml"
    threads:
        int(config["phytest"].get("threads")) if config["phytest"].get("threads") else workflow.cores
    shell:
        """
        phytest {SCRIPTS}/phytests.py \
            --tree {input.tree} \
            --data {input.data} \
            --report {output.report} \
            -n {threads}
        """
