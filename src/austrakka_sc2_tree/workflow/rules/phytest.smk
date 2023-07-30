rule phytest:
    """
    Run phytest on an upload data and tree.

    :input data:               The CSV metadata file produced by the :smk:ref:`extract_upload_metadata` rule.

    :input tree:               The Newick file produced by the :smk:ref:`ladderize_tree` rule.

    :output report:            The output HTML report file that contains the results of the phytest.

    :conda:                    Path to the Conda environment file (phytest.yaml) in the ENVS directory.

    :threads:                  Number of threads to use. This value is retrieved from the config file under the 
                            "phytest" section. If no value is provided there, the total number of cores available in 
                            the workflow are used.

    .. note::
        This rule runs the phytest using a Python script on the uploaded data and tree. The report of this test is 
        stored in an HTML file. The number of threads used for this operation can be configured in the workflow config 
        file. If no value is provided, the total number of cores available in the workflow are used. 
    """
    message: "Running phytest"
    input:
        data = rules.extract_upload_metadata.output.metadata_csv,
        tree = rules.ladderize_tree.output.newick,
    output:
        report = report("{outdir}/{name}.phytest-report.html")
    params:
        error_report = "{outdir}/{name}.error.phytest-report.html"
    conda:
        ENVS / "phytest.yaml"
    threads:
        int(config["phytest"].get("threads")) if config["phytest"].get("threads") else workflow.cores
    shell:
        """
        phytest {SCRIPTS}/phytests.py \
            -v \
            --tree {input.tree} \
            --data {input.data} \
            --report {params.error_report} \
            -n {threads}
        # prevent error report from being deleted
        mv {params.error_report} {output.report}
        """
