rule filiter_csv:
    input:
        csv_file = "data/iris.csv"
    output:
        iris = "data/iris.csv"
    conda:
        ENVS / "csvtk.yaml"
    shell:
        """
        JUST_GOOD := \
            csvtk grep -f Species_expected -p 'SARS-CoV-2' \
        | csvtk grep -f Seq_QC -p 'PASS' \
        | csvtk grep -f Seq_category -p 'SEQ' \
        | csvtk grep -f Case_status -v -r -p 'eject'

        JUST_AU := csvtk grep -f Country -p Australia
        """