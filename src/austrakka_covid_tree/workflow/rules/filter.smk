rule filter:
    input:
        alignment=rules.nextclade.output.alignment,
        at_matadata_tsv=rules.rename_columns.output.at_matadata_tsv,
    output:
        alignment_filtered = temp("{group}.filtered.afa"),
        metadata_filtered = temp("{group}.filtered.tsv"),
    params:
        id_column = "Seq_ID",
        query = " & ".join(config['filter'])
    conda:
        ENVS / "nextstrain.yaml"
    log:
        LOGS / "filter" / "filter.{group}.log"
    shell:
        '''
        augur filter \
            --sequences {input.alignment} \
            --metadata {input.at_matadata_tsv} \
            --metadata-id-columns {params.id_column} \
            --query "{params.query}" \
            --output {output.alignment_filtered} \
            --output-metadata {output.metadata_filtered}  2>&1 | tee {log}
        '''