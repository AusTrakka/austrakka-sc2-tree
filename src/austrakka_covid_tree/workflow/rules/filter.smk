rule filter:
    """
    Filters the input alignment and metadata files based on a specific query.
    
    .. note::
        The exact filter query is specified in the configuration file (config['filter']).
        Ensure that it is a valid query for the 'augur filter' command.

    :input alignment:         the aligned fasta file output from the previous 'nextclade' rule.
    :input at_matadata_tsv:   the metadata file with renamed columns from the previous 'rename_columns' rule.

    :output alignment_filtered:   the alignment file after applying the filter query.
    :output metadata_filtered:    the metadata file after applying the filter query.

    :params id_column:            the column in the metadata file used to identify sequences.
    :params query:                the filtering query to apply. This query is constructed by joining the elements of the 'filter' config with " & ".

    :conda:                       the environment used for this rule, specified by the "nextstrain.yaml" file in the "ENVS" directory.

    :log:                         the log file for this rule, located in the "LOGS/filter" directory. The log file is named "filter.{group}.log", where "{group}" is a wildcard representing the group of the sequences being processed.

    The rule is implemented by the 'augur filter' command, which takes the input alignment and metadata files, 
    applies the specified filtering query, and outputs the filtered alignment and metadata files.
    """
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