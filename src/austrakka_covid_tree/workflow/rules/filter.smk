# rule filter_metadata:
#     """
#     Filters the in fasta and metadata files based on a specific query.
#     """
#     input:
#         fasta = config["fasta"],
#         metadata = config["data"],
#     output:
#         fasta_filtered = temp("{outdir}/{name}.filtered.fasta"),
#         metadata_filtered = temp("{outdir}/{name}.filtered.tsv"),
#     params:
#         id_column = "Seq_ID",
#         query = " & ".join(config['filter'].get("metadata"))
#     conda:
#         ENVS / "nextstrain.yaml"
#     shell:
#         '''
#         augur filter \
#             --sequences {input.fasta} \
#             --metadata {input.metadata} \
#             --metadata-id-columns {params.id_column} \
#             --query "{params.query}" \
#             --output {output.alignment_filtered} \
#             --output-metadata {output.metadata_filtered}
#         '''

rule filter_nextclade:
    """
    Filters the nextcalde alignment and linage files based on a specific query.
    
    .. note::
        The exact filter query is specified in the configuration file (config['filter']).
        Ensure that it is a valid query for the 'augur filter' command.

    :input alignment:         the aligned fasta file output from the previous 'nextclade' rule.
    :input at_matadata_tsv:   the nextclade file with renamed columns from the previous 'rename_columns' rule.

    :output alignment_filtered:   the alignment file after applying the filter query.
    :output nextclade_filtered:    the nextclade file after applying the filter query.

    :params id_column:            the column in the nextclade file used to identify sequences.
    :params query:                the filtering query to apply. This query is constructed by joining the elements of the 'filter' config with " & ".

    :conda:                       the environment used for this rule, specified by the "nextstrain.yaml" file in the "ENVS" directory.

    The rule is implemented by the 'augur filter' command, which takes the input alignment and nextclade files, 
    applies the specified filtering query, and outputs the filtered alignment and nextclade files.
    """
    input:
        alignment=rules.nextclade.output.alignment,
        metadata=rules.format_nextclade.output.formated_tsv,
    output:
        alignment_filtered = temp("{outdir}/{name}.filtered.afa"),
        metadata_filtered = temp("{outdir}/{name}.nextclade.filtered.tsv"),
    params:
        id_column = "Seq_ID",
        query = " & ".join(config['filter'].get("nextclade"))
    conda:
        ENVS / "nextstrain.yaml"
    shell:
        '''
        augur filter \
            --sequences {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.id_column} \
            --query "{params.query}" \
            --output {output.alignment_filtered} \
            --output-metadata {output.metadata_filtered}
        '''