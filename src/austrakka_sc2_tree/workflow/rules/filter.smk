import datetime
today = datetime.date.today()
max_date = config.get("max_date", None)
days_ago = config.get("days_ago", None)
if days_ago:
    min_date = (today - datetime.timedelta(days=days_ago)).strftime("%Y-%m-%d")
else:
    min_date = config.get("min_date", None)


rule filter_metadata:
    """
    Filters the input fasta and metadata files based on a specific query.

    .. note::
        The exact filter query is specified in the configuration file (config['filter']).
        Ensure that it is a valid query for the 'augur filter' command.

    :input fasta:              the fasta file to be filtered.
    :input metadata:           the metadata file to be filtered.

    :output fasta_filtered:    the fasta file after applying the filter query.
    :output metadata_filtered: the metadata file after applying the filter query.
    :output tmp_metadata:      a temporary metadata file with renamed date column.

    :params id_column:         the column in the metadata file used to identify sequences.
    :params query:             the filtering query to apply. This query is constructed by joining the elements of the 'filter' config with " & ".
    :params min_date_flag:     optional flag for specifying minimum date for filtering. If not set, it's ignored.
    :params max_date_flag:     optional flag for specifying maximum date for filtering. If not set, it's ignored.

    :conda:                    the environment used for this rule, specified by the "nextstrain.yaml" file in the "ENVS" directory.

    This rule uses the 'augur filter' command to apply the specified filtering query to the input fasta and metadata files,
    and outputs the filtered fasta and metadata files. The 'sed' command is used to rename the 'Date_coll' column to 'date' in the metadata file.
    """
    input:
        fasta = config["fasta"],
        metadata = config["data"] if config["data"] else ""
    output:
        fasta_filtered = temp("{outdir}/{name}.filtered.fasta"),
        metadata_filtered = temp("{outdir}/{name}.metadata.filtered.tsv"),
        tmp_metadata = temp("{outdir}/{name}.metadata.tmp")
    params:
        id_column = "Seq_ID",
        query = " & ".join(config['filter'].get("data", [])),
        min_date_flag = f"--min-date {min_date}" if min_date else "",
        max_date_flag = f"--max-date {max_date}" if max_date else ""
    conda:
        ENVS / "nextstrain.yaml"
    shell:
        '''
        # rename date col to Date_coll for augur filter
        sed '1s/Date_coll/date/' {input.metadata} > {output.tmp_metadata}
        augur filter \
            --sequences {input.fasta} \
            --metadata {output.tmp_metadata} \
            --metadata-id-columns {params.id_column} \
            --query "{params.query}" \
            --output {output.fasta_filtered} \
            --output-metadata {output.metadata_filtered} \
            {params.min_date_flag} {params.max_date_flag}
        '''

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
        metadata=rules.extract_upload_metadata.output.metadata_csv,
    output:
        alignment_filtered = temp("{outdir}/{name}.filtered.afa"),
        nextclade_filtered = temp("{outdir}/{name}.nextclade.filtered.tsv"),
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
            --output-metadata {output.nextclade_filtered}
        '''