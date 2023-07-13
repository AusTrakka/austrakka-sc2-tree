import pandas as pd

df = pd.read_csv(snakemake.input.nextclade_tsv, sep='\t')
df = df.rename(columns={
    "seqName": "Seq_ID",
    "clade": "Clade",
    "Nextclade_pango": "Lineage",
    "coverage": "Coverage",
    "qc.overallScore": "NC_qc.overallScore",
    "qc.overallStatus": "NC_qc.overallStatus",
})
df.drop(columns=["index"], inplace=True)
df.to_csv(snakemake.output.formated_tsv, sep='\t', index=False)