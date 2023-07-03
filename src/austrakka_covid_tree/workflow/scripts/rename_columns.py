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
df.to_csv(snakemake.output.at_matadata_tsv, sep='\t', index=False)
