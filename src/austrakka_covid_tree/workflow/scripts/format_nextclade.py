import pandas as pd

COVERGE_THRESHOLD = 0.9

df = pd.read_csv(snakemake.input.nextclade_tsv, sep="\t")

# rename columns for AT proforma
df = df.rename(columns={
    "seqName": "Seq_ID",
    "clade": "Clade",
    "Nextclade_pango": "Lineage",
    "coverage": "Coverage",
    "qc.overallScore": "NC_qc.overallScore",
    "qc.overallStatus": "NC_qc.overallStatus",
})

# Set df QC to FAIL if coverage is less than COVERGE_THRESHOLD otherwise set to PASS
df["QC"] = df.apply(lambda x: "FAIL" if x["Coverage"] < COVERGE_THRESHOLD else "PASS", axis=1)
df.drop(columns=["index"], inplace=True)
df.to_csv(snakemake.output.formated_tsv, sep="\t", index=False)
