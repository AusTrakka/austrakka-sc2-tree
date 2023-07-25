import pandas as pd

COVERGE_THRESHOLD = 0.9

df = pd.read_csv(snakemake.input.nextclade_tsv, sep="\t")

# rename columns for AT proforma
df = df.rename(
    columns={
        "seqName": "Seq_ID",
        "Nextclade_pango": "Lineage",
        "coverage": "Coverage",
    }
)

# Set df QC to FAIL if coverage is less than COVERGE_THRESHOLD otherwise set to PASS
df["QC"] = df.apply(lambda x: "FAIL" if x["Coverage"] < COVERGE_THRESHOLD else "PASS", axis=1)

# Extract proforma
upload_df = df[["Seq_ID", "Coverage", "Lineage", "Lineage_family", "Lineage_full", "Lineage_note", "QC"]]
upload_df.to_csv(snakemake.output.metadata_csv, index=False)
