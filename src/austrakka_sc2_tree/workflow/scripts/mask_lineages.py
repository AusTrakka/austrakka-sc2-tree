#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv(snakemake.input.nextclade_tsv, sep="\t")

df.coverage = df.coverage.fillna(0)
df.coverage = df.coverage * 100
if snakemake.params.mask_low_coverage:
    print(type(snakemake.params.mask_low_coverage))
    print("Masking low coverage samples")
    df.loc[df.coverage < 90, "Nextclade_pango"] = "LowCoverage"
if snakemake.params.mask_unassigned:
    print("Masking unassigned samples")
    df["qc.overallStatus"] = df["qc.overallStatus"].fillna("bad")
    df.loc[df["qc.overallStatus"] == "bad", "Nextclade_pango"] = "Unassigned"
df.to_csv(snakemake.output.masked_nextclade_tsv, sep="\t", index=False)
