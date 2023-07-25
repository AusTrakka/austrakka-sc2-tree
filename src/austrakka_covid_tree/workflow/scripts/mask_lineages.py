#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv(snakemake.input.nextclade_tsv, sep="\t")

df.coverage = df.coverage.fillna(0)
df.coverage = df.coverage * 100
df.loc[df.coverage < 90, "Nextclade_pango"] = "LowCoverage"

df["qc.overallStatus"] = df["qc.overallStatus"].fillna("bad")
df.loc[df["qc.overallStatus"] == "bad", "Nextclade_pango"] = "Unassigned"
df.to_csv(snakemake.output.masked_nextclade_tsv, sep="\t", index=False)
