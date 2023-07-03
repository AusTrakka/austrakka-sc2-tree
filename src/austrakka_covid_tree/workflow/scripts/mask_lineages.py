#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv(snakemake.input.nextclade_tsv, sep="\t")
df.loc[df["qc.overallStatus"] == "bad", "Nextclade_pango"] = "Unassigned"
df.coverage = df.coverage * 100
df.loc[df.coverage < 90, "Nextclade_pango"] = "LowCoverage"
df.to_csv(snakemake.output.masked_nextclade_tsv, sep="\t", index=False)
