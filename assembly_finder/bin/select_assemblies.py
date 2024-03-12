import json
import pandas as pd
import numpy as np

with open(snakemake.input[0]) as file:
    data = json.load(file)
df = pd.json_normalize(data, record_path=["reports"])
df = df.replace(np.nan, "na")
df["entry"] = [wildcards.entry] * len(df)

# sort according to refseq category and assembly level
sort = []
try:
    df["assembly_info.refseq_category"] = pd.Categorical(
        df["assembly_info.refseq_category"],
        ["reference genome", "representative genome", "na"],
        ordered=True,
    )
    sort.append("assembly_info.refseq_category")
except KeyError:
    df["assembly_info.assembly_level"] = pd.Categorical(
        df["assembly_info.assembly_level"],
        ["Complete Genome", "Chromosome", "Scaffold", "Contig", "na"],
        ordered=True,
    )
    sort.append("assembly_info.assembly_level")

df.sort_values(
    sort,
    inplace=True,
)
df.iloc[[0]].to_csv(output[0], sep="\t", index=None)
