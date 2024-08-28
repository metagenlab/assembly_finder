import pandas as pd
import numpy as np
import json


def read_json(file):
    try:
        return pd.json_normalize(json.load(open(file)), record_path=["reports"])
    except KeyError:
        return pd.read_json(file)


# Read tables
summary_df = read_json(snakemake.input.summary)
lineage_df = pd.read_csv(snakemake.input.lineage, sep="\t")
# Params
rank = snakemake.params.rank
nrank = snakemake.params.nrank
taxon = snakemake.params.taxon

# format summary column names
summary_df.columns = (
    summary_df.columns.str.replace("assembly_info.", "")
    .str.replace("assembly_stats.", "")
    .str.replace("organism.", "")
)


# Merge lineage and genome summary
df = summary_df.merge(lineage_df, on="tax_id")
df = df.replace(np.nan, "na")
df = df[~df["accession"].str.contains(" ")]
# sort according to refseq category and assembly level
cols = ["refseq_category", "assembly_level", "checkm_version", "busco.busco_ver"]
sort = []
if "refseq_category" in df.columns:
    df["refseq_category"] = pd.Categorical(
        df["refseq_category"],
        ["reference genome", "representative genome", "na"],
        ordered=True,
    )
    sort.append("refseq_category")
if "assembly_level" in df.columns:
    df["assembly_level"] = pd.Categorical(
        df["assembly_level"],
        ["Complete Genome", "Chromosome", "Scaffold", "Contig", "na"],
        ordered=True,
    )
    sort.append("assembly_level")
if "checkm_info.completeness" in df.columns:
    sort.append("checkm_info.completeness")
if "busco.complete" in df.columns:
    sort.append("busco.complete")

df.sort_values(
    sort,
    inplace=True,
)

if rank and nrank:
    df = df.groupby(rank).head(nrank)

tax_cols = [
    "accession",
    "tax_id",
    "name",
    "rank",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

if taxon:
    tax_cols.insert(0, "taxon")

df[tax_cols].to_csv(snakemake.output.tax, sep="\t", index=None)
df[summary_df.columns].to_csv(snakemake.output.gen, sep="\t", index=None)
df[["accession"]].drop_duplicates().to_csv(
    snakemake.output.acc, sep="\t", index=None, header=False
)
