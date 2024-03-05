import pandas as pd
import numpy as np
import os
import glob
import sys
from ete3 import NCBITaxa
import json

# Download params
download = config["download"]
# Path params
outdir = config["outdir"]
taxdump = config["taxdump"]

# Assemblies params
inp = config["input"]
nb = config["nb"]
db = config["db"]

accession = "refseq_seq_accession"
if db == "genbank":
    accession = "refseq_seq_accession"

uid = config["uid"]
alvl = config["alvl"]
rcat = config["rcat"]
excl = config["excl"]
annot = config["annot"]
rank = config["rank"]
nrank = config["nrank"]

# File extensions to download
sfxs = config["sfxs"]

# Entry table colnames
colnames = [
    "entry",
    "nb",
    "db",
    "uid",
    "alvl",
    "rcat",
    "excl",
    "annot",
    "rank",
    "nrank",
]

# Entry table values
values = [inp, nb, db, uid, alvl, rcat, excl, annot, rank, nrank]
values = [str(value).split(",") for value in values]

# Get parameters keys and values
param_keys = colnames[1:]
param_values = values[1:]

# Default param value is either the first entry param or the default one
# example: if nb=['1'], then 1 is the value for all entries

empty_to_val = {}
for key, val in zip(param_keys, param_values):
    if len(val) == 1:
        empty_to_val.update({key: {np.nan: val[0]}})
    else:
        empty_to_val.update({key: {np.nan: config[key]}})

# Check if input is file
if os.path.isfile(inp):
    entry_dt = pd.read_csv(inp, sep="\t").to_dict()
    # replace empty values with default ones or first entry for the param
    entry_dt = empty_to_val | entry_dt
    entry_df = pd.DataFrame.from_dict(entry_dt).dropna(subset="entry")


# If not create the dataframe
else:
    val_dic = dict(zip(colnames, values))
    entry_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in val_dic.items()]))

# Replace empty values with default params
# Drop empty entries
entry_df = entry_df.replace(
    empty_to_val,
).dropna()
entry_df["nb"] = [int(nb) if isinstance(nb, int) else nb for nb in entry_df["nb"]]
entry_df["entry"] = [
    int(entry) if isinstance(entry, float) else entry for entry in entry_df["entry"]
]

# Set entry as index
entry_df = entry_df.astype({"entry": str}).set_index("entry")

# Get entry list
entries = list(entry_df.index)


rule download_taxdump:
    output:
        temp(os.path.join(taxdump, "taxdump.tar.gz")),
    log:
        os.path.join(outdir, "logs", "taxdump-download.log"),
    retries: 2
    shell:
        """
        curl -o {output} https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz &> {log} 
        """


rule extract_taxdump:
    input:
        os.path.join(taxdump, "taxdump.tar.gz"),
    output:
        os.path.join(taxdump, "citations.dmp"),
        os.path.join(taxdump, "delnodes.dmp"),
        os.path.join(taxdump, "division.dmp"),
        os.path.join(taxdump, "gencode.dmp"),
        os.path.join(taxdump, "images.dmp"),
        os.path.join(taxdump, "merged.dmp"),
        os.path.join(taxdump, "names.dmp"),
        os.path.join(taxdump, "nodes.dmp"),
    params:
        taxdump,
    log:
        os.path.join(outdir, "logs", "taxdump.log"),
    shell:
        """
        tar -xzvf {input} -C {params}
        """


rule get_assembly_summaries:
    output:
        temp(os.path.join(outdir, "json", "{entry}.json")),
    params:
        key=config["ncbi_key"],
    resources:
        ncbi_requests=1,
    shell:
        """
        datasets summary genome taxon {wildcards.entry} --api-key {params.key} > {output}
        sleep 0.3
        """


rule json_to_tsv:
    input:
        os.path.join(outdir, "json", "{entry}.json"),
    output:
        temp(os.path.join(outdir, "json", "{entry}.tsv")),
    run:
        with open(input[0]) as file:
            data = json.load(file)
        df = pd.json_normalize(data, record_path=["reports"])
        df = df.replace(np.nan, "na")
        df["entry"] = [wildcards.entry] * len(df)
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


rule collect_summaries:
    input:
        expand(os.path.join(outdir, "json", "{entry}.tsv"), entry=entries),
    output:
        summary=temp(os.path.join(outdir, "json", "assembly_summary.tsv")),
        accessions=temp(os.path.join(outdir, "accessions.txt")),
    run:
        df = pd.concat([pd.read_csv(tsv, sep="\t") for tsv in input])
        df = df[~df.accession.str.contains(" ")]  # filter empty accessions
        df.to_csv(output.summary, sep="\t", index=None)
        df[["accession"]].to_csv(output.accessions, sep="\t", index=None, header=False)


rule download_archive:
    input:
        os.path.join(outdir, "accessions.txt"),
    output:
        os.path.join(outdir, "archive.zip"),
    params:
        key=config["ncbi_key"],
    shell:
        """
        datasets download genome accession \\
        --inputfile {input} \\
        --include genome,seq-report \\
        --api-key {key} \\
        --dehydrated --filename {output} 
        """


rule unzip_archive:
    input:
        os.path.join(outdir, "archive.zip"),
    output:
        temp(directory(os.path.join(outdir, "archive"))),
    shell:
        """
        unzip {input} -d {output}
        """


rule rehydrate_archive:
    input:
        os.path.join(outdir, "archive"),
    output:
        temp(os.path.join(outdir, "rehydrate.flag")),
    params:
        key=config["ncbi_key"],
    shell:
        """
        datasets rehydrate --directory {input} --api-key {key}
        touch {output}
        """


rule copy_fasta:
    input:
        dir=os.path.join(outdir, "archive"),
        flag=os.path.join(outdir, "rehydrate.flag"),
    output:
        directory(os.path.join(outdir, "assemblies")),
    shell:
        """
        mkdir {output}
        cp {input.dir}/ncbi_dataset/data/*/*.fna* {output}
        """


rule cat_sequence_reports:
    input:
        os.path.join(outdir, "archive"),
        os.path.join(outdir, "rehydrate.flag"),
    output:
        temp(os.path.join(outdir, "sequence_report.jsonl")),
    params:
        os.path.join(outdir, "archive", "ncbi_dataset", "data"),
    shell:
        """
        cat {params}/*/sequence_report.jsonl > {output}
        """


rule format_sequence_reports:
    input:
        os.path.join(outdir, "sequence_report.jsonl"),
    output:
        temp(os.path.join(outdir, "sequence_report.txt")),
    shell:
        """
        dataformat tsv genome-seq --inputfile {input} > {output}
        """


rule rename_seq_report_columns:
    input:
        os.path.join(outdir, "sequence_report.txt"),
    output:
        os.path.join(outdir, "sequence_report.tsv"),
    run:
        df = pd.read_csv(input[0], sep="\t")
        df.columns = df.columns.str.lower()
        df.columns = [col.replace(" ", "_") for col in df.columns]
        df.to_csv(output[0], sep="\t", index=None)


rule add_assembly_paths:
    input:
        dir=os.path.join(outdir, "assemblies"),
        summary=os.path.join(outdir, "json", "assembly_summary.tsv"),
    output:
        os.path.join(outdir, "assembly_summary.tsv"),
    run:
        df = pd.read_csv(input.summary, sep="\t")
        df = df.rename(
            columns={
                "organism.tax_id": "taxid",
                "assembly_stats.total_sequence_length": "genome_size",
            }
        )
        df["path"] = [
            os.path.abspath(glob.glob(os.path.join(f"{input.dir}", f"{acc}*.fna*"))[0])
            for acc in df["accession"]
        ]
        df.to_csv(output[0], sep="\t", index=None)


rule get_acc2taxid:
    input:
        seq=os.path.join(outdir, "sequence_report.tsv"),
        asm=os.path.join(outdir, "assembly_summary.tsv"),
    output:
        temp(os.path.join(outdir, "acc2taxid.txt")),
    params:
        accession,
    shell:
        """
        csvtk -t join -f 1 {input.seq} {input.asm} | \
        csvtk -t cut -f {params},taxid | \
        csvtk -t rename -f 1-2 -n contig,taxid > {output}
        """


rule get_lineage:
    input:
        os.path.join(outdir, "acc2taxid.txt"),
    output:
        temp(os.path.join(outdir, "lineage.txt")),
    shell:
        """
        csvtk -t cut -f taxid {input} | \
        csvtk uniq | \
        csvtk del-header | \
        taxonkit lineage | \
        taxonkit reformat | \
        csvtk -H -t cut -f 1,3 | \
        csvtk -H -t sep -f 2 -s ';' -R | \
        csvtk add-header -t \
        -n taxid,kindom,phylum,class,order,family,genus,species > {output}
        """


rule get_tax_report:
    input:
        acc2tax=os.path.join(outdir, "acc2taxid.txt"),
        lineage=os.path.join(outdir, "lineage.txt"),
    output:
        os.path.join(outdir, "sequence_taxonomy.tsv"),
    shell:
        """
        csvtk -t join -f taxid {input.acc2tax} {input.lineage} > {output}
        """
