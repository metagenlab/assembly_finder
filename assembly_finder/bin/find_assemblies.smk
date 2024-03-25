from pathlib import Path
import json
import pandas as pd
import numpy as np
import glob
import os

# read config
outdir = config["outdir"]
inp = str(config["input"])
nbs = str(config["nb"])
# get ncbi-datasets-cli args
key = ""
if config["api_key"] != "None":
    key += f"--api-key {config['api_key']} "
args = ""
if config["annotated"]:
    args += "--annotated "
if config["assembly_level"] != "None":
    args += f"--assembly-level {config['assembly_level']} "
if config["source"]:
    args += f"--assembly-source {config['source']} "
if config["atypical"]:
    args += "--exclude-atypical "
if config["mag"]:
    args += f"--mag {config['mag']} "
if config["reference"]:
    args += "--reference "
gzip = ""
if config["compressed"]:
    gzip = "--gzip"


# functions
def read_json(file):
    try:
        return pd.json_normalize(json.load(open(file)), record_path=["reports"])
    except KeyError:
        return pd.read_json(file)


def convert_query(wildcards):
    query = wildcards.query
    if ("_" in query) and (("GCF" not in query) and ("GCA" not in query)):
        query = f'"{query.replace("_", " ")}"'
    return query


def get_limit(wildcards, nbs, dic):
    if nbs != "None":
        return dic[wildcards.query]
    else:
        return ""


# get taxa queries and genome numbers
if config["taxon"]:
    try:
        df = pd.read_csv(os.path.abspath(inp), sep="\t")
        queries = list(df["taxon"])
        if "nb" in df.columns:
            nbs = list(df["nb"])
    except (FileNotFoundError,IsADirectoryError):
        queries = inp.split(",")

    queries = [str(query) for query in queries]
    if nbs != "None":
        if type(nbs) is not list:
            nbs = nbs.split(",")
        nbs = [f"--limit {nb}" for nb in nbs]
        if len(nbs) == 1:
            nbs = nbs * len(queries)
        query2nb = dict(zip(queries, nbs))
else:
    queries = inp


if config["taxon"]:

    rule taxon_genome_summary:
        output:
            temp(os.path.join(outdir, "json", "{query}.json")),
        params:
            query=lambda wildcards: convert_query(wildcards),
            limit=lambda wildcards: get_limit(wildcards, nbs, query2nb),
            args=args,
            key=key,
        resources:
            ncbi_requests=1,
        retries: 2
        shell:
            """
            datasets \\
              summary \\
              genome \\
              taxon \\
              {params.query} \\
              {params.limit} \\
              {params.args} \\
              {params.key} \\
              > {output}
            """

    rule collect_taxa_summaries:
        input:
            expand(os.path.join(outdir, "json", "{query}.json"), query=queries),
        output:
            temp(os.path.join(outdir, "genome_summaries.json")),
        run:
            pd.concat([read_json(file) for file in input]).reset_index(
                drop=True
            ).to_json(output[0])

else:

    rule accessions_genome_summary:
        input:
            queries,
        output:
            temp(os.path.join(outdir, "genome_summaries.json")),
        params:
            args=args,
            key=key,
        resources:
            ncbi_requests=1,
        shell:
            """
            datasets \\
              summary \\
              genome \\
              accession \\
              --inputfile {input} \\
              {params.args} \\
              {params.key} \\
              > {output}
            """


rule get_taxids:
    input:
        os.path.join(outdir, "genome_summaries.json"),
    output:
        temp(os.path.join(outdir, "taxids.txt")),
    params:
        config["taxon"],
    run:
        read_json(input[0])["organism.tax_id"].drop_duplicates().to_csv(
            output[0], sep="\t", index=None, header=False
        )


rule get_lineage:
    input:
        os.path.join(outdir, "taxids.txt"),
    output:
        temp(os.path.join(outdir, "lineage.json")),
    params:
        key=key,
    resources:
        ncbi_requests=1,
    shell:
        """
        datasets \\
          summary \\
          taxonomy \\
          taxon \\
          --inputfile {input} \\
          {params.key} \\
          > {output}
        """


rule filter_genome_summaries:
    input:
        summary=os.path.join(outdir, "genome_summaries.json"),
        lineage=os.path.join(outdir, "lineage.json"),
    output:
        gen=temp(os.path.join(outdir, "assembly_summary.txt")),
        tax=os.path.join(outdir, "taxonomy.tsv"),
        acc=temp(os.path.join(outdir, "accessions.txt")),
    params:
        rank=config["rank"],
        nrank=config["nrank"],
    script:
        "select_assemblies.py"


rule archive_download:
    input:
        os.path.join(outdir, "accessions.txt"),
    output:
        os.path.join(outdir, "archive.zip"),
    params:
        key=key,
        include=config["include"],
    shell:
        """
        datasets \\
          download \\
          genome \\
          accession \\
          --inputfile {input} \\
          --include {params.include} \\
          {params.key} --dehydrated \\
          --filename {output} 
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
        key=key,
        gzip=gzip,
    shell:
        """
        datasets \\
          rehydrate \\
          --directory {input} \\
          {params.key} {params.gzip}
        touch {output}
        """


rule copy_files:
    input:
        dir=os.path.join(outdir, "archive"),
        flag=os.path.join(outdir, "rehydrate.flag"),
    output:
        directory(os.path.join(outdir, "download")),
    params:
        dir=os.path.join(outdir, "archive", "ncbi_dataset", "data", "*"),
    shell:
        """
        rsync -r {params.dir} {output}
        """


rule cat_sequence_reports:
    input:
        os.path.join(outdir, "download"),
    output:
        os.path.join(outdir, "sequence_report.tsv"),
    params:
        jsonl=os.path.join(outdir, "download", "*", "sequence_report.jsonl"),
    shell:
        """
        cat {params.jsonl} | dataformat tsv genome-seq > {output}
        """


rule add_genome_paths:
    input:
        dir=os.path.join(outdir, "download"),
        summary=os.path.join(outdir, "assembly_summary.txt"),
    output:
        os.path.join(outdir, "assembly_summary.tsv"),
    run:
        df = pd.read_csv(input.summary, sep="\t")
        df["path"] = [
            os.path.abspath(
                glob.glob(os.path.join(f"{input.dir}", "*", f"{acc}*.fna*"))[0]
            )
            for acc in df["accession"]
        ]
        df.to_csv(output[0], sep="\t", index=None)
