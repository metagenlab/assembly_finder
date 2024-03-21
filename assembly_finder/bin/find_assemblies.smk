import json
import pandas as pd
import glob
import os

outdir = config["outdir"]
accession = "refseq_seq_accession"
if config["source"] == "genbank":
    accession = "refseq_seq_accession"


args = ""
if config["annotated"]:
    args += "--annotated "
if config["api_key"] != "None":
    args += f"--api-key {config['api_key']} "
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
if config["nb"] != "None":
    args += f"--limit {config['nb']} "


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


try:
    df = pd.read_csv(os.path.abspath(config["input"]), sep="\t", header=None)
    if config["taxon"]:
        queries = list(df[0])
    else:
        queries = config["input"]

except FileNotFoundError:
    queries = config["input"].split(",")
    if config["taxon"] == False:
        pd.DataFrame.from_dict({"queries": queries}).to_csv(
            os.path.join("queries.txt"), sep="\t", header=None, index=None
        )
        queries = os.path.join("queries.txt")


if config["taxon"]:

    rule taxon_genome_summary:
        output:
            temp(os.path.join(outdir, "json", "{query}.json")),
        params:
            query=lambda wildcards: convert_query(wildcards),
            args=args,
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
              {params.args} \\
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
        config["api_key"],
    resources:
        ncbi_requests=1,
    shell:
        """
        datasets \\
          summary \\
          taxonomy \\
          taxon \\
          --inputfile {input} \\
          --api-key {params[0]} \\
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
        key=config["api_key"],
        include=config["include"],
    shell:
        """
        datasets \\
          download \\
          genome \\
          accession \\
          --inputfile {input} \\
          --include {params.include} \\
          --api-key {params.key} \\
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


gzip = ""
if config["compressed"]:
    gzip = "--gzip"


rule rehydrate_archive:
    input:
        os.path.join(outdir, "archive"),
    output:
        temp(os.path.join(outdir, "rehydrate.flag")),
    params:
        key=config["api_key"],
        gzip=gzip,
    shell:
        """
        datasets \\
          rehydrate \\
          --directory {input} \\
          --api-key {params.key} \\
          {params.gzip}
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
        dir=os.path.join(outdir, "download", "*", "sequence_report.jsonl"),
    shell:
        """
        cat {params.dir} | dataformat tsv genome-seq > {output}
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
