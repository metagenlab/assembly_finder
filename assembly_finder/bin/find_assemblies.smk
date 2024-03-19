import json
import pandas as pd
import glob
import os

outdir = config["outdir"]
accession = "refseq_seq_accession"
if config["db"] == "genbank":
    accession = "refseq_seq_accession"


def datasets_global_flags(wildcards):
    args = ""
    if config["annot"]:
        args += "--annotated "
    if config["ncbi_key"]:
        args += f"--api-key {config['ncbi_key']} "
    if config["alvl"] != "None":
        args += f"--assembly-level {config['alvl']} "
    if config["db"]:
        args += f"--assembly-source {config['db']} "
    if config["atypical"]:
        args += "--exclude-atypical "
    if config["mag"]:
        args += f"--mag {config['mag']} "
    if config["reference"]:
        args += "--reference "
    if config["nb"] != "None":
        args += f"--limit {config['nb']} "
    return args


def read_json(file):
    try:
        return pd.json_normalize(json.load(open(file)), record_path=["reports"])
    except KeyError:
        return pd.read_json(file)


checkpoint read_queries:
    output:
        temp(os.path.join(outdir, "queries.txt")),
    params:
        str(config["input"]),
    run:
        df = {"queries": []}
        if os.path.isfile(params[0]):
            queries = pd.read_csv(params[0], sep="\t", header=None)[0]
        else:
            queries = params[0].split(",")

        for query in queries:
            if ("_" in query) and (("GCF" not in query) and ("GCA" not in query)):
                query = query.replace("_", " ")
            df["queries"].append(query)
        df = pd.DataFrame.from_records(df)
        df.to_csv(output[0], sep="\t", index=None)


if config["taxon"]:

    rule taxon_genome_summary:
        output:
            temp(os.path.join(outdir, "taxon", "{query}.json")),
        params:
            args=datasets_global_flags,
        resources:
            ncbi_requests=1,
        retries: 2
        shell:
            """
            datasets \\
              summary \\
              genome \\
              taxon \\
              {wildcards.query} \\
              {params.args} \\
              > {output}
            sleep 0.3
            """

    def aggregate_summaries(wildcards):
        checkout = checkpoints.read_queries.get(**wildcards).output[0]
        queries = list(pd.read_csv(checkout, sep="\t")["queries"])
        return expand(os.path.join(outdir, "taxon", "{query}.json"), query=queries)

    rule collect_taxa_summaries:
        input:
            flag=os.path.join(outdir, "queries.txt"),
            files=aggregate_summaries,
        output:
            temp(os.path.join(outdir, "genome_summaries.json")),
        run:
            pd.concat([read_json(file) for file in input.files]).to_json(output[0])

else:

    rule accessions_genome_summary:
        input:
            os.path.join(outdir, "queries.txt"),
        output:
            temp(os.path.join(outdir, "genome_summaries.json")),
        params:
            args=datasets_global_flags,
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
        config["ncbi_key"],
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
        nb=config["nb"],
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
        key=config["ncbi_key"],
        files=config["files"],
    shell:
        """
        datasets \\
          download \\
          genome \\
          accession \\
          --inputfile {input} \\
          --include {params.files} \\
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


rule rehydrate_archive:
    input:
        os.path.join(outdir, "archive"),
    output:
        temp(os.path.join(outdir, "rehydrate.flag")),
    params:
        config["ncbi_key"],
    shell:
        """
        datasets \\
          rehydrate \\
          --directory {input} \\
          --api-key {params[0]}
        touch {output}
        """


rule copy_fasta:
    input:
        dir=os.path.join(outdir, "archive"),
        flag=os.path.join(outdir, "rehydrate.flag"),
    output:
        directory(os.path.join(outdir, "download")),
    shell:
        """
        mkdir {output}
        cp -r {input.dir}/ncbi_dataset/data/* {output}
        """


rule cat_sequence_reports:
    input:
        os.path.join(outdir, "download"),
    output:
        os.path.join(outdir, "sequence_report.tsv"),
    shell:
        """
        cat {input}/*/sequence_report.jsonl |\\
        dataformat tsv genome-seq > {output}
        """


rule add_assembly_paths:
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
