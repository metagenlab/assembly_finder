rule download_taxdump:
    output:
        temp(os.path.join(TAXONKIT, "taxdump.tar.gz")),
    log:
        os.path.join(dir.out.logs, "curl.log"),
    conda:
        os.path.join(dir.env, "utils.yml")
    params:
        link=config.links.taxdump,
    shell:
        """
        curl {params.link} -o {output} 2> {log}
        """


rule decompress_taxdump:
    input:
        os.path.join(TAXONKIT, "taxdump.tar.gz"),
    output:
        os.path.join(TAXONKIT, "names.dmp"),
        os.path.join(TAXONKIT, "nodes.dmp"),
        os.path.join(TAXONKIT, "delnodes.dmp"),
        os.path.join(TAXONKIT, "merged.dmp"),
    params:
        dir=TAXONKIT,
    log:
        os.path.join(dir.out.logs, "tar.log"),
    conda:
        os.path.join(dir.env, "utils.yml")
    shell:
        """
        tar -xzvf {input} -C {params.dir} &> {log}
        """


if TAXON:

    rule taxon_genome_summary:
        output:
            temp(os.path.join(dir.out.json, "{query}.json")),
        log:
            os.path.join(dir.out.logs, "taxons", "{query}.log"),
        params:
            query=lambda wildcards: convert_query(wildcards),
            limit=lambda wildcards: get_limit(wildcards, LIMIT, QUERY2NB),
            args=ARGS,
            key=KEY,
        resources:
            ncbi_requests=1,
        retries: 2
        conda:
            os.path.join(dir.env, "datasets.yml")
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
            > {output} 2> {log}
            """

    rule collect_taxa_summaries:
        input:
            expand(os.path.join(dir.out.json, "{query}.json"), query=QUERIES),
        output:
            temp(os.path.join(dir.out.base, "genome_summaries.json")),
        run:
            dfs = []
            for file in input:
                query = str(os.path.basename(file).split(".json")[0])
                df = read_json(file)
                df.insert(0, "taxon", [query] * len(df))
                dfs.append(df)
            pd.concat(dfs).reset_index(drop=True).to_json(output[0])


else:
    if not os.path.isfile(INPUT):

        rule get_accessions_file:
            output:
                temp(os.path.join(dir.out.base, "queries.txt")),
            params:
                queries=INPUT.split(","),
            run:
                pd.DataFrame.from_dict({"queries": params.queries}).to_csv(
                    output[0], sep="\t", index=None, header=False
                )

        INPUT = os.path.join(dir.out.base, "queries.txt")

    rule accessions_genome_summary:
        input:
            INPUT,
        output:
            temp(os.path.join(dir.out.base, "genome_summaries.json")),
        log:
            os.path.join(dir.out.logs, "accessions_summary.log"),
        params:
            args=ARGS,
            key=KEY,
        resources:
            ncbi_requests=1,
        conda:
            os.path.join(dir.env, "datasets.yml")
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
        os.path.join(dir.out.base, "genome_summaries.json"),
    output:
        temp(os.path.join(dir.out.base, "taxids.txt")),
    run:
        read_json(input[0])["organism.tax_id"].drop_duplicates().to_csv(
            output[0], sep="\t", index=None, header=False
        )


rule get_lineage:
    input:
        taxids=os.path.join(dir.out.base, "taxids.txt"),
        names=os.path.join(TAXONKIT, "names.dmp"),
    output:
        temp(os.path.join(dir.out.base, "lineage.tsv")),
    log:
        os.path.join(dir.out.logs, "lineage.log"),
    params:
        headers=config.headers.lineage,
        dir=TAXONKIT,
    resources:
        ncbi_requests=1,
    conda:
        os.path.join(dir.env, "taxonkit.yml")
    shell:
        """
        taxonkit --data-dir {params.dir} lineage -r -n {input.taxids} | \\
        taxonkit --data-dir {params.dir} reformat | \\
        csvtk -H -t cut -f 1,4,3,5 | \\
        csvtk -H -t sep -f 4 -s ';' -R | \\
        csvtk add-header -t -n {params.headers} > {output} 2> {log}
        """


rule filter_genome_summaries:
    input:
        summary=os.path.join(dir.out.base, "genome_summaries.json"),
        lineage=os.path.join(dir.out.base, "lineage.tsv"),
    output:
        gen=temp(os.path.join(dir.out.base, "assembly_summary.txt")),
        tax=os.path.join(dir.out.base, "taxonomy.tsv"),
        acc=temp(os.path.join(dir.out.base, "accessions.txt")),
    params:
        rank=RANK,
        nrank=NRANK,
    script:
        os.path.join(dir.scripts, "select_assemblies.py")


rule archive_download:
    input:
        os.path.join(dir.out.base, "accessions.txt"),
    output:
        os.path.join(dir.out.base, "archive.zip"),
    log:
        os.path.join(dir.out.logs, "archive.log"),
    params:
        key=KEY,
        include=INCLUDE,
    conda:
        os.path.join(dir.env, "datasets.yml")
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
        os.path.join(dir.out.base, "archive.zip"),
    output:
        temp(directory(os.path.join(dir.out.base, "archive"))),
    log:
        os.path.join(dir.out.logs, "unzip.log"),
    conda:
        os.path.join(dir.env, "utils.yml")
    shell:
        """
        unzip {input} -d {output} &> {log}
        """


rule rehydrate_archive:
    input:
        os.path.join(dir.out.base, "archive"),
    output:
        temp(os.path.join(dir.out.base, "rehydrate.flag")),
    params:
        key=KEY,
        gzip=GZIP,
    conda:
        os.path.join(dir.env, "datasets.yml")
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
        dir=os.path.join(dir.out.base, "archive"),
        flag=os.path.join(dir.out.base, "rehydrate.flag"),
    output:
        directory(os.path.join(dir.out.download)),
    log:
        os.path.join(dir.out.logs, "rsync.log"),
    params:
        dir=os.path.join(dir.out.base, "archive", "ncbi_dataset", "data", "*"),
    conda:
        os.path.join(dir.env, "utils.yml")
    shell:
        """
        rsync -r {params.dir} {output} 2> {log}
        """


rule cat_sequence_reports:
    input:
        os.path.join(dir.out.download),
    output:
        os.path.join(dir.out.base, "sequence_report.tsv"),
    params:
        jsonl=os.path.join(dir.out.base, "download", "*", "sequence_report.jsonl"),
    conda:
        os.path.join(dir.env, "datasets.yml")
    shell:
        """
        cat {params.jsonl} | dataformat tsv genome-seq > {output}
        """


rule add_genome_paths:
    input:
        dir=os.path.join(dir.out.download),
        summary=os.path.join(dir.out.base, "assembly_summary.txt"),
    output:
        os.path.join(dir.out.base, "assembly_summary.tsv"),
    run:
        df = pd.read_csv(input.summary, sep="\t")
        df["path"] = [
            os.path.abspath(
                glob.glob(
                    os.path.join(f"{input.dir}", "*", f"{acc}*.fna*").replace(" ", "")
                )[0]
            )
            for acc in df["accession"]
        ]
        df.to_csv(output[0], sep="\t", index=None)


rule cleanup_reports:
    input:
        os.path.join(dir.out.base, "assembly_summary.tsv"),
        os.path.join(dir.out.base, "sequence_report.tsv"),
        os.path.join(dir.out.base, "taxonomy.tsv"),
    output:
        temp(os.path.join(dir.out.base, "cleanup.flag")),
    params:
        dir.out.base,
    conda:
        os.path.join(dir.env, "utils.yml")
    shell:
        """
        find {params[0]} -name "*.json*" -print0 | xargs -0 rm
        touch {output}
        """
