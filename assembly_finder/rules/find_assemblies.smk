import pandas as pd
import os

inp = str(config["input"])
outdir = config["outdir"]
ncbi_key = config["NCBI_key"]
ncbi_email = config["NCBI_email"]
uid = config["uid"]
alvl = config["assembly_level"]
db = config["db"]
rcat = config["refseq_category"]
excl = config["exclude"]
annot = config["annotation"]
rank = config["Rank_to_filter_by"]
nrank = config["n_by_rank"]
nb = config["n_by_entry"]
dl = config["downloader"]

if os.path.isfile(inp):
    entries = list(pd.read_csv(inp, sep="\t", header=None)[0])
else:
    entries = inp.split(",")


rule check_for_update_ete3:
    output:
        temp("ete3-update.txt"),
    log:
        f"{outdir}/logs/ete3-update.log",
    script:
        "update-ete3.py"


def get_nb(entry, tb, nb):
    if os.path.isfile(tb):
        df = pd.read_csv(tb, sep="\t", names=["entry", "nb"], index_col="entry")
        return int(df.loc[int(entry)]["nb"])
    else:
        return nb


rule get_assembly_tables:
    input:
        "ete3-update.txt",
    output:
        all=f"{outdir}/tables/{{entry}}-all.tsv",
        filtered=f"{outdir}/tables/{{entry}}-filtered.tsv",
    params:
        ncbi_key=ncbi_key,
        ncbi_email=ncbi_email,
        alvl=alvl,
        db=db,
        uid=uid,
        rcat=rcat,
        excl=excl,
        annot=annot,
        rank=rank,
        n_by_rank=nrank,
        nb=lambda wildcards: get_nb(wildcards.entry, inp, config["n_by_entry"]),
    resources:
        ncbi_requests=1,
    log:
        f"{outdir}/logs/find-assemblies-{{entry}}.log",
    benchmark:
        f"{outdir}/benchmark/find-assemblies-{{entry}}.txt"
    script:
        "assembly_table.py"


rule combine_assembly_tables:
    input:
        expand(f"{outdir}/tables/{{entry}}-filtered.tsv", entry=entries),
    output:
        f"{outdir}/summary.tsv",
    run:
        pd.concat([pd.read_csv(tb, sep="\t") for tb in list(input)]).to_csv(
            output[0], sep="\t", index=None
        )


rule get_ftp_links_list:
    input:
        f"{outdir}/summary.tsv",
    output:
        temp(f"{outdir}/assemblies/ftp-links.txt"),
    params:
        db=db,
        dl=dl,
    run:
        ftplinks = pd.read_csv(input[0], sep="\t")["ftp_path"]
        links = []
        with open(output[0], "w") as ftp_links:
            for link in ftplinks:
                src_dir = link.split("/")[-1]
                if params.dl == "aspera":
                    link = link.replace("ftp://ftp.ncbi.nlm.nih.gov", "")
                ftp_links.write(link + "/" + src_dir + "_genomic.fna.gz\n")


checkpoint download_assemblies:
    input:
        f"{outdir}/assemblies/ftp-links.txt",
    output:
        temp(f"{outdir}/assemblies/checksums.txt"),
    log:
        f"{outdir}/logs/download.log",
    benchmark:
        f"{outdir}/benchmark/downloads.txt"
    params:
        dl=dl,
        asmdir=f"{outdir}/assemblies",
    run:
        if dl == "aspera":
            shell(
                """
                ascp -T -k 1 -i ${{CONDA_PREFIX}}/etc/asperaweb_id_dsa.openssh --mode=recv --user=anonftp \
                --host=ftp.ncbi.nlm.nih.gov --file-list={input} --file-checksum=sha256 --file-manifest=text \
                --file-manifest-path={params.asmdir} {params.asmdir} &> {log} 
                mv {params.asmdir}/*.manifest.txt {output}
                """
            )
        elif dl == "wget":
            shell(
                """
                wget -i {input} -P {params.asmdir} &> {log}
                """
            )


rule format_checksum:
    input:
        f"{outdir}/assemblies/checksums.txt",
    output:
        temp(f"{outdir}/assemblies/aspera-checks.txt"),
    run:
        d = {
            line.replace('"', "")
            .replace("\n", "")
            .split(", ")[1]
            .split(":")[1]: os.path.relpath(
                line.replace('"', "").replace("\n", "").split(", ")[0].split(" ")[0]
            )
            for line in open(input[0])
            if line.startswith('"')
        }
        pd.DataFrame.from_dict([d]).transpose().reset_index().to_csv(
            output[0], sep="\t", index=None, header=None
        )


def downloads(wildcards):
    checkpoint_output = checkpoints.download_assemblies.get(**wildcards).output[0]
    directory = "/".join((checkpoint_output.split("/")[0:2]))
    return expand(
        f"{outdir}/assemblies/{{i}}.fna.gz",
        i=glob_wildcards(os.path.join(directory, "{i}.fna.gz")).i,
    )


rule verify_checksums:
    input:
        f"{outdir}/assemblies/aspera-checks.txt",
        downloads,
    output:
        temp(f"{outdir}/assemblies/sha256.txt"),
    params:
        asmdir=f"{outdir}/assemblies",
    shell:
        """
        sha256sum {params.asmdir}/*.fna.gz | sed 's/ \+ /\t/g' > {output} 
        diff <(sort {output}) <(sort {input[0]})
        """
