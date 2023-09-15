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
        nb=lambda wildcards: get_nb(
            wildcards.entry, inp, config["n_by_entry"]
        )
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
        f"{outdir}/assemblies_summary.tsv",
    run:
        pd.concat([pd.read_csv(tb, sep="\t") for tb in list(input)]).to_csv(
            output[0], sep="\t", index=None
        )


rule get_ftp_links_list:
    input:
        f"{outdir}/assemblies_summary.tsv",
    output:
        temp(f"{outdir}/assemblies/ftp-links.txt"),
    params:
        db=db,
        dl=dl,
    run:
        if params.db == "refseq":
            ftplinks = pd.read_csv(input[0], sep="\t")["FtpPath_RefSeq"]
        else:
            ftplinks = pd.read_csv(input[0], sep="\t")["FtpPath_GenBank"]
        links = []
        for link in ftplinks:
            if params.dl == "aspera":
                link = link.replace("ftp://ftp.ncbi.nlm.nih.gov", "")
            link += "/" + link.split("/")[-1] + "_genomic.fna.gz\n"
            links.append(link)
            f = open(output[0], "w")
            f.writelines(links)
            f.close()


checkpoint download_assemblies:
    input:
        f"{outdir}/assemblies/ftp-links.txt",
    output:
        f"{outdir}/assemblies/dl.done",
    log:
        f"{outdir}/logs/download.log",
    benchmark:
        f"{outdir}/benchmark/downloads.txt"
    params:
        asmdir=f"{outdir}/assemblies",
        dl=dl,
    run:
        if dl == "aspera":
            shell(
                """
                ascp -T -k 1 -i ${{CONDA_PREFIX}}/etc/asperaweb_id_dsa.openssh --mode=recv --user=anonftp \
                --host=ftp.ncbi.nlm.nih.gov --file-list={input} {params.asmdir} &> {log} 
                touch {output}
                """
            )
        elif dl == "wget":
            shell(
                """
                wget -i {input} -P {params.asmdir} &> {log}
                touch {output}
                """
            )
