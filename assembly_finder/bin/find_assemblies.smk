import pandas as pd
import os
import re
import glob
import sys
from io import StringIO
from ete3 import NCBITaxa

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

if os.path.isfile(inp):
    entries = list(pd.read_csv(inp, sep="\t", header=None)[0])
    entry_to_nb = pd.read_csv(inp, sep="\t", names=["entry", "nb"], index_col="entry")

else:
    entries = inp.split(",")
    entry_to_nb = pd.DataFrame()


rule download_taxdump:
    output:
        temp("taxdump.tar.gz"),
    log:
        f"{outdir}/logs/ete3.log",
    shell:
        """
        wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz &> {log} 
        """


rule generate_ete3_NCBItaxa:
    input:
        "taxdump.tar.gz",
    output:
        "taxa.sqlite",
    log:
        f"{outdir}/logs/ete3.log",
    run:
        with open(snakemake.log[0], "w") as f:
            sys.stderr = sys.stdout = f
            NCBITaxa(dbfile=output[0], taxdump_file=input[0])


def get_nb(entry, df, default_nb):
    if df.empty or df.dropna().nb.empty:
        return default_nb
    else:
        entry = int(entry)
        return int(df.loc[entry].values)


rule get_assembly_tables:
    input:
        "taxa.sqlite",
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
        nb=lambda wildcards: get_nb(wildcards.entry, entry_to_nb, nb),
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
        temp(f"{outdir}/assemblies/assembly_summary.tsv"),
    run:
        pd.concat([pd.read_csv(tb, sep="\t") for tb in list(input)]).to_csv(
            output[0], sep="\t", index=None
        )


rule get_ftp_links_list:
    input:
        f"{outdir}/assemblies/assembly_summary.tsv",
    output:
        temp(f"{outdir}/assemblies/ftp-links.txt"),
    params:
        db=db,
    run:
        ftplinks = pd.read_csv(input[0], sep="\t")["ftp_path"]
        links = []
        with open(output[0], "w") as ftp_links:
            for link in ftplinks:
                src_dir = link.split("/")[-1]
                link = link.replace("ftp://ftp.ncbi.nlm.nih.gov", "")
                ftp_links.write(link + "/" + src_dir + "_genomic.fna.gz\n")
                ftp_links.write(link + "/" + src_dir + "_assembly_report.txt\n")


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
        asmdir=f"{outdir}/assemblies",
    shell:
        """
        ascp -T -k 1 -i ${{CONDA_PREFIX}}/etc/asperaweb_id_dsa.openssh --mode=recv --user=anonftp \
        --host=ftp.ncbi.nlm.nih.gov --file-list={input} --file-checksum=sha256 --file-manifest=text \
        --file-manifest-path={params.asmdir} {params.asmdir} &> {log} 
        mv {params.asmdir}/*.manifest.txt {output}
        """


def downloads(wildcards, extension):
    checkpoint_output = checkpoints.download_assemblies.get(**wildcards).output[0]
    directory = "/".join((checkpoint_output.split("/")[0:2]))
    return expand(
        f"{outdir}/assemblies/{{i}}{extension}",
        i=glob_wildcards(os.path.join(directory, f"{{i}}{extension}")).i,
    )


def get_reports(file):
    keys = ("# Assembly name", "# Assembly method", "# Sequencing technology")
    always_print = False
    seq_lines = ""
    asm_dic = {}
    for line in open(file):
        if line.startswith(keys):
            d = {
                re.split(r":\s{1,}", line)[0]
                .replace("# ", ""): re.split(r":\s{1,}", line)[1]
                .replace("\n", "")
            }
            asm_dic.update(d)

        elif always_print or line.startswith("# Sequence-Name"):
            seq_lines += line.replace("# ", "")
            always_print = True
    asm_df = pd.DataFrame.from_dict([asm_dic])
    seq_df = pd.read_csv(StringIO(seq_lines), sep="\t")
    seq_df.columns = seq_df.columns.str.lower().str.replace("-", "_")
    seq_df.insert(
        loc=0, value=[asm_dic["Assembly name"]] * len(seq_df), column="asm_name"
    )
    return (asm_df, seq_df)


rule get_assembly_reports:
    input:
        lambda wildcards: downloads(wildcards, "_assembly_report.txt"),
    output:
        temp(f"{outdir}/assemblies/assembly_reports.tsv"),
        temp(f"{outdir}/assemblies/sequence_reports.tsv"),
    run:
        all_reports = [get_reports(report) for report in input]
        asm_rep = pd.concat([report[0] for report in all_reports]).reset_index(
            drop=True
        )
        asm_rep.columns = ["asm_name", "asm_method", "seq_tech"]
        asm_rep.to_csv(output[0], sep="\t", index=None)
        pd.concat([report[1] for report in all_reports]).to_csv(
            output[1], sep="\t", index=None
        )


rule get_summaries:
    input:
        summary=f"{outdir}/assemblies/assembly_summary.tsv",
        asm_report=f"{outdir}/assemblies/assembly_reports.tsv",
        seq_report=f"{outdir}/assemblies/sequence_reports.tsv",
    output:
        f"{outdir}/assembly_summary.tsv",
        f"{outdir}/sequence_summary.tsv",
    params:
        asmdir=f"{outdir}/assemblies",
    run:
        df = pd.read_csv(input.summary, sep="\t")
        asm_report = pd.read_csv(input.asm_report, sep="\t")
        seq_report = pd.read_csv(input.seq_report, sep="\t")
        # replace ftp paths with absolute paths
        df.ftp_path = [
            os.path.abspath(glob.glob(f"{params.asmdir}/{acc}*.fna.gz")[0])
            for acc in df["asm_accession"]
        ]
        df.rename(columns={"ftp_path": "path"}, inplace=True)
        asm_df = df.merge(asm_report, on="asm_name")
        asm_df = asm_df[
            [
                "entry",
                "database",
                "db_uid",
                "asm_accession",
                "asm_name",
                "asm_status",
                "refseq_category",
                "contig_count",
                "contig_n50",
                "genome_size",
                "coverage",
                "asm_method",
                "seq_tech",
                "organism",
                "sub_type",
                "sub_value",
                "taxid",
                "superkingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
                "path",
            ]
        ]
        seq_df = seq_report.merge(df, on="asm_name")
        seq_df = seq_df[
            [
                "asm_name",
                "asm_accession",
                "genbank_accn",
                "refseq_accn",
                "sequence_name",
                "assigned_molecule_location/type",
                "sequence_length",
                "taxid",
            ]
        ]
        asm_df.to_csv(output[0], sep="\t", index=None)
        seq_df.to_csv(output[1], sep="\t", index=None)


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


rule verify_checksums:
    input:
        f"{outdir}/assemblies/aspera-checks.txt",
        lambda wildcards: downloads(wildcards, ".fna.gz"),
    output:
        f"{outdir}/assemblies/sha256.txt",
    params:
        asmdir=f"{outdir}/assemblies",
    shell:
        """
        sha256sum {params.asmdir}/*{{.fna.gz,assembly_report.txt}} | sed 's/ \+ /\t/g' > {output} 
        diff <(sort {output}) <(sort {input[0]})
        """
