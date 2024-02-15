import pandas as pd
import numpy as np
import itertools
import os
import re
import glob
from itertools import chain
import sys
from io import StringIO
from ete3 import NCBITaxa
from functools import reduce

# Download params
download = config["download"]
# Path params
outdir = config["outdir"]
ete_db = config["ete_db"]
asmdir = os.path.join(outdir, "assemblies")

# Assemblies params
inp = config["input"]
nb = config["nb"]
db = config["db"]
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
# entry_df.to_csv("entry.tsv", sep="\t")
# Set entry as index
entry_df = entry_df.astype({"entry": str}).set_index("entry")

# Get entry list
entries = list(entry_df.index)


rule download_taxdump:
    output:
        temp(os.path.join(ete_db, "taxdump.tar.gz")),
    log:
        os.path.join(outdir, "logs", "ete.log"),
    retries: 2
    shell:
        """
        curl -o {output} https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz &> {log} 
        """


rule generate_ete3_NCBItaxa:
    input:
        os.path.join(ete_db, "taxdump.tar.gz"),
    output:
        os.path.join(ete_db, "taxa.sqlite"),
    log:
        os.path.join(outdir, "logs", "ete.log"),
    run:
        with open(log[0], "w") as f:
            sys.stderr = sys.stdout = f
            NCBITaxa(dbfile=output[0], taxdump_file=input[0])


rule get_assembly_tables:
    input:
        os.path.join(ete_db, "taxa.sqlite"),
    output:
        all=os.path.join(outdir, "tables", "{entry}-all.tsv"),
        filtered=os.path.join(outdir, "tables", "{entry}-filtered.tsv"),
    params:
        ncbi_key=config["ncbi_key"],
        ncbi_email=config["ncbi_email"],
        alvl=lambda wildcards: entry_df.loc[str(wildcards.entry)]["alvl"],
        db=lambda wildcards: entry_df.loc[str(wildcards.entry)]["db"],
        uid=lambda wildcards: entry_df.loc[str(wildcards.entry)]["uid"],
        rcat=lambda wildcards: entry_df.loc[str(wildcards.entry)]["rcat"],
        excl=lambda wildcards: entry_df.loc[str(wildcards.entry)]["excl"],
        annot=lambda wildcards: entry_df.loc[str(wildcards.entry)]["annot"],
        rank=lambda wildcards: entry_df.loc[str(wildcards.entry)]["rank"],
        n_by_rank=lambda wildcards: entry_df.loc[str(wildcards.entry)]["nrank"],
        nb=lambda wildcards: entry_df.loc[str(wildcards.entry)]["nb"],
    resources:
        ncbi_requests=1,
    log:
        os.path.join(outdir, "logs", "{entry}.log"),
    script:
        "assembly_table.py"


rule combine_assembly_tables:
    input:
        expand(os.path.join(outdir, "tables", "{entry}-filtered.tsv"), entry=entries),
    output:
        temp(os.path.join(outdir, "assemblies", "assembly_summary.tsv")),
    run:
        pd.concat([pd.read_csv(tb, sep="\t") for tb in list(input)]).to_csv(
            output[0], sep="\t", index=None
        )


rule get_ftp_links_list:
    input:
        os.path.join(outdir, "assemblies", "assembly_summary.tsv"),
    output:
        os.path.join(outdir, "links", "ftp-links.txt"),
    params:
        download,
    run:
        ftplinks = list(pd.read_csv(input[0], sep="\t")["ftp_path"])
        if params[0] == "aspera":
            ftplinks = [
                link.replace("ftp://ftp.ncbi.nlm.nih.gov", "") for link in ftplinks
            ]
        ftplinks = [os.path.join(link, os.path.basename(link)) for link in ftplinks]
        links = [link + "_" + sfx for sfx in sfxs.split(",") for link in ftplinks]
        with open(output[0], "w") as ftp_links:
            ftp_links.write("\n".join(str(link) for link in links))


if download == "aspera":

    checkpoint download:
        input:
            os.path.join(outdir, "links", "ftp-links.txt"),
        output:
            temp(os.path.join(outdir, "assemblies", "checksums.txt")),
        log:
            os.path.join(outdir, "logs", "download.log"),
        params:
            asmdir=asmdir,
        shell:
            """
            ascp -T -k 1 -i ${{CONDA_PREFIX}}/etc/aspera/aspera_bypass_dsa.pem --mode=recv --user=anonftp \
            --host=ftp.ncbi.nlm.nih.gov --file-list={input} --file-checksum=sha256 --file-manifest=text \
            --file-manifest-path={params.asmdir} {params.asmdir} &> {log} 
            mv {params.asmdir}/*.manifest.txt {output}
            """

elif download == "ftp":

    checkpoint split_links:
        input:
            os.path.join(outdir, "links", "ftp-links.txt"),
        output:
            temp(os.path.join(outdir, "links", "split.txt")),
        params:
            os.path.join(outdir, "links"),
        run:
            with open(input[0], "r") as file:
                lines = file.readlines()
                for line in lines:
                    line = line.strip("\n")
                    file = os.path.basename(line)
                    with open(output[0], "w") as out:
                        with open(
                            os.path.join(params[0], f"{file}.link"), "w"
                        ) as link:
                            out.write(line)
                            link.write(line)

    def get_query(file):
        with open(file, "r") as file:
            return file.readlines()

    checkpoint download:
        input:
            txt=os.path.join(outdir, "links", "split.txt"),
            query=lambda wildcards: storage.ftp(
                get_query(os.path.join(outdir, "links", f"{wildcards.ftp}.link"))
            ),
        output:
            os.path.join(outdir, "assemblies", "{ftp}"),
        log:
            os.path.join(outdir, "logs", "download", "{ftp}.log"),
        shell:
            """
            cp {input.query} {output} &> {log}
            """


def downloads(wildcards):
    if config["download"] == "aspera":
        extensions = config["sfxs"].split(",")
        checkpoint_directory = os.path.dirname(
            checkpoints.download.get(**wildcards).output[0]
        )
        files = [glob.glob(f"{checkpoint_directory}/*_{e}") for e in extensions]
        return list(chain.from_iterable(files))

    elif config["download"] == "ftp":
        checkdir = os.path.dirname(checkpoints.split_links.get(**wildcards).output[0])
        return expand(
            os.path.join(outdir, "assemblies", "{ftp}"),
            ftp=glob_wildcards(os.path.join(checkdir, "{i}.link")).i,
        )


def get_reports(file):
    keys = ("# Assembly name", "# Assembly method", "# Sequencing technology")
    always_print = False
    seq_lines = ""
    asm_dic = {
        "Assembly name": np.nan,
        "Assembly method": np.nan,
        "Sequencing technology": np.nan,
    }
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
        dl=downloads,
    output:
        temp(os.path.join(outdir, "assemblies", "assembly_reports.tsv")),
        temp(os.path.join(outdir, "assemblies", "sequence_reports.tsv")),
    run:
        reports = [f for f in input.dl if "assembly_report" in f]
        all_reports = [get_reports(report) for report in reports]
        asm_rep = pd.concat([report[0] for report in all_reports]).reset_index(
            drop=True
        )
        asm_rep.columns = ["asm_name", "asm_method", "seq_tech"]
        seq_rep = pd.concat([report[1] for report in all_reports])
        asm_rep.to_csv(output[0], sep="\t", index=None)
        seq_rep.to_csv(output[1], sep="\t", index=None)


rule get_summaries:
    input:
        summary=os.path.join(outdir, "assemblies", "assembly_summary.tsv"),
        asm_report=os.path.join(outdir, "assemblies", "assembly_reports.tsv"),
        seq_report=os.path.join(outdir, "assemblies", "sequence_reports.tsv"),
    output:
        os.path.join(outdir, "assembly_summary.tsv"),
        os.path.join(outdir, "sequence_summary.tsv"),
        os.path.join(outdir, "taxonomy_summary.tsv"),
    params:
        asmdir=asmdir,
    run:
        df = pd.read_csv(input.summary, sep="\t")
        asm_report = pd.read_csv(input.asm_report, sep="\t")
        seq_report = pd.read_csv(input.seq_report, sep="\t")
        # replace ftp paths with absolute paths
        df["path"] = [
            os.path.abspath(
                glob.glob(
                    os.path.join(params.asmdir, f"{os.path.basename(ftp)}*.fna.gz")
                )[0]
            )
            for ftp in df["ftp_path"]
        ]
        dfs = [df, asm_report, seq_report]
        merge_df = reduce(lambda left, right: pd.merge(left, right, on="asm_name"), dfs)
        asm_df = merge_df[
            [
                "entry",
                "database",
                "db_uid",
                "asm_name",
                "organism",
                "taxid",
                "asm_release_date",
                "asm_status",
                "refseq_category",
                "contig_count",
                "contig_n50",
                "contig_l50",
                "genome_size",
                "coverage",
                "asm_method",
                "seq_tech",
                "path",
            ]
        ].drop_duplicates()
        tax_df = merge_df[
            [
                "asm_name",
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
            ]
        ].drop_duplicates()
        seq_df = merge_df[
            [
                "asm_name",
                "organism",
                "taxid",
                "genbank_accn",
                "refseq_accn",
                "assigned_molecule_location/type",
                "sequence_length",
            ]
        ]
        asm_df.to_csv(output[0], sep="\t", index=None)
        seq_df.to_csv(output[1], sep="\t", index=None)
        tax_df.to_csv(output[2], sep="\t", index=None)


rule clean_files:
    input:
        os.path.join(outdir, "assembly_summary.tsv"),
        os.path.join(outdir, "sequence_summary.tsv"),
        os.path.join(outdir, "taxonomy_summary.tsv"),
    output:
        temp(os.path.join(outdir, "clean.txt")),
    params:
        links=os.path.join(outdir, "links"),
    shell:
        """
        rm -r {params.links}
        touch {output}
        """


if download == "aspera":

    rule format_checksum:
        input:
            os.path.join(outdir, "assemblies", "checksums.txt"),
        output:
            temp(os.path.join(outdir, "assemblies", "aspera-checks.txt")),
        run:
            d = {
                line.replace('"', "")
                .replace("\n", "")
                .split(", ")[1]
                .split(":")[1]: os.path.relpath(
                    line.replace('"', "")
                    .replace("\n", "")
                    .split(", ")[0]
                    .split(" ")[0]
                )
                for line in open(input[0])
                if line.startswith('"')
            }
            pd.DataFrame.from_dict([d]).transpose().reset_index().to_csv(
                output[0], sep="\t", index=None, header=None
            )

    def get_ext(wildcards, asm_dir, exts):
        asm_dir = os.path.relpath(asm_dir)
        if len(exts.split(",")) > 1:
            return os.path.join(asm_dir, "*" + "{" + exts + "}")
        else:
            return os.path.join(asm_dir, "*" + exts)

    rule verify_checksums:
        input:
            os.path.join(outdir, "assemblies", "aspera-checks.txt"),
            lambda wildcards: downloads(wildcards, download, sfxs),
        output:
            temp(os.path.join(outdir, "assemblies", "sha256.txt")),
        params:
            lambda wildcards: get_ext(wildcards, asmdir, sfxs),
        shell:
            """
            sha256sum {params} | sed 's/ \+ /\t/g' > {output}
            diff <(sort {output}) <(sort {input[0]})
            """
