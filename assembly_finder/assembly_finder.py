#  Parts of this wrapper are inspired by atlas, an easy-to-use metagenomic pipeline based on snakemake.
#  Go check it out on https://github.com/metagenome-atlas/atlas
import logging
import os
import sys
import click
import subprocess


logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M",
    format="[%(asctime)s %(levelname)s] %(message)s",
)

"""
Functions
"""


def get_version():
    thisdir = os.path.abspath(os.path.dirname(__file__))
    init = os.path.join(thisdir, "__init__.py")
    return open(init).readline().split(" = ")[1].replace('"', "")


def get_snakefile():
    thisdir = os.path.abspath(os.path.dirname(__file__))
    sf = os.path.join(thisdir, "Snakefile")
    if not os.path.exists(sf):
        sys.exit(f"Unable to locate the Snakemake workflow file at {sf}")
    return sf


"""
Main
"""

version = get_version()

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
    "ignore_unknown_options": True,
}


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "-i",
    "--input",
    type=str,
    help="path to assembly_finder input table or list of queries",
    required=True,
)
@click.option("-o", "--outdir", help="output directory", type=click.Path())
@click.option(
    "-nb",
    "--n_query",
    help="number of assemblies per query",
    type=int,
    default=None,
    show_default=True,
)
@click.option(
    "-f",
    "--files",
    type=str,
    help="data files to include",
    default="genome,seq-report",
    show_default=True,
)
@click.option(
    "-n",
    "--dryrun_status",
    is_flag=True,
    default=False,
    show_default=True,
    help="snakemake dryrun to see the scheduling plan",
)
@click.option(
    "-t",
    "--threads",
    type=int,
    help="number of threads to allow for the workflow",
    default=2,
    show_default=True,
)
@click.option("-nk", "--ncbi_key", type=str, help="ncbi key for Entrez", default=None)
@click.option(
    "-db",
    "--database",
    type=click.Choice(["refseq", "genbank", "all"], case_sensitive=False),
    help="download from refseq or genbank",
    default="refseq",
    show_default=True,
)
@click.option(
    "--taxon",
    help="are inputs taxon names or ids",
    type=bool,
    default=True,
    show_default=True,
)
@click.option(
    "--reference",
    type=bool,
    help="limit to reference and representative genomes",
    default=False,
    show_default=True,
)
@click.option(
    "-al",
    "--assembly_level",
    help="select complete, chromosome, scaffold, contig",
    default=None,
    show_default=True,
)
@click.option(
    "-an",
    "--annotation",
    type=bool,
    help="select assemblies with annotation",
    default=False,
    show_default=True,
)
@click.option(
    "--atypical",
    type=bool,
    help="exclude atypical genomes",
    default=True,
    show_default=True,
)
@click.option(
    "--mag",
    type=click.Choice(["only", "exclude", "all"], case_sensitive=False),
    help="exclude or add MAGs to the dwnloaded genomes ",
    default="all",
    show_default=True,
)
@click.option(
    "-r",
    "--rank",
    help="taxonomic rank to filter by assemblies ",
    default=None,
    type=click.Choice(
        [
            "superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
        case_sensitive=False,
    ),
    show_default=True,
)
@click.option(
    "-nr",
    "--n_rank",
    help="number of genomes per taxonomic rank",
    type=int,
    default=None,
    show_default=True,
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
@click.version_option(version, "-v", "--version")
def cli(
    input,
    outdir,
    files,
    dryrun_status,
    threads,
    ncbi_key,
    database,
    taxon,
    assembly_level,
    annotation,
    atypical,
    mag,
    reference,
    rank,
    n_rank,
    n_query,
    snakemake_args,
):
    """
    \b
     ░█▀█░█▀▀░█▀▀░█▀▀░█▄█░█▀▄░█░░░█░█░░░█▀▀░▀█▀░█▀█░█▀▄░█▀▀░█▀▄
     ░█▀█░▀▀█░▀▀█░█▀▀░█░█░█▀▄░█░░░░█░░░░█▀▀░░█░░█░█░█░█░█▀▀░█▀▄
     ░▀░▀░▀▀▀░▀▀▀░▀▀▀░▀░▀░▀▀░░▀▀▀░░▀░░░░▀░░░▀▀▀░▀░▀░▀▀░░▀▀▀░▀░▀
    \b
    Snakemake pipeline to download genome assemblies from NCBI refseq/genbank
    """
    if outdir:
        outdir = os.path.abspath(outdir)
    else:
        outdir = os.path.abspath(os.path.basename(input).split(".")[0])

    if dryrun_status:
        dryrun = "-n"
    else:
        dryrun = ""

    if snakemake_args:
        args = " ".join([arg for arg in snakemake_args])
    else:
        args = ""

    cmd = (
        f"snakemake --snakefile {get_snakefile()} "
        f" --cores {threads} "
        f"all_download {dryrun} "
        f"--config ncbi_key={ncbi_key} "
        f"input={input} "
        f"nb={n_query} "
        f"files={files} "
        f"db={database} "
        f"taxon={taxon} "
        f"alvl={assembly_level} "
        f"atypical={atypical} "
        f"mag={mag} "
        f"reference={reference} "
        f"annot={annotation} "
        f"rank={rank} "
        f"nrank={n_rank} "
        f"outdir={outdir} "
        f"{args}"
    )
    logging.info(f"Executing: {cmd}")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)


if __name__ == "__main__":
    cli()
