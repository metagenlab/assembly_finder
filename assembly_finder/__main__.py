#  Parts of this wrapper are inspired by atlas, an easy-to-use metagenomic pipeline based on snakemake.
#  Go check it out on https://github.com/metagenome-atlas/atlas
import logging
import os
import sys
import rich_click as click
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

click.rich_click.OPTION_GROUPS = {
    "assembly_finder": [
        {
            "name": "Options",
            "options": [
                "--input",
                "--outdir",
                "--number",
                "--threads",
                "--taxon",
                "--rank",
                "--nrank",
                "--dryrun",
            ],
        },
        {
            "name": "NCBI datasets options",
            "options": [
                "--api-key",
                "--compressed",
                "--source",
                "--include",
                "--reference",
                "--assembly-level",
                "--annotated",
                "--atypical",
                "--mag",
            ],
        },
        {
            "name": "Help",
            "options": [
                "--help",
                "--version",
            ],
        },
    ]
}

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
    "--number",
    help="number of assemblies per query",
    type=int,
    default=None,
    show_default=True,
)
@click.option(
    "--dryrun",
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
@click.option("--api-key", type=str, help="NCBI api-key", default=None)
@click.option(
    "--compressed",
    type=bool,
    help="Download compressed files",
    default=True,
    show_default=True,
)
@click.option(
    "--include",
    type=str,
    help="Comma seperated files to download : genome,rna,protein,cds,gff3,gtf,gbff,seq-report,none",
    default="genome,seq-report",
    show_default=True,
)
@click.option(
    "--source",
    type=click.Choice(["refseq", "genbank", "all"], case_sensitive=False),
    help="Download from refseq or genbank",
    default="all",
    show_default=True,
)
@click.option(
    "--taxon",
    help="Are inputs taxa names or ids",
    type=bool,
    default=True,
    show_default=True,
)
@click.option(
    "--reference",
    type=bool,
    help="Limit to reference and representative genomes",
    default=False,
    show_default=True,
)
@click.option(
    "--assembly-level",
    help="Comma seperated list of assembly level: complete,chromosome,scaffold,contig",
    default=None,
    show_default=True,
)
@click.option(
    "--annotated",
    type=bool,
    help="Select annotated genomes only",
    default=False,
    show_default=True,
)
@click.option(
    "--atypical",
    type=bool,
    help="Exclude atypical genomes",
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
    "--rank",
    help="Select genomes at taxonomic rank",
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
    "--nrank",
    help="Number of genomes per taxonomic rank",
    type=int,
    default=None,
    show_default=True,
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
@click.version_option(version, "-v", "--version")
def cli(
    input,
    outdir,
    include,
    dryrun,
    threads,
    api_key,
    compressed,
    source,
    taxon,
    assembly_level,
    annotated,
    atypical,
    mag,
    reference,
    rank,
    nrank,
    number,
    snakemake_args,
):
    """
    \b
     ░█▀█░█▀▀░█▀▀░█▀▀░█▄█░█▀▄░█░░░█░█░░░█▀▀░▀█▀░█▀█░█▀▄░█▀▀░█▀▄
     ░█▀█░▀▀█░▀▀█░█▀▀░█░█░█▀▄░█░░░░█░░░░█▀▀░░█░░█░█░█░█░█▀▀░█▀▄
     ░▀░▀░▀▀▀░▀▀▀░▀▀▀░▀░▀░▀▀░░▀▀▀░░▀░░░░▀░░░▀▀▀░▀░▀░▀▀░░▀▀▀░▀░▀
    \b
    Snakemake-cli to download genomes with NCBI datasets
    """
    if outdir:
        outdir = os.path.abspath(outdir)
    else:
        outdir = os.path.abspath(os.path.basename(input).split(".")[0])

    if dryrun:
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
        f"--config api_key={api_key} "
        f"compressed={compressed} "
        f"input={input} "
        f"nb={number} "
        f"include={include} "
        f"source={source} "
        f"taxon={taxon} "
        f"assembly_level={assembly_level} "
        f"atypical={atypical} "
        f"mag={mag} "
        f"reference={reference} "
        f"annotated={annotated} "
        f"rank={rank} "
        f"nrank={nrank} "
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
