"""
Entrypoint for Assembly finder

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import rich_click as click

from snaketool_utils.cli_utils import run_snakemake, echo_click


def snake_base(rel_path):
    """Get the filepath to a Snaketool system file (relative to __main__.py)"""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    """Read and print the version from the version file"""
    with open(snake_base("assembly_finder.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    """Read and print the Citation information from the citation file"""
    with open(snake_base("assembly_finder.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_outdir(ctx, param, value):
    """Callback for click options; gets default output dir as input basename"""
    if not value:
        return os.path.splitext(ctx.params["input"])[0]
    return value


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "-o",
            "--output",
            type=click.Path(),
            default=None,
            callback=default_outdir,
            help="Output directory",
        ),
        click.option(
            "--print-versions/--no-print-versions",
            default=False,
            help="Print all tool versions at workflow end",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/config.yaml]",
        ),
        click.option(
            "--threads", help="Number of threads to use", default=1, show_default=True
        ),
        click.option(
            "--profile",
            default=None,
            help="Snakemake profile to use",
            show_default=False,
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("workflow", "conda")),
            help="Default conda env prefix directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--log",
            default="assembly_finder.log",
            callback=default_to_output,
            hidden=True,
        ),
        click.option(
            "--system-config",
            default=snake_base(os.path.join("config", "config.yaml")),
            hidden=True,
        ),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.OPTION_GROUPS = {
    "assembly_finder": [
        {
            "name": "Options",
            "options": [
                "--input",
                "--output",
                "--taxonkit",
                "--threads",
                "--requests",
                "--taxon",
                "--rank",
                "--nrank",
                "--print-versions",
            ],
        },
        {
            "name": "NCBI datasets options",
            "options": [
                "--api-key",
                "--limit",
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
            "name": "Snakemake options",
            "options": [
                "--configfile",
                "--profile",
                "--use-conda",
                "--conda-prefix",
                "--snake-default",
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
@common_options
@click.version_option(get_version(), "-v", "--version", is_flag=True)
@click.option(
    "-i",
    "--input",
    type=str,
    help="Comma seperated queries or input file",
    required=True,
)
@click.option(
    "--taxonkit",
    default=lambda: os.path.join(os.getcwd(), ".taxonkit"),
    help="Define path to taxonkit data-dir",
    type=click.Path(),
)
@click.option(
    "-nb",
    "--limit",
    help="Limit number of genomes per query",
    type=str,
    default=None,
)
@click.option("--api-key", type=str, help="NCBI api-key", default=None)
@click.option(
    "--requests",
    type=int,
    help="Number of NCBI datasets commands to run in parallel",
    default=1,
    show_default=True,
)
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
    help="Comma seperated files to download : genome,rna,protein,cds,gff3,gtf,gbff,seq-report",
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
    "--taxon/--accession",
    help="Are queries taxa names or accession",
    type=bool,
    default=True,
    show_default=True,
)
@click.option(
    "--reference",
    type=bool,
    help="Limit to reference and representative genomes",
    default=True,
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
    type=click.Choice(["exclude", "all", "only"], case_sensitive=False),
    help="Exclude, include or limit to metagenome assembled genomes",
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
def main(**kwargs):
    """
    \b
     ░█▀█░█▀▀░█▀▀░█▀▀░█▄█░█▀▄░█░░░█░█░░░█▀▀░▀█▀░█▀█░█▀▄░█▀▀░█▀▄
     ░█▀█░▀▀█░▀▀█░█▀▀░█░█░█▀▄░█░░░░█░░░░█▀▀░░█░░█░█░█░█░█▀▀░█▀▄
     ░▀░▀░▀▀▀░▀▀▀░▀▀▀░▀░▀░▀▀░░▀▀▀░░▀░░░░▀░░░▀▀▀░▀░▀░▀▀░░▀▀▀░▀░▀
    \b
    Snakemake-powered cli to download genomes with NCBI datasets
    """

    merge_config = {"args": kwargs}

    run_snakemake(
        snakefile_path=snake_base(os.path.join("workflow", "Snakefile")),
        merge_config=merge_config,
        **kwargs,
    )


if __name__ == "__main__":
    main()
