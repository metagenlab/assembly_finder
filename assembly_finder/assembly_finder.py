
#  Parts of this wrapper are inspired by atlas, an easy-to-use metagenomic pipeline based on snakemake.
#  Go check it out on https://github.com/metagenome-atlas/atlas
import logging
import os
import sys
import click
import subprocess
from assembly_finder import __version__

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M",
    format="[%(asctime)s %(levelname)s] %(message)s",
)


def get_snakefile():
    thisdir = os.path.abspath(os.path.dirname(__file__))
    sf = os.path.join(thisdir, 'Snakefile')
    if not os.path.exists(sf):
        sys.exit(f"Unable to locate the Snakemake workflow file at {sf}")
    return sf


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """
    ┌─┐┌─┐┌─┐┌─┐┌┬┐┌┐ ┬ ┬ ┬  ┌─┐┬┌┐┌┌┬┐┌─┐┬─┐
    ├─┤└─┐└─┐├┤ │││├┴┐│ └┬┘  ├┤ ││││ ││├┤ ├┬┘
    ┴ ┴└─┘└─┘└─┘┴ ┴└─┘┴─┘┴   └  ┴┘└┘─┴┘└─┘┴└─
    github: https://github.com/metagenlab/assembly_finder
    """
@cli.command(
    'run',
    context_settings=dict(ignore_unknown_options=True),
    short_help="run assembly_finder with command line arguments"
)
@click.option(
    "-i",
    "--input",
    type=str,
    help="path to assembly_finder input_table_path or list of entries",
)
@click.option(
    "-o",
    "--output",
    help="Output directory",
    type=str,
)
@click.option(
    "-p",
    "--conda-prefix",
    type=click.Path(exists=True, resolve_path=True),
    help="path to conda environment",
)
@click.option(
    "-n",
    "--dryrun_status",
    is_flag=True,
    default=False,
    show_default=True,
    help="Snakemake dryrun to see the scheduling plan",
)
@click.option(
    "-t",
    "--threads",
    type=int,
    help="number of threads to allow for the workflow",
    default=2,
    show_default=True
)
@click.option(
    "-nk",
    "--ncbi_key",
    type=str,
    help="ncbi key for Entrez",
    default="",
)
@click.option(
    "-ne",
    "--ncbi_email",
    type=str,
    help="ncbi email for Entrez",
    default="",
)
@click.option(
    "-db",
    "--database",
    type=str,
    help="download from refseq or genbank",
    default='refseq',
    show_default=True
)
@click.option(
    "-id",
    "--uid",
    help="are inputs UIDs",
    is_flag=True,
    default=False,
    show_default=True
)
@click.option(
    "-rc",
    "--refseq_category",
    type=str,
    help="select reference and/or representative genomes",
    default='all',
    show_default=True
)
@click.option(
    "-al",
    "--assembly_level",
    type=str,
    help="select complete_genome, chromosome, scaffold or contig level assemblies",
    default='complete_genome',
    show_default=True
)
@click.option(
    "-an",
    "--annotation",
    help="select assemblies with annotation",
    is_flag=True,
    default=True,
    show_default=True
)
@click.option(
    "-ex",
    "--exclude",
    type=str,
    help="exclude genomes",
    default='metagenome',
    show_default=True
)

@click.option(
    "-f",
    "--filter_rank",
    help="Rank to filter by (example: species)",
    default='none',
    is_flag=False,
    type=str,
    show_default=True
)
@click.option(
    "-nr",
    "--n_by_rank",
    help="Max number of genome by target rank (example: 1 per species)",
    type=int,
    default=1,
)
@click.option(
    "-nb",
    "--n_by_entry",
    help="Number of genomes per entry",
    type=str,
    default='all',
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)

def run_workflow(conda_prefix, 
                 input,
                 output,
                 dryrun_status, 
                 threads,
                 ncbi_key,
                 ncbi_email,
                 database,
                 uid,
                 refseq_category,
                 assembly_level,
                 annotation,
                 exclude,
                 filter_rank,
                 n_by_rank,
                 n_by_entry,
                 snakemake_args):
    import datetime
    if not output:
        output = datetime.datetime.today().strftime("%Y-%m-%d")

    if dryrun_status:
        dryrun = '-n'
    else:
        dryrun = ''
    if conda_prefix is None:
        conda_prefix = os.environ['CONDA_PREFIX']
    if not os.path.exists(conda_prefix):
        logging.critical(f"conda env path not found: {conda_prefix}")
        sys.exit(1)
    cmd = (
        f"snakemake --snakefile {get_snakefile()} --use-conda --conda-prefix {conda_prefix} "
        f" --cores {threads} "
        f"all_download {dryrun} "
        f"--config input={input} " 
        f"NCBI_key={ncbi_key} "
        f"NCBI_email={ncbi_email} "
        f"outdir={output} "
        f"db={database} "
        f"uid={uid} "
        f"refseq_category={refseq_category} "
        f"assembly_level={assembly_level} "
        f"annotation={annotation} "
        f"exclude={exclude} "
        f"Rank_to_filter_by={filter_rank} "
        f"n_by_rank={n_by_rank} "
        f"n_by_entry={n_by_entry} "
        f"{' '.join(snakemake_args)}")
    logging.info(f"Executing: {cmd}")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)

if __name__ == "__main__":
    cli()