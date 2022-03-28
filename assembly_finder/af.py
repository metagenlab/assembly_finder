
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
    assembly_finder is a snakemake workflow used to download genome assemblies from RefSeq and Genbank
    """

@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True),
    short_help="run assembly_finder with command line arguments"
)
@click.option(
    "-i",
    "--input",
    type=click.Path(exists=True, resolve_path=True),
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
    "--complete_assemblies",
    help="download only complete assemblies (default=False)",
    default=True
)

@click.option(
    "-gr",
    "--reference_assemblies",
    help="download only reference assemblies",
    default=False
)
@click.option(
    "-gre",
    "--representative_assemblies",
    help="download only representative assemblies",
    default=False
)
@click.option(
    "-gb",
    "--genbank_assemblies",
    help="download genbank assemblies (default True)",
    is_flag=True
)
@click.option(
    "-rs",
    "--refseq_assemblies",
    help="download refseq assemblies (default True)",
    type='str'
)
@click.option(
    "-ex",
    "--exclude",
    help="exclude assemblies",
    is_flag=True,
    type=str
)
@click.option(
    "-f",
    "--filter_rank",
    help="Rank filter",
    default='',
    is_flag=False,
    type=str
)
@click.option(
    "-nr",
    "--n_by_rank",
    help="Max number of genome by target rank (eg 1/species)",
    type=int,
    default=1,
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)

def run_workflow(conda_prefix, 
                 input,
                 output,
                 dryrun_status, 
                 threads,
                 ncbi_key,
                 ncbi_email,
                 complete_assemblies,
                 reference_assemblies,
                 representative_assemblies,
                 exclude_from_metagenomes,
                 genbank_assemblies,
                 refseq_assemblies,
                 filter_rank,
                 n_by_rank,
                 snakemake_args):
    """
    Runs assembly_finder pipeline with all steps

    
    input_table: path/to/input_table
    ncbi_key: your_ncbi_api_key
    ncbi_email: your_ncbi_email
    ##Parameters for search_assemblies function
    #This set of parameters is to search all possible assemblies
    complete_assemblies: False
    reference_assemblies: False
    representative_assemblies: False
    exclude_from_metagenomes: True
    Genbank_assemblies: True
    Refseq_assemblies: True
    ##Parameters for the filtering function
    Rank_to_filter_by: False
        
    """
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
        f"complete_assemblies={complete_assemblies} "
        f"reference_assemblies={reference_assemblies} "
        f"representative_assemblies={representative_assemblies} "
        f"exclude_from_metagenomes={exclude_from_metagenomes} "
        f"Genbank_assemblies={genbank_assemblies} "
        f"Refseq_assemblies={refseq_assemblies} "
        f"Rank_to_filter_by={filter_rank} "
        f"n_by_rank={n_by_rank} "
        f"{' '.join(snakemake_args)}")
    logging.info("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)


@cli.command(
    "conf",
    context_settings=dict(ignore_unknown_options=True),
    short_help="run assembly_finder using config file"
)
@click.option(
    "-t",
    "--config-file",
    type=click.Path(exists=True, resolve_path=True),
    help="path to config file",
)
@click.option(
    "-t",
    "--threads",
    type=int,
    help="number of threads to allow for the workflow",
    default=2,
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
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)

def run_workflow_conf(config_file, threads, dryrun_status, conda_prefix, snakemake_args):
    """
    Runs assembly_finder using config file

    config: path/to/config.yml
        
    """

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
        f" --configfile {config_file} "
        f"all_download {dryrun} "
        f"{' '.join(snakemake_args)}")
    logging.info("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)


if __name__ == "__main__":
    cli()