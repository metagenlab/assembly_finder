import subprocess
from pathlib import Path
import pytest
import shutil
import os


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


test_data_path = Path("assembly_finder/test_data")
outdir = Path("test_out")
threads = 2


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_cli():
    exec_command("assembly_finder -h")
    exec_command("assembly_finder -v")


def test_taxons():
    """download genomes from taxons"""
    exec_command(
        f"assembly_finder --threads {threads} -i {test_data_path}/taxons.tsv --output {outdir}"
    )
    remove_directory(outdir)


def test_accessions():
    """download genomes from accessions"""
    exec_command(
        f"assembly_finder --threads {threads} -i {test_data_path}/accessions.txt --accession --output {outdir}"
    )
    remove_directory(outdir)
