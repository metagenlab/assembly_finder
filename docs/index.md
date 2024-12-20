# assembly_finder
[![tests](https://github.com/metagenlab/assembly_finder/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/metagenlab/assembly_finder/actions/workflows/unit-tests.yml)
[![docs](https://github.com/metagenlab/assembly_finder/actions/workflows/build-docs.yml/badge.svg)](https://github.com/metagenlab/assembly_finder/actions/workflows/build-docs.yml)
[![docker](https://github.com/metagenlab/assembly_finder/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/metagenlab/assembly_finder/actions/workflows/docker-publish.yml)

[![snaketool](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
[![license](https://img.shields.io/github/license/metagenlab/assembly_finder.svg)](https://github.com/metagenlab/assembly_finder/blob/main/LICENSE)
[![version](https://img.shields.io/conda/vn/bioconda/assembly_finder)](http://bioconda.github.io/recipes/assembly_finder/README.html)
[![downloads](https://img.shields.io/conda/dn/bioconda/assembly_finder)](https://anaconda.org/bioconda/assembly_finder)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13353494.svg)](https://zenodo.org/doi/10.5281/zenodo.13353494)

`assembly_finder` is a [Snakemake](https://github.com/snakemake/snakemake)-powered cli, written in [Snaketool](https://github.com/beardymcjohnface/Snaketool), to download genomes with [NCBI datasets](https://github.com/ncbi/datasets).  

## Installation

=== "Conda <small>(recommended)</small>" 

    ```sh
    conda create -n assembly_finder assembly_finder
    ```
    !!! note
        Requires a [Miniforge](https://github.com/conda-forge/miniforge) installation

=== "Apptainer" 

    ```sh
    apptainer pull docker://ghcr.io/metagenlab/assembly_finder:latest
    ```

    !!! note
        Add `--no-use-conda` when using the container

=== "git" 

    ```sh
    git clone https://github.com/metagenlab/assembly_finder.git
    pip install -e assembly_finder
    ```
    !!! note
        Requires a [Miniforge](https://github.com/conda-forge/miniforge) installation

## Usage 
### Command

=== "standard"

    ```sh
    assembly_finder -i staphylococcus_aureus -nb 1 
    ```

=== "container"

    ```sh
    apptainer run docker://ghcr.io/metagenlab/assembly_finder:latest \
    assembly_finder -i staphylococcus_aureus -nb 1 --no-use-conda
    ```

### Output
```sh
📂staphylococcus_aureus
 ┣ 📂download
 ┃ ┣ 📂GCF_000013425.1
 ┃ ┃ ┗ 📜GCF_000013425.1_ASM1342v1_genomic.fna.gz
 ┃ ┗ 📜.snakemake_timestamp
 ┣ 📂logs
 ┃ ┣ 📂taxons
 ┃ ┃ ┗ 📜staphylococcus_aureus.log
 ┃ ┣ 📜archive.log
 ┃ ┣ 📜lineage.log
 ┃ ┣ 📜rsync.log
 ┃ ┗ 📜unzip.log
 ┣ 📜archive.zip
 ┣ 📜assembly_finder.log
 ┣ 📜assembly_summary.tsv
 ┣ 📜config.yaml
 ┗ 📜taxonomy.tsv
```

## Command-line options

![`assembly_finder -h`](images/af-help.svg)