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

## :zap: Quick start 
### Installation
#### [Conda](https://github.com/conda-forge/miniforge)

```sh
conda create -n assembly_finder assembly_finder
```

> [!NOTE]  
> Miniforge is the recommended conda-based distribution


#### [Apptainer](https://github.com/apptainer/apptainer)
```sh
apptainer pull docker://ghcr.io/metagenlab/assembly_finder:latest
```

### Usage
#### Conda
```sh
assembly_finder -i staphylococcus_aureus -nb 1
```
#### Apptainer
```sh
apptainer run docker://ghcr.io/metagenlab/assembly_finder:latest \
assembly_finder -i staphylococcus_aureus -nb 1 --no-use-conda
```
> [!NOTE]  
> set --no-use-conda when running in a container

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
 ┣ 📜sequence_report.tsv
 ┗ 📜taxonomy.tsv
```

## :books: Documentation

You can find more information on assembly_finder's inputs, outputs and example commands in the [documentation](https://metagenlab.github.io/assembly_finder/)

## :scroll: Help

![`assembly_finder -h`](docs/images/af-help.svg)
