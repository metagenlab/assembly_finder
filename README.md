# Assembly Finder

[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
[![license](https://img.shields.io/github/license/metagenlab/assembly_finder.svg)](https://github.com/metagenlab/assembly_finder/blob/main/LICENSE)
[![version](https://img.shields.io/conda/v/bioconda/assembly_finder?label=version)](http://bioconda.github.io/recipes/assembly_finder/README.html)
[![downloads](https://img.shields.io/conda/dn/bioconda/assembly_finder)](https://anaconda.org/bioconda/assembly_finder)
[![tests](https://github.com/metagenlab/assembly_finder/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/metagenlab/assembly_finder/actions/workflows/linux-tests.yml)
[![docs](https://github.com/metagenlab/assembly_finder/actions/workflows/build-docs.yml/badge.svg)](https://github.com/metagenlab/assembly_finder/actions/workflows/build-docs.yml)

Assembly finder is a Snakemake-powered cli to download genomes with [NCBI datasets](https://github.com/ncbi/datasets).

## Installation

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/assembly_finder/README.html)

```sh
mamba create -n assembly_finder assembly_finder
```

## Example

### Command

```sh
assembly_finder -i staphylococcus_aureus -nb 1
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
 ┣ 📜sequence_report.tsv
 ┗ 📜taxonomy.tsv
```

## Usage

![`assembly_finder -h`](docs/images/af-help.svg)
