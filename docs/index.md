# Assembly Finder

[![version](https://img.shields.io/conda/v/bioconda/assembly_finder?label=version)](http://bioconda.github.io/recipes/assembly_finder/README.html)
[![downloads](https://img.shields.io/conda/dn/bioconda/assembly_finder)](https://anaconda.org/bioconda/assembly_finder)

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
 ┃ ┣ 📂GCF_000418345.1
 ┃ ┃ ┗ 📜GCF_000418345.1_ASM41834v1_genomic.fna.gz
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

![`assembly_finder -h`](images/af-help.svg)