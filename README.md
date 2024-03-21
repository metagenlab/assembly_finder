# Assembly Finder

[![version](https://img.shields.io/conda/v/bioconda/assembly_finder?label=version)](http://bioconda.github.io/recipes/assembly_finder/README.html)
[![downloads](https://img.shields.io/conda/dn/bioconda/assembly_finder)](https://anaconda.org/bioconda/assembly_finder)

Assembly finder is a Snakemake cli to download genomes with [NCBI datasets](https://github.com/ncbi/datasets).

## Installation

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/assembly_finder/README.html)

```sh
mamba create -n assembly_finder assembly_finder
```

## Usage

```sh
assembly_finder -i <input>
```

### Quick usage

```sh
assembly_finder -i taxons.txt
```

or

```sh
assembly_finder -i 1290,staphylococcus_aureus,287 -nb 1 -o taxons
```

## Output

This is how the taxons download directory looks like

```sh
📂taxons
 ┣ 📂download
 ┃ ┣ 📂GCF_000233495.1
 ┃ ┃ ┣ 📜GCF_000233495.1_Pseudomonas_sp_2_1_26_V1_genomic.fna.gz
 ┃ ┃ ┗ 📜sequence_report.jsonl
 ┃ ┣ 📂GCF_000418345.1
 ┃ ┃ ┣ 📜GCF_000418345.1_ASM41834v1_genomic.fna.gz
 ┃ ┃ ┗ 📜sequence_report.jsonl
 ┃ ┣ 📂GCF_003812505.1
 ┃ ┃ ┣ 📜GCF_003812505.1_ASM381250v1_genomic.fna.gz
 ┃ ┃ ┗ 📜sequence_report.jsonl
 ┃ ┣ 📜.snakemake_timestamp
 ┃ ┣ 📜assembly_data_report.jsonl
 ┃ ┗ 📜dataset_catalog.json
 ┣ 📜archive.zip
 ┣ 📜assembly_summary.tsv
 ┣ 📜sequence_report.tsv
 ┗ 📜taxonomy.tsv
```

## Examples

### Download the top 1 assemblies for each bacterial species

```sh
assembly_finder -i bacteria --api-key <api-key> --rank species --nrank 1
```

### Donwload all refseq bacteria viruses and archaea complete genomes (exclude MAGs and anomalous genomes)

```sh
assembly_finder -i bacteria,viruses,archaea -o outdir --api-key <api-key> --assembly-level complete --exclude-atypical --mag exclude
```

### Download specific genomes with accessions

```sh
assembly_finder -i GCF_003812505.1,GCF_000418345.1,GCF_000233495.1 -o accessions
```

## Parameters

```sh

 Usage: assembly_finder [OPTIONS] [SNAKEMAKE_ARGS]...


  ░█▀█░█▀▀░█▀▀░█▀▀░█▄█░█▀▄░█░░░█░█░░░█▀▀░▀█▀░█▀█░█▀▄░█▀▀░█▀▄
  ░█▀█░▀▀█░▀▀█░█▀▀░█░█░█▀▄░█░░░░█░░░░█▀▀░░█░░█░█░█░█░█▀▀░█▀▄
  ░▀░▀░▀▀▀░▀▀▀░▀▀▀░▀░▀░▀▀░░▀▀▀░░▀░░░░▀░░░▀▀▀░▀░▀░▀▀░░▀▀▀░▀░▀

 Snakemake-cli to download genomes with NCBI datasets

╭─ Options ──────────────────────────────────────────────────────────────────────────╮
│ *  --input    -i   TEXT                            path to assembly_finder input   │
│                                                    table or list of queries        │
│                                                    [required]                      │
│    --outdir   -o   PATH                            output directory                │
│    --number   -nb  INTEGER                         number of assemblies per query  │
│    --threads  -t   INTEGER                         number of threads to allow for  │
│                                                    the workflow                    │
│                                                    [default: 2]                    │
│    --taxon                                         are inputs taxon names or ids   │
│                                                    [default: True]                 │
│    --rank          [superkingdom|phylum|class|ord  taxonomic rank to filter by     │
│                    er|family|genus|species]        assemblies                      │
│    --nrank         INTEGER                         number of genomes per taxonomic │
│                                                    rank                            │
│    --dryrun                                        snakemake dryrun to see the     │
│                                                    scheduling plan                 │
╰────────────────────────────────────────────────────────────────────────────────────╯
╭─ NCBI datasets options ────────────────────────────────────────────────────────────╮
│ --api-key           TEXT                  NCBI api-key                             │
│ --compressed                              Download compressed files                │
│ --source            [refseq|genbank|all]  download from refseq or genbank          │
│                                           [default: all]                           │
│ --include           TEXT                  Comma seperated files to download :      │
│                                           genome,rna,protein,cds,gff3,gtf,gbff,se… │
│                                           [default: genome,seq-report]             │
│ --reference                               limit to reference and representative    │
│                                           genomes                                  │
│                                           [default: False]                         │
│ --assembly-level    TEXT                  select complete, chromosome, scaffold,   │
│                                           contig                                   │
│ --annotated                               select assemblies with annotation        │
│                                           [default: False]                         │
│ --atypical                                exclude atypical genomes [default: True] │
│ --mag               [only|exclude|all]    exclude or add MAGs to the dwnloaded     │
│                                           genomes                                  │
│                                           [default: all]                           │
╰────────────────────────────────────────────────────────────────────────────────────╯
╭─ Help ─────────────────────────────────────────────────────────────────────────────╮
│ --help     -h    Show this message and exit.                                       │
│ --version  -v    Show the version and exit.                                        │
╰────────────────────────────────────────────────────────────────────────────────────╯
```
