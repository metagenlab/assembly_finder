# Assembly Finder

[![version](https://img.shields.io/conda/v/bioconda/assembly_finder?label=version)](http://bioconda.github.io/recipes/assembly_finder/README.html)
[![downloads](https://img.shields.io/conda/dn/bioconda/assembly_finder)](https://anaconda.org/bioconda/assembly_finder)

Assembly finder is a Snakemake-powered cli to download genomes with [NCBI datasets](https://github.com/ncbi/datasets).

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
         v0.4.0

        Snakemake pipeline to download genome assemblies from NCBI

        github: https://github.com/metagenlab/assembly_finder

Options:
  -i, --input TEXT                path to assembly_finder input table or list
                                  of entries  [required]
  -nb, --n_by_entry TEXT          number of assemblies per entry  [default:
                                  all]
  -s, --suffixes TEXT             suffix of files to download from NCBI's ftp
                                  [default:
                                  assembly_report.txt,genomic.fna.gz]
  -o, --outdir TEXT               output directory
  -n, --dryrun_status             snakemake dryrun to see the scheduling plan
  -t, --threads INTEGER           number of threads to allow for the workflow
                                  [default: 2]
  -nk, --ncbi_key TEXT            ncbi key for Entrez
  -ne, --ncbi_email TEXT          ncbi email for Entrez
  -db, --database [refseq|genbank]
                                  download from refseq or genbank  [default:
                                  refseq]
  -id, --uid TEXT                 are inputs UIDs or assembly names  [default:
                                  False]
  -rc, --refseq_category TEXT     select reference, representative or all
                                  [default: all]
  -al, --assembly_level TEXT      select complete, chromosome, scaffold,
                                  contig or all  [default: complete]
  -an, --annotation [False|True]  select assemblies with annotation  [default:
                                  False]
  -ex, --exclude TEXT             filter to exclude assemblies (example:
                                  exclude from metagenomes)  [default:
                                  metagenome]
  -r, --rank [superkingdom|phylum|class|order|family|genus|species|none]
                                  taxonomic rank to filter by assemblies
                                  [default: none]
  -nr, --n_by_rank TEXT           max number of genome by target rank
                                  (example: 1 per species)  [default: none]
  -et, --ete_db TEXT              path where to save/find ete taxa.sqlite file
                                  [default: /home/fchaaban/.etetoolkit]
  -v, --version                   Show the version and exit.
  -h, --help                      Show this message and exit.
```

### Input

Input can be a table with entries and their respective parametes as columns (nb, rank, refseq category ...). See [minimal](minimal.tsv) and [full](full.tsv) table examples.

Additionally, the input can be a string of entries (taxids, taxonomic names or other).

#### Entry examples

taxid

```sh
assembly_finder -i 114185
```

species name

```sh
assembly_finder -i candidatus-carsonella
```

assembly accession

```sh
assembly_finder -i GCF_000287275.1
```

assembly name

```sh
assembly_finder -i ASM28727v1
```

assembly uid

```sh
assembly_finder -i 421728 -id True
```

:warning: Make sure to add the id flag, because 421728 is also a taxid !

ATCC number

```sh
assembly_finder -i ATCC_13985
```

:warning: Using entries such as ATCC strain number is not as precise as using taxids, taxonomic names or assembly accessions/names and can give unexpected results.

### Suffixes

Option to set which files to download from NCBI's ftp.

#### Suffix examples

Download assembly reports only

```sh
assembly_finder -i 114185 -s assembly_report.txt
```

Download reports, fasta and gff files

```sh
assembly_finder -i 114185 -s genomic.fna.gz,genomic.gff.gz,assembly_report.txt
```

### Assembly parameters

#### Assembly level

Select genomes assembled at the chromsome level only

```sh
assembly_finder -i 114185 -al chromosome
```

You can combine multiple levels using underscores

```sh
assembly_finder -i 114185 -al chromosome_scaffold_contig
```

For more information on [assembly levels](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/glossary/)

> - Complete genome: All chromosomes are gapless and contain runs of nine or less ambiguous bases (Ns), there are no unplaced or unlocalized scaffolds, and all the expected chromosomes are present (i.e., the assembly is not noted as having partial genome representation). Plasmids and organelles may or may not be included in the assembly, but if they are present, the sequences are gapless.

> - Chromosome: There is a sequence for one or more chromosomes. This may be a completely sequenced chromosome without gaps or a chromosome containing scaffolds or contigs with gaps between them. There may also be unplaced or unlocalized scaffolds.\*

> - Contig: Nothing is assembled beyond the level of sequence contigs.

> - Scaffold: Some sequence contigs have been connected across gaps to create scaffolds, but the scaffolds are all unplaced or unlocalized.

#### Refseq Category

Select reference genomes only

```sh
assembly_finder -i 114185 -rc reference
```

Select reference and representative genomes only

```sh
assembly_finder -i 114185 -rc reference_representative
```

No refseq category selection

```sh
assembly_finder -i 114185 -rc all
```

More on [refseq categories](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/glossary/) :

> - Reference genome: a manually selected high quality genome assembly that NCBI and the community have identified as being important as a standard against which other data are compared
> - Representative genome: a genome computationally or manually selected as a representative from among the best genomes available for a species or clade that does not have a designated reference genome

#### Exclude

Option to use exclude filters.

Exclude anomalous, metagenome and low gene count assemblies

```sh
assembly_finder -i 114185 -ex metagenome_anomalous_low-gene-count
```

### Taxonomy parameters

#### rank and number of assemblies per rank

Options to select n assemblies at a specific taxonomic rank

Download all complete assemblies for each chlamydia sepcies

```sh
assembly_finder -i chlamydia -r species
```

Download the top 1 assemby per chlamydia species

```sh
assembly_finder -i chlamydia -r species -nr 1
```
