# Assembly Finder

[![version](https://img.shields.io/conda/v/bioconda/assembly_finder?label=bioconda)](http://bioconda.github.io/recipes/assembly_finder/README.html)
[![Downloads](https://img.shields.io/conda/dn/bioconda/assembly_finder)](https://anaconda.org/bioconda/assembly_finder)

Assembly finder is a snakemake workflow for downloading genomes from [NCBI assembly](https://www.ncbi.nlm.nih.gov/assembly).

## Table of contents

- [Installation](#installation)
- [Usage](#quick-usage)
- [Output](#output)
- [Examples](#examples)
- [Parameters](#parameters)

## Installation

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/assembly_finder/README.html)

```sh
mamba create -n assembly_finder assembly_finder
```

## Usage

```sh
assembly_finder -i <input> -o <outdir> -ne <ncbi_email> -nk <ncbi_key>
```

### Quick usage
```sh
assembly_finder -i 1290,1813735,114185 -nb 3,2,1 -o test
```

## Output

Compressed fasta files are saved in the assemblies directory with assembly, sequence and taxonomy summary tables in assemblies' parent directory:

```sh
📂 test
 ┣ 📂 assemblies
 ┃ ┣ 📜GCF_000287275.1_ASM28727v1_genomic.fna.gz
 ┃ ┣ 📜GCF_001618865.1_ASM161886v1_genomic.fna.gz
 ┃ ┣ 📜GCF_003812505.1_ASM381250v1_genomic.fna.gz
 ┃ ┣ 📜GCF_009730115.1_ASM973011v1_genomic.fna.gz
 ┃ ┣ 📜GCF_009730135.1_ASM973013v1_genomic.fna.gz
 ┃ ┗ 📜GCF_016865485.1_ASM1686548v2_genomic.fna.gz
 ┣ 📂 benchmark
 ┃ ┣ 📜downloads.txt
 ┃ ┣ 📜find-assemblies-114185.txt
 ┃ ┣ 📜find-assemblies-1290.txt
 ┃ ┗ 📜find-assemblies-1813735.txt
 ┣ 📂 logs
 ┃ ┣ 📜download.log
 ┃ ┣ 📜find-assemblies-114185.log
 ┃ ┣ 📜find-assemblies-1290.log
 ┃ ┗ 📜find-assemblies-1813735.log
 ┣ 📂 tables
 ┃ ┣ 📜114185-all.tsv
 ┃ ┣ 📜114185-filtered.tsv
 ┃ ┣ 📜1290-all.tsv
 ┃ ┣ 📜1290-filtered.tsv
 ┃ ┣ 📜1813735-all.tsv
 ┃ ┗ 📜1813735-filtered.tsv
 ┣ 📜assembly_summary.tsv
 ┣ 📜sequence_summary.tsv
 ┗ 📜taxonomy_summary.tsv
```

### Assembly summary

| entry   | database | db_uid   | asm_name     | organism                                                            | asm_release_date | asm_status      | refseq_category       | contig_count | contig_n50 | genome_size | coverage | asm_method                 | seq_tech                                             | path                                                            |
| :------ | :------- | :------- | :----------- | :------------------------------------------------------------------ | :--------------- | :-------------- | :-------------------- | :----------- | :--------- | :---------- | :------- | :------------------------- | :--------------------------------------------------- | :-------------------------------------------------------------- |
| 1290    | refseq   | 7934998  | ASM381250v1  | Staphylococcus hominis (firmicutes)                                 | 2018/11/26 00:00 | Complete Genome | representative genome | 3            | 2220494    | 2257431     | 19.64    | SMRT v. 2.3.0, HGAP v. 3.0 | PacBio; Illumina                                     | /path/to/assemblies/GCF_003812505.1_ASM381250v1_genomic.fna.gz  |
| 1290    | refseq   | 15589428 | ASM973011v1  | Staphylococcus hominis (firmicutes)                                 | 2019/12/09 00:00 | Complete Genome | na                    | 4            | 2178824    | 2242116     | 1980.61  | canu v. 1.4                | Pacbio; Illumina                                     | /path/to/assemblies/GCF_009730115.1_ASM973011v1_genomic.fna.gz  |
| 1290    | refseq   | 15589668 | ASM973013v1  | Staphylococcus hominis (firmicutes)                                 | 2019/12/09 00:00 | Complete Genome | na                    | 3            | 2323613    | 2372192     | 1875.96  | SMRT v. 2.3.0, HGAP v. 3   | Pacbio; Illumina                                     | /path/to/assemblies/GCF_009730135.1_ASM973013v1_genomic.fna.gz  |
| 114185  | refseq   | 421728   | ASM28727v1   | Candidatus Carsonella ruddii HC isolate Thao2000 (g-proteobacteria) | 2012/09/07 00:00 | Complete Genome | representative genome | 1            | 166163     | 166163      |          |                            |                                                      | /path/to/assemblies/GCF_000287275.1_ASM28727v1_genomic.fna.gz   |
| 1813735 | refseq   | 6612678  | ASM161886v1  | Luteitalea pratensis (bacteria)                                     | 2018/06/07 00:00 | Complete Genome | representative genome | 1            | 7480314    | 7480314     | 112.0    | HGAP SMRTPortal v. 2.3.0   | PacBio; Illumina HiSeq 2500                          | /path/to/assemblies/GCF_001618865.1_ASM161886v1_genomic.fna.gz  |
| 1813735 | refseq   | 31969758 | ASM1686548v2 | Luteitalea sp. TBR-22 (bacteria)                                    | 2022/02/22 00:00 | Complete Genome | na                    | 1            | 6468984    | 6468984     | 54.0     | HybridSPAdes v. 3.13.0     | Illumina Miseq system; Oxford Nanopore MinION system | /path/to/assemblies/GCF_016865485.1_ASM1686548v2_genomic.fna.gz |

### Taxonomy summary

| asm_name     | organism                                                            | sub_type | sub_value    | taxid   | superkingdom | phylum          | class               | order              | family              | genus                 | species                      |
| :----------- | :------------------------------------------------------------------ | :------- | :----------- | :------ | :----------- | :-------------- | :------------------ | :----------------- | :------------------ | :-------------------- | :--------------------------- |
| ASM381250v1  | Staphylococcus hominis (firmicutes)                                 | strain   | FDAARGOS_575 | 1290    | Bacteria     | Bacillota       | Bacilli             | Bacillales         | Staphylococcaceae   | Staphylococcus        | Staphylococcus hominis       |
| ASM973011v1  | Staphylococcus hominis (firmicutes)                                 | strain   | FDAARGOS_747 | 1290    | Bacteria     | Bacillota       | Bacilli             | Bacillales         | Staphylococcaceae   | Staphylococcus        | Staphylococcus hominis       |
| ASM973013v1  | Staphylococcus hominis (firmicutes)                                 | strain   | FDAARGOS_746 | 1290    | Bacteria     | Bacillota       | Bacilli             | Bacillales         | Staphylococcaceae   | Staphylococcus        | Staphylococcus hominis       |
| ASM28727v1   | Candidatus Carsonella ruddii HC isolate Thao2000 (g-proteobacteria) | strain   | HC           | 1202538 | Bacteria     | Pseudomonadota  | Gammaproteobacteria | Oceanospirillales  | Halomonadaceae      | Candidatus Carsonella | Candidatus Carsonella ruddii |
| ASM161886v1  | Luteitalea pratensis (bacteria)                                     | strain   | DSM 100886   | 1855912 | Bacteria     | Acidobacteriota | Vicinamibacteria    | Vicinamibacterales | Vicinamibacteraceae | Luteitalea            | Luteitalea pratensis         |
| ASM1686548v2 | Luteitalea sp. TBR-22 (bacteria)                                    | strain   | TBR-22       | 2802971 | Bacteria     | Acidobacteriota | Vicinamibacteria    | Vicinamibacterales | Vicinamibacteraceae | Luteitalea            | Luteitalea sp. TBR-22        |

### Sequence summary

| asm_name     | organism                                                            | taxid   | genbank_accn | refseq_accn   | assigned_molecule_location/type | sequence_length |
| :----------- | :------------------------------------------------------------------ | :------ | :----------- | :------------ | :------------------------------ | :-------------- |
| ASM381250v1  | Staphylococcus hominis (firmicutes)                                 | 1290    | CP033732.1   | NZ_CP033732.1 | Chromosome                      | 2220494         |
| ASM381250v1  | Staphylococcus hominis (firmicutes)                                 | 1290    | CP033731.1   | NZ_CP033731.1 | Plasmid                         | 32498           |
| ASM381250v1  | Staphylococcus hominis (firmicutes)                                 | 1290    | CP033733.1   | NZ_CP033733.1 | Plasmid                         | 4439            |
| ASM973011v1  | Staphylococcus hominis (firmicutes)                                 | 1290    | CP046301.1   | NZ_CP046301.1 | Chromosome                      | 2178824         |
| ASM973011v1  | Staphylococcus hominis (firmicutes)                                 | 1290    | CP046300.1   | NZ_CP046300.1 | Plasmid                         | 14283           |
| ASM973011v1  | Staphylococcus hominis (firmicutes)                                 | 1290    | CP046302.1   | NZ_CP046302.1 | Plasmid                         | 46602           |
| ASM973011v1  | Staphylococcus hominis (firmicutes)                                 | 1290    | CP046303.1   | NZ_CP046303.1 | Plasmid                         | 2407            |
| ASM973013v1  | Staphylococcus hominis (firmicutes)                                 | 1290    | CP046306.1   | NZ_CP046306.1 | Chromosome                      | 2323613         |
| ASM973013v1  | Staphylococcus hominis (firmicutes)                                 | 1290    | CP046304.1   | NZ_CP046304.1 | Plasmid                         | 23366           |
| ASM973013v1  | Staphylococcus hominis (firmicutes)                                 | 1290    | CP046305.1   | NZ_CP046305.1 | Plasmid                         | 25213           |
| ASM28727v1   | Candidatus Carsonella ruddii HC isolate Thao2000 (g-proteobacteria) | 1202538 | CP003543.1   | NC_018416.1   | Chromosome                      | 166163          |
| ASM161886v1  | Luteitalea pratensis (bacteria)                                     | 1855912 | CP015136.1   | NZ_CP015136.1 | Chromosome                      | 7480314         |
| ASM1686548v2 | Luteitalea sp. TBR-22 (bacteria)                                    | 2802971 | AP024452.2   | NZ_AP024452.1 | Chromosome                      | 6468984         |

## Examples

### Download the top 1 assemblies for each bacterial species

```sh
assembly_finder -i bacteria -o bacteria -ne <ncbi_email> -nk <ncbi_key> -r species -nr 1
```

### Donwload all refseq bacteria viruses and archaea complete genomes (exclude metagenome and anomalous)

```sh
assembly_finder -i bacteria,viruses,archaea -o outdir -ne <ncbi_email> -nk <ncbi_key> -al complete -ex metagenome_anomalous -t 3
```

### Download specific assemblies from genbank

```sh
assembly_finder -i <UID1,UID2,UID3> -o <outdir> -db genbank -uid True
```

## Parameters

```sh
Usage: assembly_finder.py [OPTIONS] [SNAKEMAKE_ARGS]...

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
assembly_finder -i 114185 -al reference_representative
```

No refseq category selection

```sh
assembly_finder -i 114185 -al all
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
