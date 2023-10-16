# Assembly Finder

Assembly finder is a snakemake workflow for downloading genomes from [NCBI assembly](https://www.ncbi.nlm.nih.gov/assembly).

## Table of contents

- [Install](#installation)
- [Usage](#quick-usage)
- [Output](#output)
- [Examples](#examples)
- [Parameters](#parameters)

## Installation

Install with [mamba](https://github.com/mamba-org/mamba)

```sh
mamba create -c bioconda -n assembly_finder assembly_finder
```

## Usage

```sh
assembly_finder -i <input> -o <outdir> -ne <ncbi_email>
```

### Quick usage

```sh
assembly_finder -i test.tsv -o test
```

or

```sh
assembly_finder -i 1290,1813735,114185 -nb 3,2,1 -o test
```

### Output

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

#### Assembly summary

| entry   | database | db_uid   | asm_name     | organism                                                            | asm_release_date | asm_status      | refseq_category       | contig_count | contig_n50 | genome_size | coverage | asm_method                 | seq_tech                                             | path                                                            |
| :------ | :------- | :------- | :----------- | :------------------------------------------------------------------ | :--------------- | :-------------- | :-------------------- | :----------- | :--------- | :---------- | :------- | :------------------------- | :--------------------------------------------------- | :-------------------------------------------------------------- |
| 1290    | refseq   | 7934998  | ASM381250v1  | Staphylococcus hominis (firmicutes)                                 | 2018/11/26 00:00 | Complete Genome | representative genome | 3            | 2220494    | 2257431     | 19.64    | SMRT v. 2.3.0, HGAP v. 3.0 | PacBio; Illumina                                     | /path/to/assemblies/GCF_003812505.1_ASM381250v1_genomic.fna.gz  |
| 1290    | refseq   | 15589428 | ASM973011v1  | Staphylococcus hominis (firmicutes)                                 | 2019/12/09 00:00 | Complete Genome | na                    | 4            | 2178824    | 2242116     | 1980.61  | canu v. 1.4                | Pacbio; Illumina                                     | /path/to/assemblies/GCF_009730115.1_ASM973011v1_genomic.fna.gz  |
| 1290    | refseq   | 15589668 | ASM973013v1  | Staphylococcus hominis (firmicutes)                                 | 2019/12/09 00:00 | Complete Genome | na                    | 3            | 2323613    | 2372192     | 1875.96  | SMRT v. 2.3.0, HGAP v. 3   | Pacbio; Illumina                                     | /path/to/assemblies/GCF_009730135.1_ASM973013v1_genomic.fna.gz  |
| 114185  | refseq   | 421728   | ASM28727v1   | Candidatus Carsonella ruddii HC isolate Thao2000 (g-proteobacteria) | 2012/09/07 00:00 | Complete Genome | representative genome | 1            | 166163     | 166163      |          |                            |                                                      | /path/to/assemblies/GCF_000287275.1_ASM28727v1_genomic.fna.gz   |
| 1813735 | refseq   | 6612678  | ASM161886v1  | Luteitalea pratensis (bacteria)                                     | 2018/06/07 00:00 | Complete Genome | representative genome | 1            | 7480314    | 7480314     | 112.0    | HGAP SMRTPortal v. 2.3.0   | PacBio; Illumina HiSeq 2500                          | /path/to/assemblies/GCF_001618865.1_ASM161886v1_genomic.fna.gz  |
| 1813735 | refseq   | 31969758 | ASM1686548v2 | Luteitalea sp. TBR-22 (bacteria)                                    | 2022/02/22 00:00 | Complete Genome | na                    | 1            | 6468984    | 6468984     | 54.0     | HybridSPAdes v. 3.13.0     | Illumina Miseq system; Oxford Nanopore MinION system | /path/to/assemblies/GCF_016865485.1_ASM1686548v2_genomic.fna.gz |

#### Taxonomy summary

| asm_name     | organism                                                            | sub_type | sub_value    | taxid   | superkingdom | phylum          | class               | order              | family              | genus                 | species                      |
| :----------- | :------------------------------------------------------------------ | :------- | :----------- | :------ | :----------- | :-------------- | :------------------ | :----------------- | :------------------ | :-------------------- | :--------------------------- |
| ASM381250v1  | Staphylococcus hominis (firmicutes)                                 | strain   | FDAARGOS_575 | 1290    | Bacteria     | Bacillota       | Bacilli             | Bacillales         | Staphylococcaceae   | Staphylococcus        | Staphylococcus hominis       |
| ASM973011v1  | Staphylococcus hominis (firmicutes)                                 | strain   | FDAARGOS_747 | 1290    | Bacteria     | Bacillota       | Bacilli             | Bacillales         | Staphylococcaceae   | Staphylococcus        | Staphylococcus hominis       |
| ASM973013v1  | Staphylococcus hominis (firmicutes)                                 | strain   | FDAARGOS_746 | 1290    | Bacteria     | Bacillota       | Bacilli             | Bacillales         | Staphylococcaceae   | Staphylococcus        | Staphylococcus hominis       |
| ASM28727v1   | Candidatus Carsonella ruddii HC isolate Thao2000 (g-proteobacteria) | strain   | HC           | 1202538 | Bacteria     | Pseudomonadota  | Gammaproteobacteria | Oceanospirillales  | Halomonadaceae      | Candidatus Carsonella | Candidatus Carsonella ruddii |
| ASM161886v1  | Luteitalea pratensis (bacteria)                                     | strain   | DSM 100886   | 1855912 | Bacteria     | Acidobacteriota | Vicinamibacteria    | Vicinamibacterales | Vicinamibacteraceae | Luteitalea            | Luteitalea pratensis         |
| ASM1686548v2 | Luteitalea sp. TBR-22 (bacteria)                                    | strain   | TBR-22       | 2802971 | Bacteria     | Acidobacteriota | Vicinamibacteria    | Vicinamibacterales | Vicinamibacteraceae | Luteitalea            | Luteitalea sp. TBR-22        |

#### Sequence summary

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
assembly_finder -i bacteria,viruses,archaea -o bacteria_viruses_archaea -ne <ncbi_email> -nk <ncbi_key> -al complete -ex metagenome,anomalous -t 3
```

### Download specific assemblies from genbank

```sh
assembly_finder -i <UID1,UID2,UID3> -o <outdir> -db genbank -uid
```

## Parameters

```sh
Usage: assembly_finder.py [OPTIONS] [SNAKEMAKE_ARGS]...

         ░█▀█░█▀▀░█▀▀░█▀▀░█▄█░█▀▄░█░░░█░█░░░█▀▀░▀█▀░█▀█░█▀▄░█▀▀░█▀▄
         ░█▀█░▀▀█░▀▀█░█▀▀░█░█░█▀▄░█░░░░█░░░░█▀▀░░█░░█░█░█░█░█▀▀░█▀▄
         ░▀░▀░▀▀▀░▀▀▀░▀▀▀░▀░▀░▀▀░░▀▀▀░░▀░░░░▀░░░▀▀▀░▀░▀░▀▀░░▀▀▀░▀░▀
         v0.3.0

        Snakemake pipeline to download genome assemblies from NCBI

        github: https://github.com/metagenlab/assembly_finder

Options:
  -i, --input TEXT             path to assembly_finder input table or list of
                               entries  [required]
  -o, --output TEXT            Output directory
  -n, --dryrun_status          Snakemake dryrun to see the scheduling plan
  -t, --threads INTEGER        number of threads to allow for the workflow
                               [default: 2]
  -nk, --ncbi_key TEXT         ncbi key for Entrez
  -ne, --ncbi_email TEXT       ncbi email for Entrez
  -db, --database TEXT         download from refseq or genbank  [default:
                               refseq]
  -id, --uid                   are inputs UIDs
  -rc, --refseq_category TEXT  select reference and/or representative genomes
                               [default: all]
  -al, --assembly_level TEXT   select complete, chromosome, scaffold or contig
                               level assemblies  [default: complete]
  -an, --annotation            select assemblies with annotation
  -ex, --exclude TEXT          exclude genomes  [default: metagenome]
  -r, --filter_rank TEXT       Rank to filter by (example: species)  [default:
                               none]
  -nr, --n_by_rank TEXT        Max number of genome by target rank (example: 1
                               per species)  [default: none]
  -nb, --n_by_entry TEXT       Number of genomes per entry  [default: all]
  -et, --ete_db TEXT           path where to save/find ete taxa.sqlite file
  -v, --version                Show the version and exit.
  -h, --help                   Show this message and exit.

```

#### Input type

assembly_finder assumes that inputs are either scientific names or taxids. If you want to download specific assemblies, you have to provide their UID and the -id flag.

#### Number of genomes per entry

assembly_finder downloads all assemblies per entries. less can be selected by modifying the -nb flag.
(-nb 1 to select only one genome per entry)

#### NCBI filters

By default assembly_finder downloads from the refseq database: reference, representative (and na) complete, annotated genomes, excluding genomes from metagenomes.
assembly_finder does not select assemblies with annotations, to do so add the -an flag.

#### Taxonomy filters

To filter n assemblies from taxonomic rank (species, genus, etc...).
Example: filter the best 10 assemblies from a species:

```sh
-r species -nr 10
```
