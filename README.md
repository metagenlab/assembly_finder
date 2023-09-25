# Assembly Finder

Assembly finder is a snakemake workflow used for downloading genomes from NCBI's assembly database.

## Installation

Install with [mamba](https://github.com/mamba-org/mamba)

```sh
mamba create -c bioconda -c conda-forge -c hcc -c metagenlab -n assembly_finder assembly_finder
```

## Usage

```sh
assembly_finder -i <input> -o <outdir> -ne <ncbi_email>
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
                               [default: False]
  -t, --threads INTEGER        number of threads to allow for the workflow
                               [default: 2]
  -nk, --ncbi_key TEXT         ncbi key for Entrez
  -ne, --ncbi_email TEXT       ncbi email for Entrez
  -db, --database TEXT         download from refseq or genbank  [default:
                               refseq]
  -id, --uid                   are inputs UIDs  [default: False]
  -rc, --refseq_category TEXT  select reference and/or representative genomes
                               [default: all]
  -al, --assembly_level TEXT   select complete_genome, chromosome, scaffold or
                               contig level assemblies  [default:
                               complete_genome]
  -an, --annotation            select assemblies with annotation  [default:
                               False]
  -ex, --exclude TEXT          exclude genomes  [default: metagenome]
  -r, --filter_rank TEXT       Rank to filter by (example: species)  [default:
                               none]
  -nr, --n_by_rank TEXT        Max number of genome by target rank (example: 1
                               per species)
  -nb, --n_by_entry TEXT       Number of genomes per entry  [default: all]
  -v, --version                Show the version and exit.
  -h, --help                   Show this message and exit.

```

## Defaults

### Input type

assembly_finder assumes that inputs are either scientific names or taxids. If you want to download specific assemblies, you have to provide their UID and the -id flag.

### Number of genomes per entry

assembly_finder downloads all assemblies per entries. less can be selected by modifying the -nb flag.
(-nb 1 to select only one genome per entry)

### NCBI filters

By default assembly_finder downloads from the refseq database: reference, representative (and na) complete, annotated genomes, excluding genomes from metagenomes.
assembly_finder does not select assemblies with annotations, to do so add the -an flag.

### Taxonomy filters

To filter n assemblies from taxonomic rank (species, genus, etc...) :

```sh
-r <rank> -nr <n>
```

## Examples

### Using the [example](test.tsv)

```sh
assembly_finder -i test.tsv -o test
```

or

```sh
assembly_finder -i 1290,1813735,114185 -o test -nb 1
```

## Download the top 1 assemblies for each bacterial species

```sh
assembly_finder -i bacteria -o <outdir> -ne <ncbi_email> -nk <ncbi_key> -r species -nr 1
```

## Donwload all refseq bacteria viruses and archaea complete genomes (exclude metagenome and anomalous)

```sh
assembly_finder -i bacteria,viruses,archaea -o <outdir> -ne <ncbi_email> -nk <ncbi_key> -al complete_genome -ex metagenome,anomalous -t 3
```

## Download specific assemblies from genbank

```sh
assembly_finder -i <UID1,UID2,UID3> -o <outdir> -db genbank -uid
```

## Output

Compressed fasta files are saved in the assemblies directory, and a summary.tsv report file is generated

```sh
📂test
 ┣ 📂assemblies
 ┃ ┣ 📜GCF_001274515.1_ASM127451v1_genomic.fna.gz
 ┃ ┣ 📜GCF_001618865.1_ASM161886v1_genomic.fna.gz
 ┃ ┗ 📜GCF_003812505.1_ASM381250v1_genomic.fna.gz
 ┣ 📂benchmark
 ┃ ┣ 📜downloads.txt
 ┃ ┣ 📜find-assemblies-114185.txt
 ┃ ┣ 📜find-assemblies-1290.txt
 ┃ ┗ 📜find-assemblies-1813735.txt
 ┣ 📂logs
 ┃ ┣ 📜download.log
 ┃ ┣ 📜ete3-update.log
 ┃ ┣ 📜find-assemblies-114185.log
 ┃ ┣ 📜find-assemblies-1290.log
 ┃ ┗ 📜find-assemblies-1813735.log
 ┣ 📂tables
 ┃ ┣ 📜114185-all.tsv
 ┃ ┣ 📜114185-filtered.tsv
 ┃ ┣ 📜1290-all.tsv
 ┃ ┣ 📜1290-filtered.tsv
 ┃ ┣ 📜1813735-all.tsv
 ┃ ┗ 📜1813735-filtered.tsv
 ┗ 📜summary.tsv
 ```

summary.tsv :
| entry   | database | uid      | asm_accession   | asm_name    | path                                                                                                     | asm_status      | refseq_category       | contig_count | contig_l50 | contig_n50 | coverage | genome_size | taxid   | organism                                        | sub_type | sub_value    | superkingdom | phylum          | class               | order              | family              | genus                 | species                      |
| :------ | :------- | :------- | :-------------- | :---------- | :------------------------------------------------------------------------------------------------------- | :-------------- | :-------------------- | :----------- | :--------- | :--------- | :------- | :---------- | :------ | :---------------------------------------------- | :------- | :----------- | :----------- | :-------------- | :------------------ | :----------------- | :------------------ | :-------------------- | :--------------------------- |
| 1290    | refseq   | 7934998  | GCF_003812505.1 | ASM381250v1 | absolute/path/to/assemblies/GCF_003812505.1_ASM381250v1_genomic.fna.gz | Complete Genome | representative genome | 3            | 1          | 2220494    | 19.64    | 2257431     | 1290    | Staphylococcus hominis (firmicutes)             | strain   | FDAARGOS_575 | Bacteria     | Bacillota       | Bacilli             | Bacillales         | Staphylococcaceae   | Staphylococcus        | Staphylococcus hominis       |
| 1813735 | refseq   | 6612678  | GCF_001618865.1 | ASM161886v1 | absolute/path/to/assemblies/GCF_001618865.1_ASM161886v1_genomic.fna.gz | Complete Genome | representative genome | 1            | 1          | 7480314    | 112.0    | 7480314     | 1855912 | Luteitalea pratensis (bacteria)                 | strain   | DSM 100886   | Bacteria     | Acidobacteriota | Vicinamibacteria    | Vicinamibacterales | Vicinamibacteraceae | Luteitalea            | Luteitalea pratensis         |
| 114185  | refseq   | 15546308 | GCF_001274515.1 | ASM127451v1 | absolute/path/to/assemblies/GCF_001274515.1_ASM127451v1_genomic.fna.gz | Complete Genome | na                    | 1            | 1          | 174018     | 85.24    | 174018      | 114186  | Candidatus Carsonella ruddii (g-proteobacteria) | strain   | YCCR         | Bacteria     | Pseudomonadota  | Gammaproteobacteria | Oceanospirillales  | Halomonadaceae      | Candidatus Carsonella | Candidatus Carsonella ruddii |
