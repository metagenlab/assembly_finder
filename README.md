# Assembly Finder

Assembly finder is a snakemake workflow used for downloading genomes from NCBI's assembly database.

# Installation
Install with [mamba](https://github.com/mamba-org/mamba) 
```bash
mamba create -c bioconda -c conda-forge -c hcc -c metagenlab -n assembly_finder assembly_finder
```

# Usage
```bash
assembly_finder run -i <input> -o <outdir> -ne <ncbi_email> -nk <ncbi_key> -t <threads>
```

# Parameters
```bash
  █████╗ ███████╗███████╗███████╗███╗   ███╗██████╗ ██╗  ██╗   ██╗    ███████╗██╗███╗   ██╗██████╗ ███████╗██████╗
 ██╔══██╗██╔════╝██╔════╝██╔════╝████╗ ████║██╔══██╗██║  ╚██╗ ██╔╝    ██╔════╝██║████╗  ██║██╔══██╗██╔════╝██╔══██╗
 ███████║███████╗███████╗█████╗  ██╔████╔██║██████╔╝██║   ╚████╔╝     █████╗  ██║██╔██╗ ██║██║  ██║█████╗  ██████╔╝
 ██╔══██║╚════██║╚════██║██╔══╝  ██║╚██╔╝██║██╔══██╗██║    ╚██╔╝      ██╔══╝  ██║██║╚██╗██║██║  ██║██╔══╝  ██╔══██╗
 ██║  ██║███████║███████║███████╗██║ ╚═╝ ██║██████╔╝███████╗██║       ██║     ██║██║ ╚████║██████╔╝███████╗██║  ██║
 ╚═╝  ╚═╝╚══════╝╚══════╝╚══════╝╚═╝     ╚═╝╚═════╝ ╚══════╝╚═╝       ╚═╝     ╚═╝╚═╝  ╚═══╝╚═════╝ ╚══════╝╚═╝  ╚═╝
                                                                                         version 1.1.0
Usage: assembly_finder run [OPTIONS] [SNAKEMAKE_ARGS]...

Options:
  -i, --input TEXT             path to assembly_finder input table or list of
                               entries
  -o, --output TEXT            Output directory
  -p, --conda-prefix PATH      path to conda environment
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
  -al, --assembly_level TEXT   select complete_genome, chromosome, scaffold or
                               contig level assemblies  [default:
                               complete_genome]
  -an, --annotation            select assemblies with annotation
  -ex, --exclude TEXT          exclude genomes  [default: metagenome]
  -f, --filter_rank TEXT       Rank to filter by (example: species)  [default:
                               none]
  -nr, --n_by_rank INTEGER     Max number of genome by target rank (example: 1
                               per species)
  -nb, --n_by_entry TEXT       Number of genomes per entry
  -dl, --downloader TEXT       Use wget or aspera to download genomes
  -h, --help                   Show this message and exit.
```
# Defaults
## Input type
assembly_finder assumes that inputs are either scientific names or taxids. If you want to download specific assemblies, you have to provide their UID and the -id flag.
## Number of genomes per entry
assembly_finder downloads all assemblies per entries. less can be selected by modifying the -nb flag.
(-nb 1 to select only one genome per entry)
## NCBI filters
By default assembly_finder downloads from the refseq database: reference, representative (and na) complete, annotated genomes, excluding genomes from metagenomes. 
assembly_finder does not select assemblies with annotations, to do so add the -an flag.
## Taxonomy filters
To filter n assemblies from taxonomic rank (species, genus, etc...) :
```bash
-f <rank> -nr <n>
```
## Downloader
Batch download of files is done via [aspera-cli](https://github.com/IBM/aspera-cli).
[wget](https://www.gnu.org/software/wget) can be used instead by specifying -dl wget

# Examples
## Download the top 1 assemblies for each bacterial species
```bash
assembly_finder run -i bacteria -o <outdir> -ne <ncbi_email> -nk <ncbi_key> -f species -nr 1 
```
## Donwload all refseq bacteria viruses and archaea complete genomes (exclude metagenome and anomalous)
```bash
assembly_finder run -i 'bacteria,viruses,archaea' -o <outdir> -ne <ncbi_email> -nk <ncbi_key> -al complete_genome -ex 'metagenome,anomalous' -t 3 
```
## Download specific assemblies from genbank
```bash
assembly_finder run -i <'UID1,UID2','UID3'> -o <outdir> -db genbank -uid
```

# Outputs
Compressed fasta files are saved in the assemblies directory, and an assemblies_summary.tsv report file is generated.

assembly_summary.tsv example :
|AssemblyName|RefSeq_category|AssemblyStatus|ContigN50| AssemblyLength|ContigCount|Taxid|species
|-|-|-|-|-|-|-|-|
|GCF_000898155.1_ViralProj167578|reference genome|Complete Genome|9908|9908|1|861561|Cyrtanthus elatus virus A
|GCF_000864765.1_ViralProj15476|reference genome|Complete Genome|9181|9181|1|11676|Human immunodeficiency virus 1
