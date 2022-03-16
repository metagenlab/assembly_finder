# Assembly Finder
Assembly finder is a snakemake workflow used for downloading genomes from NCBI's assembly database.

## Installation
#### Conda
Download and install miniconda 3 (Linux 64-bit)
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
#### Assembly Finder
Install with conda
```bash
conda install -c bioconda -c conda-forge -c hcc assembly_finder
```

## Required files
Assembly finder requires as input a tsv table containing queries, and a configuration file for setting the workflow's parameters.

### Input table example
Assembly finder searches NCBI's assembly database for a number of genomes at a given taxonomic rank or taxid inputted in a tsv table as shown below.
```
TaxnonomyInput           NbGenomes
1813735                  1
114185                   1
staphylococcus_aureus    1
```
The workflow accepts taxonomy identifiers, taxonomic ranks (like species names) and strain numbers as queries as long as they are found in NCBI's taxonomy database. 

### Config file example
```
input_table_path: path/to/input_table

##E-utilities parameters
NCBI_key: your_ncbi_api_key
NCBI_email: your_ncbi_email

##Search filter parameters
complete_assemblies: False
reference_assemblies: False
representative_assemblies: False
exclude_from_metagenomes: True
Genbank_assemblies: True
Refseq_assemblies: True

##Filtering function parameter
Rank_to_filter_by: 'None'
```

#### E-utilities parameters
Assembly finder uses NCBI's Entrez utilities to search for genomes in the assembly database. Thus, the user is required to create an NCBI account with an API key to avoid IP address bans and have higher requests rate per second  (Go to http://www.ncbi.nlm.nih.gov/account/ to sign-up).

#### Search filter parameters
The user can expand or narrow down the set of assemblies found per query by modifying search filters. For example, parameters shown above include all Refseq and Genbank genomes except assemblies from metagenomes.

#### Filtering function parameter
For each query, assembly informations are retrieved and stored in tables which are then sorted according to Refseq category, assembly status, contig count and Genbank release date by a filtering function.
The filtering function can then be set to keep one assembly per taxonomic rank. For example, after ranking the best assemblies at the top of the table, the user can choose one assembly per species by setting Rank_to_filter_by to 'species'. 
By default, no specific taxonomic rank is used to select assemblies.

## Running Assembly Finder using config file

```bash
Usage: af conf [OPTIONS] [SNAKEMAKE_ARGS]...

  Runs assembly_finder using config file

  config: path/to/config.yml

Options:
  -c, --config-file PATH   path to config file
  -c, --cores INTEGER      number of cores to allow for the workflow
  -p, --conda-prefix PATH  path to conda environment
  -n, --dryrun_status      Snakemake dryrun to see the scheduling plan
                           [default: False]
  -h, --help               Show this message and exit.
```

Below is an example command to run assembly finder with 10 cores and a maximum of 3 NCBI requests per second using the previously described required files.

```bash
af conf --configfile config.yml --resources ncbi_requests=3 --cores 10
```

Assemblies are saved in the assembly_gz/ directory.

## Running Assembly Finder in command line 

Assembly Finder can also be executed from the command line without configuration file. 

```
Usage: af run [OPTIONS] [SNAKEMAKE_ARGS]...

  Runs assembly_finder pipeline with all steps

  input_table_path: path/to/input_table ncbi_key: your_ncbi_api_key
  ncbi_email: your_ncbi_email ##Parameters for search_assemblies function
  #This set of parameters is to search all possible assemblies
  complete_assemblies: False reference_assemblies: False
  representative_assemblies: False exclude_from_metagenomes: True
  Genbank_assemblies: True Refseq_assemblies: True ##Parameters for the
  filtering function Rank_to_filter_by: False

Options:
  -i, --input-table PATH          path to assembly_finder input_table_path
  -p, --conda-prefix PATH         path to conda environment
  -n, --dryrun_status             Snakemake dryrun to see the scheduling plan
                                  [default: False]

  -c, --cores INTEGER             number of cores to allow for the workflow
  -nk, --ncbi_key TEXT            ncbi key for Entrez
  -ne, --ncbi_email TEXT          ncbi email for Entrez
  -gc, --complete_assemblies TEXT
                                  download only complete assemblies
                                  (default=False)

  -gr, --reference_assemblies TEXT
                                  download only reference assemblies
  -gre, --representative_assemblies TEXT
                                  download only representative assemblies
  -gb, --genbank_assemblies       download genbank assemblies (default True)
  -rs, --refseq_assemblies        download refseq assemblies (default True)
  -rs, --exclude_from_metagenomes
                                  exclude metagnomes (default True)
  -f, --filter_rank TEXT          Rank filter
  -nr, --n_by_rank INTEGER        Max number of genome by target rank (eg
                                  1/species)

  -h, --help                      Show this message and exit.

  ```
