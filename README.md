# Assembly Finder
Assembly finder is a snakemake workflow used for downloading genomes from NCBI's assembly database.

## Installation
#### Conda
Download and install miniconda 3 (Linux 64-bit)
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
#### Snakemake
Create a conda environment and install the latest version of snakemake
```bash
conda create -c bioconda -c conda-forge -n snakemake snakemake
```
#### Assembly Finder
Clone the github repository
```bash
git clone https://github.com/metagenlab/assembly_finder.git
```

## Required files
Assembly finder requires as input a tsv table containing queries, and a configuration file for setting the workflow's parameters.
### Input table example
Assembly finder searches NCBI's assembly database for a number of genomes at a given taxonomic rank or taxid inputted in a tsv table as shown below.
```
UserInputNames           nb_genomes
1813735                  1
114185                   1
ATCC_13985               1
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
####E-utilities parameters
Assembly finder uses NCBI's Entrez utilities to search for genomes in the assembly database. Thus, the user is required to create an NCBI account with an API key to avoid IP address bans and have higher requests rate per second  (Go to http://www.ncbi.nlm.nih.gov/account/ to sign-up).
####Search filter parameters
The user can expand or narrow down the set of assemblies found per query by modifying search filters. For example, parameters shown above include all Refseq and Genbank genomes except assemblies from metagenomes.
####Filtering function parameter
For each query, assembly informations are retrieved and stored in tables which are then sorted according to Refseq category, assembly status, contig count and Genbank release date by a filtering function.
The filtering function can then be set to keep one assembly per taxonomic rank. For example, after ranking the best assemblies at the top of the table, the user can choose one assembly per species by setting Rank_to_filter_by to 'species'. 
By default, no specific taxonomic rank is used to select assemblies.


## Running Assembly Finder
Below is an example command to run assembly finder with 10 cores and a maximum of 3 NCBI requests per second using the previously described required files.
```bash
conda activate snakemake
snakemake --snakefile pipeline/path/Snakefile --configfile config.yml --use-conda --conda-prefix path/to/conda/envs --resources ncbi_requests=3 --cores 10 all_download
```
Assemblies are saved in the assembly_gz/ directory.
