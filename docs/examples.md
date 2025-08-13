## Download summary tables

Starting from [v0.8.0](https://github.com/metagenlab/assembly_finder/releases/tag/v0.8.0), you can restrict outputs to `assembly_summary.tsv` and `taxonomy.tsv`

* Command

```sh
assembly_finder -i staphylococcus_aureus --reference --summary
```

* Output

```sh
📂staphylococcus_aureus
 ┣ 📂logs
 ┃ ┣ 📂taxons
 ┃ ┃ ┗ 📜staphylococcus_aureus.log
 ┃ ┗📜lineage.log
 ┣ 📜assembly_finder.log
 ┣ 📜assembly_summary.tsv
 ┣ 📜config.yaml
 ┗ 📜taxonomy.tsv
```

## Download genomes
### Small datasets
* *Staphylococcus aureus* complete genomes

```sh
assembly_finder -i staphylococcus_aureus 
```
!!! note
    By default, assembly_finder searches assembly levels in the following order: **complete**, **chromosome**, **scaffold**, and **contig**.
    
    The search stops at the first assembly level where genomes are found.
    
    This behavior was introduced in [v0.9.0](https://github.com/metagenlab/assembly_finder/releases/tag/v0.9.0) to allow finding the best genomes available for each taxon  

* All *Staphylococcus aureus* genomes 

```sh
assembly_finder -i staphylococcus_aureus --all 
```

!!! note
    The --all option disables the default iteration over assembly levels.
    When used, all genomes for the specified taxon are downloaded, regardless of their assembly level.

* Any *Staphylococcus aureus* complete genome

```sh
assembly_finder -i staphylococcus_aureus -nb 1 
```


### Big datasets

!!! warning

    These examples are for big datasets downloads, so using an NCBI api-key is highly recommended

* Download all chlamydia genomes

```sh
assembly_finder -i chlamydia --all --api-key <api-key>
```

* Best ranking complete genome per bacteria species

```sh
assembly_finder -i eubacteria --api-key <api-key> --rank species --nrank 1
```

* Complete bacteria viruses and archaea genomes from RefSeq <small>(excluding MAGs and atypical)</small>

```sh
assembly_finder -i eubacteria,viruses,archaea \
--api-key <api-key> \
--source refseq \
--mag exclude \
-o outdir
```

* Specific bioproject

```sh
assembly_finder -i PRJNA289059 --api-key <api-key> --accession
```
## Download other files <small>(cds, proteins, gff3 ...)</small>
```sh
assembly_finder -i staphylococcus_aureus --reference \
--include rna,protein,cds,gff3,gtf,gbff,seq-report
```
Output:
```sh
📂staphylococcus_aureus
 ┣ 📂download
 ┃ ┣ 📂GCF_000013425.1
 ┃ ┃ ┣ 📜GCF_000013425.1_ASM1342v1_genomic.fna.gz
 ┃ ┃ ┣ 📜cds_from_genomic.fna.gz
 ┃ ┃ ┣ 📜genomic.gbff.gz
 ┃ ┃ ┣ 📜genomic.gff.gz
 ┃ ┃ ┣ 📜genomic.gtf.gz
 ┃ ┃ ┗ 📜protein.faa.gz
 ┃ ┃ ┗ 📜sequence_report.jsonl 
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
 ┗ 📜taxonomy.tsv
```
