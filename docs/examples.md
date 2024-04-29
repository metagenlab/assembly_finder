## Small datasets

### *Staphylococcus aureus* reference genome

!!! note
    By default, assembly_finder limits genomes to reference or representative

```sh
assembly_finder -i staphylococcus_aureus -nb 1
```
### Any *Staphylococcus aureus* genome
```sh
assembly_finder -i staphylococcus_aureus -nb 1 --reference False
```
### Download from a list of taxa

```sh
assembly_finder -i 1290,1813735,114185 -o taxa -nb 1
```
### Download from a table of taxa

??? info "**taxa.tsv**"
        

    | taxon | nb |
    | :-------------------- | :-- |
    | 1290 | 1 |
    | 1813735 | 1 |
    | 114185 | 1 |

```sh
assembly_finder -i taxa.tsv
```

## Big datasets

!!! warning

    These examples are for big datasets downloads, so using an NCBI api-key is highly recommended

### Download all chlamydia genomes

```sh
assembly_finder -i chlamydia --api-key <api-key>
```

### Best ranking genome for each bacteria species

```sh
assembly_finder -i bacteria --api-key <api-key> --rank species --nrank 1
```

### Complete RefSeq bacteria viruses and archaea <small>(excluding MAGs and atypical)</small>

```sh
assembly_finder -i bacteria,viruses,archaea -o outdir --api-key <api-key> --source refseq --assembly-level complete --mag exclude --atypical
```

### Specific bioproject

```sh
assembly_finder -i PRJNA289059 --api-key <api-key> --accession
```
