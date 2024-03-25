## Download the top 1 genome per bacterial species

```sh
assembly_finder -i bacteria --api-key <api-key> --rank species --nrank 1
```

## Donwload all refseq bacteria viruses and archaea complete genomes (exclude metagenome assembled and atypical genomes)

```sh
assembly_finder -i bacteria,viruses,archaea -o outdir --api-key <api-key> --source refseq --assembly-level complete --mag exclude --atypical
```

## Download sequences from a specific bioproject 

```sh
assembly_finder -i PRJNA289059 --accession
```
