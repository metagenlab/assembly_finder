!!! warning
    These examples are for big datasets downloads, so using an NCBI api-key is highly recommended

## Best ranking genome for each bacteria species

```sh
assembly_finder -i bacteria --api-key <api-key> --rank species --nrank 1
```

## Complete RefSeq bacteria viruses and archaea (excluding metagenomes and atypical genomes)

```sh
assembly_finder -i bacteria,viruses,archaea -o outdir --api-key <api-key> --source refseq --assembly-level complete --mag exclude --atypical
```

## Download sequences from a specific bioproject 

```sh
assembly_finder -i PRJNA289059 --accession
```
