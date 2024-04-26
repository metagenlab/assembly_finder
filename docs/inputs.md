# Inputs

Input can be either a string or a table, and queries can be either taxa or accession as shown in [NCBI datasets docs](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install).

## Strings

=== "Taxons"

    !!! note

        Taxons can be either taxids or taxon names

    ```sh
    assembly_finder -i 1290,staphylococcus_aureus,562 -nb 1 -o taxons
    ```

=== "Accessions"

    ```sh
    assembly_finder --accession -i GCF_003812505.1,GCF_000418345.1,GCF_000157115.2 -o accessions
    ```

## Tables

=== "Taxons"

    !!! note
        You can set the number of genomes per taxon in the table

    | taxon | nb |
    | :-------------------- | :-- |
    | 1290 | 1 |
    | staphylococcus_aureus | 1 |
    | 562 | 1 |

=== "Accessions"

    !!! note

        The accession table does not have a header

    | GCF_003812505.1 |
    | :-------------- |
    | GCF_000418345.1 |
    | GCF_000157115.2 |
