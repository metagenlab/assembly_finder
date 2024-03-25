# Inputs
Input can be either a string or a file, and queries can be either taxa or accession as shown in [ncbi-datasets-cli](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install).
## String
### Taxa
``` sh
assembly_finder -i 1290,staphylococcus_aureus,562 -nb 1 -o taxons
```
### Accessions
``` sh
assembly_finder -i GCF_003812505.1,GCF_000418345.1,GCF_000157115.2 --accession -o accessions
```

## File
### Taxa table
| taxon                 | nb  |
| :-------------------- | :-- |
| 1290                  | 1   |
| staphylococcus_aureus | 1   |
| 562                   | 1   |

### Accessions list
|  GCF_003812505.1  |
| :----------------|
| GCF_000418345.1 |
| GCF_000157115.2 |