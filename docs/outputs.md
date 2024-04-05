## Directory
This is how the download directory looks like when using the [example](./inputs.md/#taxa-table)
``` sh
taxons
 ┣ download
 ┃ ┣ GCF_000157115.2
 ┃ ┃ ┣ GCF_000157115.2_Escherichia_sp_3_2_53FAA_V2_genomic.fna.gz
 ┃ ┃ ┗ sequence_report.jsonl
 ┃ ┣ GCF_000418345.1
 ┃ ┃ ┣ GCF_000418345.1_ASM41834v1_genomic.fna.gz
 ┃ ┃ ┗ sequence_report.jsonl
 ┃ ┣ GCF_003812505.1
 ┃ ┃ ┣ GCF_003812505.1_ASM381250v1_genomic.fna.gz
 ┃ ┃ ┗ sequence_report.jsonl
 ┃ ┣ .snakemake_timestamp
 ┃ ┣ assembly_data_report.jsonl
 ┃ ┗ dataset_catalog.json
 ┣ archive.zip
 ┣ assembly_summary.tsv
 ┣ sequence_report.tsv
 ┗ taxonomy.tsv
```

## Assembly summary
Table with assembly informations such as assembly level, reference category, checkM and BUSCO completeness, sequencing technology, number of contigs ...
!!! note

    I removed comments, biosample.sample_ids, biosample.attributes and biosample.lineage for visual clarity

| accession       | current_accession | paired_accession | source_database        | annotation_info.method                          | annotation_info.name                               | annotation_info.pipeline                           | annotation_info.provider | annotation_info.release_date | annotation_info.software_version | annotation_info.stats.gene_counts.non_coding | annotation_info.stats.gene_counts.protein_coding | annotation_info.stats.gene_counts.pseudogene | annotation_info.stats.gene_counts.total | assembly_level  | assembly_method                        | assembly_name               | assembly_status | assembly_type | bioproject_accession | biosample.accession | biosample.bioprojects          | biosample.description.organism_name | biosample.description.tax_id | biosample.description.title                                              | biosample.last_updated  | biosample.models | biosample.owner.contacts | biosample.owner.name            | biosample.package | biosample.publication_date | biosample.status.status | biosample.status.when   | biosample.submission_date | paired_assembly.accession | paired_assembly.annotation_name                    | paired_assembly.status | refseq_category       | release_date | sequencing_tech     | submitter                       | contig_l50 | contig_n50 | gc_count | gc_percent | genome_coverage | number_of_component_sequences | number_of_contigs | number_of_scaffolds | scaffold_l50 | scaffold_n50 | total_number_of_chromosomes | total_sequence_length | total_ungapped_length | average_nucleotide_identity.best_ani_match.ani | average_nucleotide_identity.best_ani_match.assembly | average_nucleotide_identity.best_ani_match.assembly_coverage | average_nucleotide_identity.best_ani_match.category | average_nucleotide_identity.best_ani_match.organism_name | average_nucleotide_identity.best_ani_match.type_assembly_coverage | average_nucleotide_identity.category | average_nucleotide_identity.comment | average_nucleotide_identity.match_status | average_nucleotide_identity.submitted_ani_match.ani | average_nucleotide_identity.submitted_ani_match.assembly | average_nucleotide_identity.submitted_ani_match.assembly_coverage | average_nucleotide_identity.submitted_ani_match.category | average_nucleotide_identity.submitted_ani_match.organism_name | average_nucleotide_identity.submitted_ani_match.type_assembly_coverage | average_nucleotide_identity.submitted_organism | average_nucleotide_identity.submitted_species | average_nucleotide_identity.taxonomy_check_status | checkm_info.checkm_marker_set | checkm_info.checkm_marker_set_rank | checkm_info.checkm_species_tax_id | checkm_info.checkm_version | checkm_info.completeness | checkm_info.completeness_percentile | checkm_info.contamination | infraspecific_names.strain | organism_name          | tax_id | biosample.description.comment | common_name | wgs_info.master_wgs_url                             | wgs_info.wgs_contigs_url                       | wgs_info.wgs_project_accession | path                                                                                |
| :-------------- | :---------------- | :--------------- | :--------------------- | :---------------------------------------------- | :------------------------------------------------- | :------------------------------------------------- | :----------------------- | :--------------------------- | :------------------------------- | :------------------------------------------- | :----------------------------------------------- | :------------------------------------------- | :-------------------------------------- | :-------------- | :------------------------------------- | :-------------------------- | :-------------- | :------------ | :------------------- | :------------------ | :----------------------------- | :---------------------------------- | :--------------------------- | :----------------------------------------------------------------------- | :---------------------- | :--------------- | :----------------------- | :------------------------------ | :---------------- | :------------------------- | :---------------------- | :---------------------- | :------------------------ | :------------------------ | :------------------------------------------------- | :--------------------- | :-------------------- | :----------- | :------------------ | :------------------------------ | :--------- | :--------- | :------- | :--------- | :-------------- | :---------------------------- | :---------------- | :------------------ | :----------- | :----------- | :-------------------------- | :-------------------- | :-------------------- | :--------------------------------------------- | :-------------------------------------------------- | :----------------------------------------------------------- | :-------------------------------------------------- | :------------------------------------------------------- | :---------------------------------------------------------------- | :----------------------------------- | :---------------------------------- | :--------------------------------------- | :-------------------------------------------------- | :------------------------------------------------------- | :---------------------------------------------------------------- | :------------------------------------------------------- | :------------------------------------------------------------ | :--------------------------------------------------------------------- | :--------------------------------------------- | :-------------------------------------------- | :------------------------------------------------ | :---------------------------- | :--------------------------------- | :-------------------------------- | :------------------------- | :----------------------- | :---------------------------------- | :------------------------ | :------------------------- | :--------------------- | :----- | :---------------------------- | :---------- | :-------------------------------------------------- | :--------------------------------------------- | :----------------------------- | :---------------------------------------------------------------------------------- |
| GCF_003812505.1 | GCF_003812505.1   | GCA_003812505.1  | SOURCE_DATABASE_REFSEQ | Best-placed reference protein set; GeneMarkS-2+ | NCBI Prokaryotic Genome Annotation Pipeline (PGAP) | NCBI Prokaryotic Genome Annotation Pipeline (PGAP) | NCBI RefSeq              | 2023-03-24                   | 6.4                              | 85                                           | 2141                                             | 34                                           | 2260                                    | Complete Genome | SMRT v. 2.3.0, HGAP v. 3.0             | ASM381250v1                 | current         | haploid       | PRJNA231221          | SAMN10163251        | [{'accession': 'PRJNA231221'}] | Staphylococcus hominis              | 1290                         | Pathogen: clinical or host-associated sample from Staphylococcus hominis | 2019-05-14T13:08:20.304 | ['Pathogen.cl']  | [{}]                     | US Food and Drug Administration | Pathogen.cl.1.0   | 2018-10-02T00:00:00.000    | live                    | 2018-10-02T12:23:11.101 | 2018-10-02T12:23:11.100   | GCA_003812505.1           | NCBI Prokaryotic Genome Annotation Pipeline (PGAP) | current                | representative genome | 2018-11-21   | PacBio; Illumina    | US Food and Drug Administration | 1          | 2220494    | 713682   | 31.5       | 19.6x           | 3                             | 3                 | 3                   | 1            | 2220494      | 3.0                         | 2257431               | 2257431               | 99.99                                          | GCA_900458635.1                                     | 98.99                                                        | type                                                | Staphylococcus hominis                                   | 99.01                                                             | category_na                          | na                                  | species_match                            | 99.99                                               | GCA_900458635.1                                          | 98.99                                                             | type                                                     | Staphylococcus hominis                                        | 99.01                                                                  | Staphylococcus hominis                         | Staphylococcus hominis                        | OK                                                | Staphylococcus hominis        | species                            | 1290                              | v1.2.2                     | 90.97                    | 58.15603                            | 2.63                      | FDAARGOS_575               | Staphylococcus hominis | 1290   | na                            | na          | na                                                  | na                                             | na                             | /path/to/GCF_003812505.1/GCF_003812505.1_ASM381250v1_genomic.fna.gz                 |
| GCF_000418345.1 | GCF_000418345.1   | GCA_000418345.1  | SOURCE_DATABASE_REFSEQ | Best-placed reference protein set; GeneMarkS-2+ | NCBI Prokaryotic Genome Annotation Pipeline (PGAP) | NCBI Prokaryotic Genome Annotation Pipeline (PGAP) | NCBI RefSeq              | 2023-12-14                   | 6.6                              | 80                                           | 2825                                             | 109                                          | 3014                                    | Complete Genome | Newbler v. 2.5.3; Celera software v.7  | ASM41834v1                  | current         | haploid       | PRJNA196937          | SAMN02603524        | [{'accession': 'PRJNA196937'}] | Staphylococcus aureus Bmb9393       | 1321369                      | Sample from Staphylococcus aureus Bmb9393                                | 2015-05-18T13:18:49.507 | ['Generic']      | na                       | NCBI                            | Generic.1.0       | 2014-01-30T14:21:36.850    | live                    | 2014-01-30T14:21:36.850 | 2014-01-30T14:21:36.850   | GCA_000418345.1           | Annotation submitted by LNCC                       | current                | na                    | 2013-07-05   | 454 GS FLX Titanium | LNCC                            | 1          | 2980548    | 981903   | 33.0       | 25.0x           | 2                             | 2                 | 2                   | 1            | 2980548      | 2.0                         | 2983456               | 2983456               | 99.47                                          | GCA_006364675.1                                     | 89.23                                                        | type                                                | Staphylococcus aureus                                    | 95.52                                                             | category_na                          | na                                  | species_match                            | 99.47                                               | GCA_006364675.1                                          | 89.23                                                             | type                                                     | Staphylococcus aureus                                         | 95.52                                                                  | Staphylococcus aureus                          | Staphylococcus aureus                         | OK                                                | Staphylococcus aureus         | species                            | 1280                              | v1.2.2                     | 98.37                    | 55.30531                            | 0.47                      | Bmb9393                    | Staphylococcus aureus  | 1280   | na                            | na          | na                                                  | na                                             | na                             | /path/to/GCF_000418345.1/GCF_000418345.1_ASM41834v1_genomic.fna.gz                  |
| GCF_000157115.2 | GCF_000157115.2   | GCA_000157115.2  | SOURCE_DATABASE_REFSEQ | Best-placed reference protein set; GeneMarkS-2+ | NCBI Prokaryotic Genome Annotation Pipeline (PGAP) | NCBI Prokaryotic Genome Annotation Pipeline (PGAP) | NCBI RefSeq              | 2023-05-18                   | 6.5                              | 98                                           | 4773                                             | 265                                          | 5136                                    | Scaffold        | Newbler v. 2.0.0-PreRelease-04/10/2008 | Escherichia_sp_3_2_53FAA_V2 | current         | haploid       | PRJNA32465           | SAMN02463704        | [{}]                           | Escherichia coli                    | 562                          | MIGS Cultured Bacterial/Archaeal sample from Escherichia coli            | 2019-05-16T02:06:52.537 | ['MIGS.ba']      | na                       | NCBI                            | MIGS.ba.6.0       | 2013-12-18T00:00:00.000    | live                    | 2013-12-18T12:35:26     | 2013-12-18T11:11:28.597   | GCA_000157115.2           | Annotation submitted by Broad Institute            | current                | na                    | 2015-07-15   | 454                 | Broad Institute                 | 8          | 213698     | 2600342  | 50.5       | 19.0x           | 104                           | 104               | 12                  | 2            | 953648       | na                          | 5153453               | 5144253               | 99.77                                          | GCA_000013265.1                                     | 93.44                                                        | claderef                                            | Escherichia coli                                         | 92.8                                                              | category_na                          | na                                  | species_match                            | 99.77                                               | GCA_000013265.1                                          | 93.44                                                             | claderef                                                 | Escherichia coli UTI89                                        | 92.8                                                                   | Escherichia coli                               | Escherichia coli                              | OK                                                | Escherichia coli              | species                            | 562                               | v1.2.2                     | 98.56                    | 27.096497                           | 0.89                      | 3_2_53FAA                  | Escherichia coli       | 562    | Keywords: GSC:MIxS;MIGS:6.0   | E. coli     | https://www.ncbi.nlm.nih.gov/nuccore/ACAC00000000.2 | https://www.ncbi.nlm.nih.gov/Traces/wgs/ACAC02 | ACAC02                         | /path/to/GCF_000157115.2/GCF_000157115.2_Escherichia_sp_3_2_53FAA_V2_genomic.fna.gz |

## Taxonomy

Table containing the full lineage from superkingdom to species of each tax_id

| accession       | current_scientific_name | tax_id | rank    | lineage_id                        | superkingdom | phylum         | class               | order            | family             | genus          | species                |
| :-------------- | :---------------------- | :----- | :------ | :-------------------------------- | :----------- | :------------- | :------------------ | :--------------- | :----------------- | :------------- | :--------------------- |
| GCF_003812505.1 | Staphylococcus hominis  | 1290   | SPECIES | 2,1239,91061,1385,90964,1279,1290 | Bacteria     | Bacillota      | Bacilli             | Bacillales       | Staphylococcaceae  | Staphylococcus | Staphylococcus hominis |
| GCF_000418345.1 | Staphylococcus aureus   | 1280   | SPECIES | 2,1239,91061,1385,90964,1279,1280 | Bacteria     | Bacillota      | Bacilli             | Bacillales       | Staphylococcaceae  | Staphylococcus | Staphylococcus aureus  |
| GCF_000157115.2 | Escherichia coli        | 562    | SPECIES | 2,1224,1236,91347,543,561,562     | Bacteria     | Pseudomonadota | Gammaproteobacteria | Enterobacterales | Enterobacteriaceae | Escherichia    | Escherichia coli       |

## Sequence summary

Table with sequence summary information like GenBank or RefSeq accession, molecule type (Chromosome, Plasmid), sequence length ...

| Assembly Accession | Assembly Unplaced Count | Assembly-unit accession | Chromosome name | GC Count | GC Percent | GenBank seq accession | Molecule type | Ordering | RefSeq seq accession | Role               | Seq length | UCSC style name | Unlocalized Count |
| :----------------- | :---------------------- | :---------------------- | :-------------- | :------- | :--------- | :-------------------- | :------------ | :------- | :------------------- | :----------------- | :--------- | :-------------- | :---------------- |
| GCF_000157115.2    |                         | Primary Assembly        | chromosome      |          |            | KQ235739.1            | Chromosome    |          | NZ_KQ235739.1        | unplaced-scaffold  | 953648     |                 |                   |
| GCF_000157115.2    |                         | Primary Assembly        | chromosome      |          |            | KQ235740.1            | Chromosome    |          | NZ_KQ235740.1        | unplaced-scaffold  | 35389      |                 |                   |
| GCF_000157115.2    |                         | Primary Assembly        | chromosome      |          |            | KQ235741.1            | Chromosome    |          | NZ_KQ235741.1        | unplaced-scaffold  | 782863     |                 |                   |
| GCF_000157115.2    |                         | Primary Assembly        | chromosome      |          |            | KQ235742.1            | Chromosome    |          | NZ_KQ235742.1        | unplaced-scaffold  | 4343       |                 |                   |
| GCF_000157115.2    |                         | Primary Assembly        | chromosome      |          |            | KQ235743.1            | Chromosome    |          | NZ_KQ235743.1        | unplaced-scaffold  | 5458       |                 |                   |
| GCF_000157115.2    |                         | Primary Assembly        | chromosome      |          |            | KQ235744.1            | Chromosome    |          | NZ_KQ235744.1        | unplaced-scaffold  | 176042     |                 |                   |
| GCF_000157115.2    |                         | Primary Assembly        | chromosome      |          |            | KQ235745.1            | Chromosome    |          | NZ_KQ235745.1        | unplaced-scaffold  | 529000     |                 |                   |
| GCF_000157115.2    |                         | Primary Assembly        | chromosome      |          |            | KQ235746.1            | Chromosome    |          | NZ_KQ235746.1        | unplaced-scaffold  | 2543925    |                 |                   |
| GCF_000157115.2    |                         | Primary Assembly        | chromosome      |          |            | KQ235747.1            | Chromosome    |          | NZ_KQ235747.1        | unplaced-scaffold  | 1815       |                 |                   |
| GCF_000157115.2    |                         | Primary Assembly        | chromosome      |          |            | KQ235748.1            | Chromosome    |          | NZ_KQ235748.1        | unplaced-scaffold  | 1451       |                 |                   |
| GCF_000157115.2    |                         | Primary Assembly        | chromosome      |          |            | KQ235749.1            | Chromosome    |          | NZ_KQ235749.1        | unplaced-scaffold  | 113382     |                 |                   |
| GCF_000157115.2    |                         | Primary Assembly        | chromosome      |          |            | KQ235750.1            | Chromosome    |          | NZ_KQ235750.1        | unplaced-scaffold  | 6137       |                 |                   |
| GCF_000418345.1    |                         | Primary Assembly        | chromosome      | 981048   |            | CP005288.1            | Chromosome    |          | NC_021670.1          | assembled-molecule | 2980548    |                 |                   |
| GCF_000418345.1    |                         | Primary Assembly        | pBmb9393        | 855      |            | CP005289.1            | Plasmid       |          | NC_021657.1          | assembled-molecule | 2908       |                 |                   |
| GCF_003812505.1    |                         | Primary Assembly        | chromosome      | 702792   |            | CP033732.1            | Chromosome    |          | NZ_CP033732.1        | assembled-molecule | 2220494    |                 |                   |
| GCF_003812505.1    |                         | Primary Assembly        | unnamed1        | 9555     |            | CP033731.1            | Plasmid       |          | NZ_CP033731.1        | assembled-molecule | 32498      |                 |                   |
| GCF_003812505.1    |                         | Primary Assembly        | unnamed2        | 1335     |            | CP033733.1            | Plasmid       |          | NZ_CP033733.1        | assembled-molecule | 4439       |                 |                   |		