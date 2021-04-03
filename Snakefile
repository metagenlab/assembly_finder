import pandas as pd
include: 'rules/find_assemblies.rules'

def downloaded_list(wildcards):
    f = checkpoints.combine_assembly_tables.get(**wildcards).output[0]
    assemblynames = pd.read_csv(f,sep='\t')['AssemblyNames']
    return expand("assembly_gz/{assemblyname}_genomic.fna.gz", assemblyname=assemblynames)

community_name=config["community_name"]
rule all_download:
    input: f"{community_name}-assemblies-summary.tsv",
           downloaded_list