import pandas as pd
def parse_summary_tb(wildcards):
    checkpoint_output = checkpoints.combine_assembly_tables.get(**wildcards).output[0]
    tb=pd.read_csv(checkpoint_output,delimiter='\t')
    filenames=tb['AssemblyNames']
    expd=expand('assembly_gz/{assemblyname}.fna.gz',assemblyname=filenames)
    return expd


rule all_download:
    input: parse_summary_tb

include: 'rules/find_assemblies.rules'

