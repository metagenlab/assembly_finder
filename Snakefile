rule all:
    input: 'table_combined.tsv'

include: 'rules/find_assemblies.rules'

