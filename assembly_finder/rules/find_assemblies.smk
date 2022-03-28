import pandas as pd
community_name = config['community_name']
try:
    nrank = config['n_by_rank']
    ncbi_key = config['NCBI_key'] 
    ncbi_email = config['NCBI_email']
    alvl = config['assembly_level']
    db = config['db']
    rcat = config['refseq_category']
    excl = config['exclude']
    annot = config['annotation']
    rank = config['Rank_to_filter_by']
except KeyError:
    ncbi_key = ''
    ncbi_email = ''
    alvl = ['complete genome']
    db = 'refseq'
    rcat = ['reference', 'representative']
    excl = ['metagenome']
    annot = True
    rank = None
    nrank = 1
try:
    entries = list(pd.read_csv(config['input_table_path'],delimiter='\t')['Taxonomy'])
    isassembly = False
    col = 'Taxonomy'
except KeyError:
    entries = list(pd.read_csv(config['input_table_path'],delimiter='\t')['Assembly'])
    isassembly = True
    col = 'Assembly'

rule check_for_update_ete3:
    container: "docker://metagenlab/assemblyfinder:v.1.1"

    output: temp('ete3-update.txt')

    log: f'{community_name}/logs/ete3-update.log'

    script: 'update-ete3.py'

rule get_assembly_tables:
    container: "docker://metagenlab/assemblyfinder:v.1.1"

    input: config["input_table_path"],
           'ete3-update.txt'

    output: all=f'{community_name}/tables/{{entry}}-all.tsv',
            filtered=f'{community_name}/tables/{{entry}}-filtered.tsv'

    params: ncbi_key=ncbi_key, ncbi_email=ncbi_email,
            alvl=alvl, db=db, rcat=rcat, excl=excl,
            annot = annot, rank=rank, n_by_rank=nrank,
            assembly=isassembly, column=col

    resources: ncbi_requests=1

    log: f'{community_name}/logs/find-assemblies-{{entry}}.log'

    benchmark: f"{community_name}/benchmark/find-assemblies-{{entry}}.txt"

    script: 'assembly_table.py'


rule combine_assembly_tables:
    container: "docker://metagenlab/assemblyfinder:v.1.1"

    input: expand(f'{community_name}/tables/{{entry}}-filtered.tsv',entry=entries)

    output: f'{community_name}/assemblies-summary.tsv'

    params: column=col

    script: 'combine_tables.py'


rule get_ftp_links_list:
    container: "docker://metagenlab/assemblyfinder:v.1.1"

    input: f'{community_name}/assemblies-summary.tsv'

    output: temp(f"{community_name}/assemblies/ftp-links.txt")

    params: db=db
    
    script: "concat-ftp-links.py"


checkpoint download_assemblies:
    container: "docker://metagenlab/aspera-cli-conda:v.1.0"

    input: f"{community_name}/assemblies/ftp-links.txt"

    output: f"{community_name}/assemblies/dl.done"

    log: f"{community_name}/logs/download.log"

    benchmark: f"{community_name}/benchmark/downloads.txt"

    params: f"{community_name}/assemblies"

    shell:
        """
        ascp -T -k 1 -i ${{CONDA_PREFIX}}/etc/asperaweb_id_dsa.openssh --mode=recv --user=anonftp \
        --host=ftp.ncbi.nlm.nih.gov --file-list={input} {params} 1> {log}
        touch {output}
        """