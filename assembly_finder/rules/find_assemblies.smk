import pandas as pd

outdir = config['outdir']
ncbi_key = config['NCBI_key'] 
ncbi_email = config['NCBI_email']
uid = config['uid']
alvl = config['assembly_level']
db = config['db']
rcat = config['refseq_category']
excl = config['exclude']
annot = config['annotation']
rank = config['Rank_to_filter_by']
nrank = config['n_by_rank']
nb = config['n_by_entry']

try:
    entries = list(pd.read_csv(config['input'], delimiter='\t')[0])

except FileNotFoundError:
    entries = config['input'].split(',')

rule check_for_update_ete3:
    output: temp('ete3-update.txt')

    log: f'{outdir}/logs/ete3-update.log'

    script: 'update-ete3.py'

rule get_assembly_tables:
    input: 'ete3-update.txt'

    output: all=f'{outdir}/tables/{{entry}}-all.tsv',
            filtered=f'{outdir}/tables/{{entry}}-filtered.tsv'

    params: ncbi_key=ncbi_key,ncbi_email=ncbi_email,
            alvl=alvl,db=db,uid=uid,rcat=rcat,excl=excl,
            annot=annot,rank=rank,n_by_rank=nrank,nb=nb

    resources: ncbi_requests=1

    log: f'{outdir}/logs/find-assemblies-{{entry}}.log'

    benchmark: f"{outdir}/benchmark/find-assemblies-{{entry}}.txt"

    script: 'assembly_table.py'


rule combine_assembly_tables:
    input: expand(f'{outdir}/tables/{{entry}}-filtered.tsv',entry=entries)

    output: f'{outdir}/assemblies-summary.tsv'

    run:
        pd.concat([pd.read_csv(tb, sep='\t') for tb in list(input)]).to_csv(output[0], sep='\t', index=None)
        


rule get_ftp_links_list:
    input: f'{outdir}/assemblies-summary.tsv'

    output: temp(f"{outdir}/assemblies/ftp-links.txt")

    params: db=db
    
    script: "concat-ftp-links.py"


checkpoint download_assemblies:
    input: f"{outdir}/assemblies/ftp-links.txt"

    output: f"{outdir}/assemblies/dl.done"

    log: f"{outdir}/logs/download.log"

    benchmark: f"{outdir}/benchmark/downloads.txt"

    params: f"{outdir}/assemblies"

    shell:
        """
        ascp -T -k 1 -i ${{CONDA_PREFIX}}/etc/asperaweb_id_dsa.openssh --mode=recv --user=anonftp \
        --host=ftp.ncbi.nlm.nih.gov --file-list={input} {params} 1> {log}
        touch {output}
        """