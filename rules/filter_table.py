import logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
logger.addHandler(logging.FileHandler('assembly_finder.log', 'a'))
print = logger.info

import pandas as pd
def select_assemblies(table, nb=5, rank_to_select='None'):

    fact_table = table.replace(['reference genome', 'representative genome', 'Complete Genome', 'Chromosome',
                                    'Scaffold','Contig','na'], [0, 1, 2, 3, 4, 5, 6])
    sorted_table = fact_table.sort_values(['Refseq_category', 'AssemblyStatus'], ascending=[True, True])

    if rank_to_select != 'None':
        print('Filtering according to {0}, sorting {1} and {2}'.format(rank_to_select,'assembly status',
                                                                       'Refseq category'))
        select_index = []
        unique_list = list(set(sorted_table[rank_to_select]))

        if len(unique_list) > 1:
            for i in unique_list:
                select_index.append(sorted_table[sorted_table[rank_to_select] == i].sample(1).index[0])
                #randomly select one assembly ID for each unique selected rank (species for example)
            sorted_table=sorted_table.loc[select_index, :]
        if len(unique_list)==1:
            print('Same {0} for all assemblies, no filtering'.format(rank_to_select))
        if len(unique_list)==0:
            print('{0} is not a target rank'.format(rank_to_select))
    else :
        print('No filter specified, sorting assembly status and Refseq category')

    sorted_table = sorted_table[0:nb]
    return sorted_table

'''
Main
'''
input_tb=pd.read_csv(snakemake.params['tb_path'],
                          delimiter='\t',index_col=0)
input_tb.index=input_tb.index.astype('str')
nb_assemblies=input_tb.loc[snakemake.wildcards.entry]['nb_genomes']

tb=pd.read_csv(snakemake.input[0],delimiter='\t')
filtered_table=select_assemblies(tb,nb=nb_assemblies,rank_to_select=snakemake.params['rank_filter'])
filtered_table.to_csv(snakemake.output[0],sep='\t',index=False)