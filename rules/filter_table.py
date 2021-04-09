import logging
import pandas as pd
logging.basicConfig(format='%(asctime)s %(message)s', filename=snakemake.log[0], level=logging.DEBUG)


def select_assemblies(table, nb, rank_to_select='None'):

    fact_table = table.replace({'Refseq_category': {'reference genome': 0, 'representative genome': 1, 'na': 6},
                                'AssemblyStatus':{'Complete Genome': 2, 'Chromosome': 3,
                                    'Scaffold': 4, 'Contig': 5, 'na': 6}})
    sorted_table = fact_table.sort_values(['Refseq_category', 'AssemblyStatus', 'Contig_count', 'Release_date_Genbank'],
                                          ascending=[True, True, True, False])

    if rank_to_select != 'None':
        logging.info(f'Filtering according to {rank_to_select}, Refseq categories, assembly status, '
                     f'contig count and release date')
        select_index = []
        unique_list = list(set(sorted_table[rank_to_select]))
        if len(unique_list) > 1:
            for i in unique_list:
                select_index.append(sorted_table[sorted_table[rank_to_select] == i].sample(1).index[0])
                # randomly select one assembly ID for each unique selected rank (species for example)
            sorted_table=sorted_table.loc[select_index, :]
        if len(unique_list) == 1:
            logging.info(f'Same {rank_to_select} for all assemblies, no filtering')
        if len(unique_list) == 0:
            logging.error(f'{rank_to_select} is not a target rank')
    else :
        logging.info('No taxonomic rank specified, sorting according to Refseq category, '
                     'assembly status, contig count and release date')

    if len(sorted_table) >= nb:
        logging.info(f'Selecting {nb} sorted assemblies out of {len(sorted_table)}')
        sorted_table = sorted_table[0:nb]
    if len(sorted_table) < nb:
        logging.warning(f'Found less than {nb} assemblies in total, returning {len(sorted_table)} instead')
    return sorted_table


'''
Main
'''
input_tb = pd.read_csv(snakemake.params['tb_path'], sep='\t', index_col=0)
input_tb.index = input_tb.index.astype('str')
nb_assemblies = int(input_tb.loc[snakemake.wildcards.entry]['nb_genomes'])

tb=pd.read_csv(snakemake.input[0], sep='\t')
filtered_table = select_assemblies(tb, nb=nb_assemblies,rank_to_select=snakemake.params['rank_filter'])
filtered_table.to_csv(snakemake.output[0], sep='\t', index=False)