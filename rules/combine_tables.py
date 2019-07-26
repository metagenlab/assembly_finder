'''
Main
'''

import pandas as pd
df_list=[pd.read_csv(i,delimiter='\t') for i in snakemake.input]
df=pd.concat(df_list,sort=False)
df=df.replace([0, 1, 2, 3, 4, 5, 6],['reference genome', 'representative genome','Complete Genome',
                                     'Chromosome','Scaffold','Contig','na'])
df.to_csv('table_combined.tsv',sep='\t',index=False)