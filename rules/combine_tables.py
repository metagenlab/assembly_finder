'''
Main
'''

import pandas as pd
df_list=[pd.read_csv(i,delimiter='\t') for i in snakemake.input]
df=pd.concat(df_list,sort=False)
df=df.replace({'Refseq_category':{0:'reference genome',1:'representative genome',6:'na'},'AssemblyStatus':{2:'Complete Genome',3:
                                     'Chromosome',4:'Scaffold',5:'Contig',6:'na'}})
links=list(df['FtpPath_Genbank'])
filenames=[]
for i in links:
    splt=i.split('/')
    filenames.append(splt[len(splt)-1])
df.insert(loc=0,column='AssemblyNames',value=filenames)
df.to_csv('table_combined.tsv',sep='\t',index=False)