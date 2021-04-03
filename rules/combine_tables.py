import pandas as pd
"""
Main
"""
df_list = []
for file in snakemake.input:
    entry = file.split('/')[2].split('.tsv')[0]
    tb = pd.read_csv(file, sep='\t')
    tb.insert(loc=0, column='UserInputNames', value=[entry]*len(tb))
    df_list.append(tb)
df = pd.concat(df_list,sort=False)
df = df.replace({'Refseq_category': {0: 'reference genome', 1: 'representative genome', 6: 'na'},
                 'AssemblyStatus': {2: 'Complete Genome', 3: 'Chromosome', 4: 'Scaffold', 5: 'Contig', 6: 'na'}})
links = list(df['FtpPath_Genbank'])
filenames = [link.split('/')[-1] for link in links]
df.insert(loc=0, column='AssemblyNames', value=filenames)
df.to_csv(snakemake.output[0], sep='\t', index=None)
