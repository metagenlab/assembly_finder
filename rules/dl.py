import os
import ftplib
import pandas as pd

def dl_fna_ftp(ftp_login,input_table):
        table=pd.read_csv(input_table,delimiter='\t',index_col=0)
        full_path=table.loc[snakemake.wildcards.sample]['FtpPath_Genbank']
        split_path=full_path.split('/')
        #print(split_path)
        redundant_name=split_path[-1]
        print(redundant_name+'_genomic.fna.gz')
        path=full_path.split('ftp://ftp.ncbi.nlm.nih.gov/')[1]
        print(path)
        ftp_login.cwd(path)
        #filename=snakemake.wildcards.sample
        print('Downloading from this path: %s'%path)
        print(snakemake.output[0])
        print(os.path.join(os.getcwd(), snakemake.output[0]))
        ftp_login.retrbinary("RETR "+redundant_name+'_genomic.fna.gz',
                             open(os.path.join(os.getcwd(),snakemake.output[0]),"wb").write)
        print('Done')

'''
Main
'''
ftp=ftplib.FTP(host='ftp.ncbi.nih.gov',user='anonymous',passwd=snakemake.params['NCBI_email'])
print('Downloading assemblies')

dl_fna_ftp(ftp,snakemake.input[0])
