import logging
import os
import ftplib
import pandas as pd
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', datefmt='%d %b %Y %H:%M:%S',
                    filename=snakemake.log[0], level=logging.DEBUG)


def dl_fna_ftp(ftp_login, input_table):
        table = pd.read_csv(input_table,sep='\t',index_col=0)
        full_path = table.loc[snakemake.wildcards.assemblyname]['FtpPath_Genbank']
        split_path = full_path.split('/')
        redundant_name = split_path[-1]
        path = full_path.split('ftp://ftp.ncbi.nlm.nih.gov/')[1]
        ftp_login.cwd(path)
        logging.info(f'Downloading from : {path}')
        logging.info(f'Downloading to {snakemake.output[0]}')
        ftp_login.retrbinary("RETR " + redundant_name + '_genomic.fna.gz',
                             open(os.path.join(os.getcwd(), snakemake.output[0]), "wb").write)
        logging.info('Done')


"""
Main
"""
ftp = ftplib.FTP(host='ftp.ncbi.nih.gov', user='anonymous', passwd=snakemake.params['NCBI_email'])
dl_fna_ftp(ftp, snakemake.input[0])