
from Bio import Entrez

import pandas as pd

import logging

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',datefmt='%d %b %Y %H:%M:%S',filename=snakemake.log[0], level=logging.DEBUG)

def taxid_find(name_input):
    '''
    Function to parse input_list to convert scientific names to taxid
    returns dictionary with taxid found",
    '''

    logging.info('> Searching for taxIDs {0} ...'.format(name_input))
    try:
        int(name_input)
        logging.info('Query is a taxID')
        taxid = name_input

    except ValueError:
        logging.warning('Query is not a taxID, enter taxID to be more precise')
        logging.info('Search term: {0}[all Names]'.format(name_input))
        taxid_list = Entrez.read(Entrez.esearch(db='taxonomy', term='{0}[all Names]'.format(name_input), retmax=100))['IdList']
        if len(taxid_list) == 1:
            taxid = taxid_list[0]
            logging.info('TaxID:{0} found'.format(taxid))
        if len(taxid_list) > 1:
            taxid = taxid_list[0]
            logging.warning('{0} TaxIDfound, change query (taking first one : {1})'.format(len(taxid_list)), taxid)
        if len(taxid_list) == 0:
            taxid=None
            logging.error('TaxID not found! Change search term!')

    return taxid

def search_assemblies(TaxID, Genbank=False, Refseq=True,
                      representative=True, reference=True, complete=True, exclude_metagenomes=True):
    search_term = 'txid{0}[Organism:exp] AND (latest[filter]) AND (all[filter] NOT "derived from surveillance project"[filter] AND all[filter] NOT anomalous[filter])) '

    if Refseq and not Genbank:
        search_term += 'AND("latest refseq"[filter]) '
    if Genbank and not Refseq:
        search_term += 'AND("latest genbank"[filter]) '
    if Genbank and Refseq:
        search_term += 'AND ("latest genbank"[filter] OR "latest refseq"[filter]) '
    if complete:
        search_term += 'AND ("complete genome"[filter]) '
    if representative and reference:
        search_term += 'AND ("representative genome"[filter] OR "reference genome"[filter]) '
    if representative and not reference:
        search_term += 'AND ("representative genome"[filter]) '
    if reference and not representative:
        search_term += 'AND ("reference genome"[filter]) '
    if exclude_metagenomes:
        search_term += 'AND (all[filter] NOT "derived from metagenome"[filter])'

    logging.info('> Search term: {0}'.format(search_term.format(TaxID)))

    assembly_dic = Entrez.read(Entrez.esearch(db='assembly', term=search_term.format(TaxID), retmax=100000))['IdList']

    logging.info('Found {0} assemblies'.format(len(assembly_dic)))

    if len(assembly_dic) == 0:
        assembly_dic = None
        logging.error('No assemblies found ! Change search term!')

    return assembly_dic

def _chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def generate_table(assembly_list):
    assembly_list=','.join(assembly_list)
    assembly_tax_link = Entrez.read(Entrez.elink(dbfrom='assembly', db='taxonomy', id=assembly_list))
    all_taxid = [i['Id'] for i in assembly_tax_link[0]['LinkSetDb'][0]['Link']]
    taxonomy_summary = Entrez.read(Entrez.efetch(db='taxonomy', id=','.join(all_taxid)))
    assembly_summary = Entrez.read(Entrez.esummary(db='assembly', id=assembly_list))

    target_ranks = ['species', 'genus', 'family', 'order', 'phylum', 'superkingdom']
    dic = {}
    for i in range(len(taxonomy_summary)):
        dic[all_taxid[i]] = {}
        if taxonomy_summary[i]['Rank'] in target_ranks:
            dic[all_taxid[i]][taxonomy_summary[i]['Rank']] = taxonomy_summary[i]['ScientificName']
        for j in taxonomy_summary[i]['LineageEx']:
            rank = j['Rank']
            if rank in target_ranks:
                names = j['ScientificName']
                dic[all_taxid[i]][rank] = names

        for n, rank in enumerate(target_ranks):
            if rank not in dic[all_taxid[i]]:  # If a target rank is lacking
                previous_rank = dic[all_taxid[i]][target_ranks[n - 1]].split(" ")[
                    0]  # Get previous rank from the position of the lacking rank
                placeholder = '%s_' % rank[0] + previous_rank  # First letter of the lacking rank + previous rank
                dic[all_taxid[i]][rank] = placeholder

    dico = {}
    assembly_ID = [assembly_summary['DocumentSummarySet']['DocumentSummary'][i].attributes['uid'] for i in
                   range(len(assembly_summary['DocumentSummarySet']['DocumentSummary']))]

    for i in range(len(assembly_ID)):
        dico[assembly_ID[i]] = {}
        taxonomy_ID = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['Taxid']
        Refseq_cat = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['RefSeq_category']
        Genbank_path = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['FtpPath_GenBank']
        assembly_stat = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['AssemblyStatus']
        Refseq_path = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['FtpPath_RefSeq']
        dico[assembly_ID[i]]['AssemblyStatus'] = assembly_stat
        dico[assembly_ID[i]]['Refseq_category'] = Refseq_cat
        dico[assembly_ID[i]]['FtpPath_Refseq'] = Refseq_path
        dico[assembly_ID[i]]['FtpPath_Genbank'] = Genbank_path
        for rank in dic[taxonomy_ID]:
            dico[assembly_ID[i]][rank] = dic[taxonomy_ID][rank]

    assembly_table = pd.DataFrame.from_dict(dico, orient='index')
    return assembly_table

'''
Main
'''
Entrez.email = snakemake.params['NCBI_email']



Entrez.api_key = snakemake.params['NCBI_key']

comp=snakemake.params['comp']

ref=snakemake.params['ref']

rep=snakemake.params['rep']

met=snakemake.params['met']

gb=snakemake.params['gb']

rs=snakemake.params['rs']


taxid=taxid_find(snakemake.wildcards.entry)


assemblies_found=search_assemblies(taxid,Genbank=gb,Refseq=rs,complete=comp,
                                   reference=ref,representative=rep,exclude_metagenomes=met)

if len(assemblies_found)>300:
    assemblies_chunks=_chunks(assemblies_found,300)#Divide assembly lists by chunks of 300 if more than 300 found
    logging.info('Parsing assemblies by chucks of 300')
    table_chunks=[]
    for n,chunks in enumerate(assemblies_chunks):
        logging.info('chunk n°{0}'.format(n))
        tb=generate_table(chunks)
        table_chunks.append(tb)
        #time.sleep(1)
    all_tb = pd.concat(table_chunks,sort=False)
else:
    all_tb=generate_table(assemblies_found)

all_tb.index.name='AssemblyID'

all_tb.to_csv(snakemake.output[0],sep='\t')


