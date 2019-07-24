
from Bio import Entrez
#import random as rd
import pandas as pd


def taxid_find(name_input):
    '''
    Function to return a taxonomy id based on a query
    '''

    print('\n> Searching for taxIDs', name_input, '...')
    try:
        int(name_input)
        print('Query is a taxID')
        taxid = name_input

    except ValueError:
        print('Query is not a taxID\nSearching for TaxID')
        taxid = Entrez.read(Entrez.esearch(db='taxonomy', term=name_input, retmax=50))['IdList']
        if len(taxid) == 1:
            print('One TaxID:%s found' % taxid[0])
            # TaxID.append(taxid[0])
        if len(taxid) > 1:
            print('%sTaxIDfound' % len(taxid))
        if len(taxid) == 0:
            print('\nERROR: TaxID not found! \nChange search term!')
    return taxid[0]


def search_assemblies(TaxID, Genbank=False, Refseq=True,
                      representative=True, reference=True, complete=True, exclude_metagenomes=True):
    search_term = 'txid%s[Organism:exp] AND (latest[filter]) AND (all[filter] NOT "derived from surveillance project"[filter] AND all[filter] NOT anomalous[filter])) '

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

    print('> Search term \n%s' % search_term)

    assembly_dic = Entrez.read(Entrez.esearch(db='assembly', term=search_term % TaxID, retmax=100000))['IdList']

    print('\nFound %i assemblies' % len(assembly_dic))
    print(assembly_dic)

    if len(assembly_dic) == 0:
        assembly_dic = None
        print('\nError: change search term')

    return assembly_dic

def generate_table(assembly_list):
    assembly_tax_link = Entrez.read(Entrez.elink(dbfrom='assembly', db='taxonomy', id=assembly_list))
    all_taxid = [i['Id'] for i in assembly_tax_link[0]['LinkSetDb'][0]['Link']]
    taxonomy_summary = Entrez.read(Entrez.efetch(db='taxonomy', id=all_taxid))
    assembly_summary = Entrez.read(Entrez.esummary(db='assembly', id=assembly_list))
    target_ranks = ['superkingdom', 'phylum', 'order', 'family', 'genus', 'species']

    dic = {}#dic = dictionary regrouping info on tax ids
    for i in range(len(taxonomy_summary)):
        dic[all_taxid[i]] = {}
        if taxonomy_summary[i]['Rank'] in target_ranks:
            dic[all_taxid[i]][taxonomy_summary[i]['Rank']] = taxonomy_summary[i]['ScientificName']
        for j in taxonomy_summary[i]['LineageEx']:

            rank = j['Rank']
            if rank in target_ranks:
                names = j['ScientificName']
                dic[all_taxid[i]][rank] = names

    dico = {}# dico = dictionary regrouping info on the assemblies
    for i in range(len(assembly_summary['DocumentSummarySet']['DocumentSummary'])):
        assembly_ID = assembly_summary['DocumentSummarySet']['DocumentSummary'][i].attributes['uid']
        taxonomy_ID = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['Taxid']
        dico[assembly_ID] = dic[taxonomy_ID]# get the tax info for all assemblies
        Refseq_cat = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['RefSeq_category']
        assembly_stat = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['AssemblyStatus']
        ftp_path = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['FtpPath_RefSeq']
        dico[assembly_ID]['Refseq_cat'] = Refseq_cat
        dico[assembly_ID]['AssemblySatus'] = assembly_stat
        dico[assembly_ID]['FtpPath_RefSeq'] = ftp_path

    assembly_table = pd.DataFrame.from_dict(dico, orient='index')
    return assembly_table

'''
Main
'''
Entrez.email = snakemake.params['NCBI_email']



Entrez.api_key = snakemake.params['NCBI_key']

complete_or_not=snakemake.params['complete_assemblies']

print(complete_or_not)

print(snakemake.wildcards.entry)

taxid=taxid_find(snakemake.wildcards.entry)


print(taxid)
assemblies=search_assemblies(taxid,complete=complete_or_not)
tb=generate_table(assemblies)
print(snakemake.output[0])
tb.to_csv(snakemake.output[0],sep='\t')


