from Bio import Entrez
import pandas as pd
import logging
from ete3 import NCBITaxa
ncbi = NCBITaxa()
if snakemake.params.update:
    ncbi.update_taxonomy_database()
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',datefmt='%d %b %Y %H:%M:%S',
                    filename=snakemake.log[0], level=logging.DEBUG)


def taxid_find(name_input):
    """
    Function to parse input_list to convert scientific names to taxid
    returns dictionary with taxid found
    """
    logging.info(f'> Searching for taxIDs {name_input} ...')
    try:
        int(name_input)
        logging.info('Query is a taxID')
        taxid = name_input

    except ValueError:
        logging.warning('Query is not a taxID, enter taxID to be more precise')
        logging.info(f'Search term: {name_input}[all Names]')
        taxid_list = Entrez.read(Entrez.esearch(db='taxonomy', term=f'{name_input}[all Names]', retmax=100))['IdList']
        if len(taxid_list) == 1:
            taxid = taxid_list[0]
            logging.info(f'TaxID:{taxid} found')
        if len(taxid_list) > 1:
            taxid = taxid_list[0]
            logging.warning(f'{len(taxid_list)} TaxIDs found, change query (taking first one : {taxid})')
        if len(taxid_list) == 0:
            taxid = None
            logging.error('TaxID not found! Change search term!')
    return taxid


def search_assemblies(TaxID, Genbank=False, Refseq=True, representative=True, reference=True, complete=True,
                      exclude_metagenomes=True):
    search_term = f'txid{TaxID}[Organism:exp] '

    if Refseq and not Genbank:
        search_term += 'AND("latest refseq"[filter]) '
    if Genbank and not Refseq:
        search_term += 'AND("latest genbank"[filter]) '
    if Genbank and Refseq:
        search_term += 'AND ("latest genbank"[filter] OR "latest refseq"[filter]) '
    if complete and not representative and not reference:
        search_term += 'AND ("complete genome"[filter]) '
    if complete and representative and not reference:
        search_term += 'AND ("complete genome"[filter]) OR ("representative genome"[filter]) '
    if complete and representative and reference:
        search_term += 'AND ("complete genome"[filter]) OR ("representative genome"[filter]) OR ' \
                       '("reference genome"[filter]) '
    if representative and not reference:
        search_term += 'AND ("representative genome"[filter]) '
    if reference and not representative:
        search_term += 'AND ("reference genome"[filter]) '
    if representative and reference:
        search_term += 'AND ("representative genome"[filter] OR "reference genome"[filter]) '
    if exclude_metagenomes:
        search_term += 'AND (all[filter] NOT "derived from metagenome"[filter])'

    assembly_dic = Entrez.read(Entrez.esearch(db='assembly', term=search_term, retmax=100000))['IdList']
    logging.info(f'> Search term: {search_term}')
    logging.info(f'found {len(assembly_dic)} assemblies')
    if not assembly_dic:
        logging.error('No assemblies found ! Change search term!')
    return assembly_dic


def _chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]


def get_lineage(taxid, target_ranks):
    tax = {}
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    ranks = ncbi.get_rank(lineage)
    rank_list = list(ranks.values())
    name_list = list(names.values())
    rank2names = dict(zip(rank_list, name_list))
    for sub_taxid in lineage:
        rank = ranks[sub_taxid]
        if rank in target_ranks:
            tax[f"{rank}_taxid"] = str(int(sub_taxid))
            for n, rank in enumerate(target_ranks):
                if rank in rank2names.keys():
                    tax[rank] = rank2names[rank]
                else:
                    previous_name = 'root'
                    for ranks_to_skip in range(1, len(target_ranks) + 1):
                        previous_rank = target_ranks[n - ranks_to_skip]
                        if previous_rank in rank2names.keys():
                            previous_name = rank2names[previous_rank]
                            break
                        if previous_rank not in rank2names.keys():
                            continue
                    tax[rank] = f'{previous_name}_{rank[0:1]}'
    return tax


def generate_table(assembly_list):
    assembly_list = ','.join(assembly_list)
    assembly_summary = Entrez.read(Entrez.esummary(db='assembly', id=assembly_list), validate=False)
    target_ranks = ['species', 'genus', 'family', 'order', 'phylum', 'superkingdom']
    dico = {}
    assembly_ID = [assembly_summary['DocumentSummarySet']['DocumentSummary'][i].attributes['uid'] for i in
                   range(len(assembly_summary['DocumentSummarySet']['DocumentSummary']))]

    for i in range(len(assembly_ID)):
        dico[assembly_ID[i]] = {}
        taxonomy_ID = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['Taxid']
        Refseq_cat = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['RefSeq_category']
        Genbank_path = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['FtpPath_GenBank']
        contig_count = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['Meta'].split(
            ' <Stat category="contig_count" sequence_tag="all">')[1].split('</Stat>')[0]
        assembly_length = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['Meta'].split(
            ' <Stat category="total_length" sequence_tag="all">')[1].split('</Stat>')[0]
        Genbank_release_date = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['AsmReleaseDate_GenBank']
        assembly_stat = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['AssemblyStatus']
        Refseq_path = assembly_summary['DocumentSummarySet']['DocumentSummary'][i]['FtpPath_RefSeq']
        dico[assembly_ID[i]]['AssemblyStatus'] = assembly_stat
        dico[assembly_ID[i]]['Refseq_category'] = Refseq_cat
        dico[assembly_ID[i]]['Contig_count'] = contig_count
        dico[assembly_ID[i]]['Assembly_length'] = assembly_length
        dico[assembly_ID[i]]['Release_date_Genbank'] = Genbank_release_date
        dico[assembly_ID[i]]['FtpPath_Refseq'] = Refseq_path
        dico[assembly_ID[i]]['FtpPath_Genbank'] = Genbank_path
        dico[assembly_ID[i]]['Taxid'] = taxonomy_ID
        taxdic = get_lineage(taxonomy_ID, target_ranks)
        for rank in taxdic:
            dico[assembly_ID[i]][rank] = taxdic[rank]

    assembly_table = pd.DataFrame.from_dict(dico, orient='index')
    assemtb = assembly_table[['AssemblyStatus', 'Refseq_category', 'Contig_count', 'Assembly_length',
                              'Release_date_Genbank','FtpPath_Refseq', 'FtpPath_Genbank', 'superkingdom',
                              'phylum', 'order', 'family', 'genus', 'species', 'Taxid']]
    return assemtb


'''
Main
'''
Entrez.email = snakemake.params['NCBI_email']
Entrez.api_key = snakemake.params['NCBI_key']
comp = snakemake.params['comp']
ref = snakemake.params['ref']
rep = snakemake.params['rep']
met = snakemake.params['met']
gb = snakemake.params['gb']
rs = snakemake.params['rs']
taxid = taxid_find(snakemake.wildcards.entry)
assemblies_found=search_assemblies(taxid, Genbank=gb, Refseq=rs,complete=comp,
                                   reference=ref, representative=rep, exclude_metagenomes=met)
if len(assemblies_found) > 300:
    assemblies_chunks = _chunks(assemblies_found, 300)  # Divide assembly lists by chunks of 300 if more than 300 found
    logging.info('Parsing assemblies by chucks of 300')
    table_chunks = []
    for n, chunks in enumerate(assemblies_chunks):
        logging.info(f'chunk n°{n}')
        tb = generate_table(chunks)
        table_chunks.append(tb)
        # time.sleep(1)
    all_tb = pd.concat(table_chunks, sort=False)
else:
    all_tb = generate_table(assemblies_found)

all_tb.index.name = 'AssemblyID'
all_tb.to_csv(snakemake.output[0], sep='\t')


