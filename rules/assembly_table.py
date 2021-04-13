from Bio import Entrez
import pandas as pd
import logging
from ete3 import NCBITaxa
ncbi = NCBITaxa()


class AssemblyFinder:
    def __init__(self, name, genbank=False, refseq=True, representative=True, reference=True, complete=True,
                 exclude_metagenomes=True, nb=1, rank_to_select='None', outf='f.tsv', outnf='nf.tsv'):
        self.name = name
        self.genbank = genbank
        self.refseq = refseq
        self.representative = representative
        self.reference = reference
        self.complete = complete
        self.exclude_metagenomes = exclude_metagenomes
        self.target_ranks = ['species', 'genus', 'family', 'order', 'phylum', 'superkingdom']
        self.nchunks = 10000
        self.rank_to_select = rank_to_select
        self.nb = nb
        self.outf = outf
        self.outnf = outnf
        logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', datefmt='%d %b %Y %H:%M:%S',
                            filename=snakemake.log[0], level=logging.DEBUG)
        logging.captureWarnings(True)

    def taxid_find(self):
        """
        Function to parse input_list to convert scientific names to taxid
        returns dictionary with taxid found
        """
        logging.info(f'> Searching for taxIDs {self.name} ...')
        try:
            int(self.name)
            logging.info('Query is a taxID')
            taxid = self.name

        except ValueError:
            logging.warning('Query is not a taxID, enter taxID to be more precise')
            logging.info(f'Search term: {self.name}[all Names]')
            taxid_list = Entrez.read(Entrez.esearch(db='taxonomy', term=f'{self.name}[all Names]', retmax=100))[
                'IdList']
            if len(taxid_list) == 1:
                taxid = taxid_list[0]
                logging.info(f'TaxID:{taxid} found')
            if len(taxid_list) > 1:
                taxid = taxid_list[0]
                logging.warning(f'{len(taxid_list)} TaxIDs found, change query (taking first one : {taxid})')
            if len(taxid_list) == 0:
                raise Exception('TaxID not found! Change search term!')
        return taxid

    def search_assemblies(self):
        taxid = self.taxid_find()
        search_term = f'txid{taxid}[Organism:exp] '
        if self.refseq and not self.genbank:
            search_term += 'AND("latest refseq"[filter]) '
        if self.genbank and not self.refseq:
            search_term += 'AND("latest genbank"[filter]) '
        if self.genbank and self.refseq:
            search_term += 'AND ("latest genbank"[filter] OR "latest refseq"[filter]) '
        if self.complete and not self.representative and not self.reference:
            search_term += 'AND ("complete genome"[filter]) '
        if self.complete and self.representative and not self.reference:
            search_term += 'AND ("complete genome"[filter]) OR ("representative genome"[filter]) '
        if self.complete and self.representative and self.reference:
            search_term += 'AND ("complete genome"[filter]) OR ("representative genome"[filter]) OR ' \
                           '("reference genome"[filter]) '
        if self.representative and not self.reference:
            search_term += 'AND ("representative genome"[filter]) '
        if self.reference and not self.representative:
            search_term += 'AND ("reference genome"[filter]) '
        if self.representative and self.reference:
            search_term += 'AND ("representative genome"[filter] OR "reference genome"[filter]) '
        if self.exclude_metagenomes:
            search_term += 'AND (all[filter] NOT "derived from metagenome"[filter])'
        assembly_ids = Entrez.read(Entrez.esearch(db='assembly', term=search_term, retmax=500000))['IdList']
        logging.info(f'> Search term: {search_term}')
        logging.info(f'found {len(assembly_ids)} assemblies')
        if not assembly_ids:
            raise Exception('No assemblies found ! Change search term!')
        return assembly_ids

    def get_lineage(self, txid, target_ranks):
        tax = {txid: {'Taxid': txid}}
        lineage = ncbi.get_lineage(txid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        rank_list = list(ranks.values())
        name_list = list(names.values())
        rank2names = dict(zip(rank_list, name_list))
        for lineage_taxid in lineage:
            rank = ranks[lineage_taxid]
            if rank in target_ranks:
                tax[txid][f"{rank}_taxid"] = str(int(lineage_taxid))
                for n, rank in enumerate(target_ranks):
                    if rank in rank2names.keys():
                        tax[txid][rank] = rank2names[rank]
                    else:
                        previous_name = 'root'
                        for ranks_to_skip in range(1, len(target_ranks) + 1):
                            previous_rank = target_ranks[n - ranks_to_skip]
                            if previous_rank in rank2names.keys():
                                previous_name = rank2names[previous_rank]
                                break
                            if previous_rank not in rank2names.keys():
                                continue
                        tax[txid][rank] = f'{previous_name}_{rank[0:1]}'
        return tax

    def generate_table(self, assemblies):
        assembly_list = ','.join(assemblies)
        assembly_summary = Entrez.read(Entrez.esummary(db='assembly', id=assembly_list), validate=False)
        summaries = assembly_summary['DocumentSummarySet']['DocumentSummary']
        uids = [summ.attributes['uid'] for summ in summaries]
        taxids = [summ['Taxid'] for summ in summaries]
        names = [summ['AssemblyName'] for summ in summaries]
        cat = [summ['RefSeq_category'] for summ in summaries]
        status = [summ['AssemblyStatus'] for summ in summaries]
        gbftp = [summ['FtpPath_GenBank'] for summ in summaries]
        rsftp = [summ['FtpPath_RefSeq'] for summ in summaries]
        unique_taxids = list(set(taxids))
        linear_taxonomy = [self.get_lineage(atxid, self.target_ranks) for atxid in unique_taxids]
        contigs = [int(summ['Meta'].split(
            ' <Stat category="contig_count" sequence_tag="all">')[1].split(
            '</Stat>')[0]) for summ in summaries]
        lens = [int(summ['Meta'].split(
            ' <Stat category="total_length" sequence_tag="all">')[1].split(
            '</Stat>')[0]) for summ in summaries]
        dates = [summ['AsmReleaseDate_GenBank'] for summ in summaries]
        assemblies_tb = pd.DataFrame(list(zip(uids, names, status, cat, contigs, lens, dates, rsftp, gbftp, taxids)),
                                     columns=['AssemblyID', 'AssemblyNames', 'AssemblyStatus', 'Refseq_category',
                                              'Contig_count', 'Assembly_length', 'Release_date_Genbank',
                                              'FtpPath_Refseq', 'FtpPath_Genbank', 'Taxid'])
        assemblies_tax = [pd.DataFrame.from_dict(lint, orient='index') for lint in linear_taxonomy]
        assemblies_tax = pd.concat(assemblies_tax)
        merged_table = pd.merge(assemblies_tb, assemblies_tax, on='Taxid')
        return merged_table

    def chunks(self, lst, n):
        return [lst[i:i + n] for i in range(0, len(lst), n)]

    def select_assemblies(self, table):

        fact_table = table.replace({'Refseq_category': {'reference genome': 0, 'representative genome': 1, 'na': 6},
                                    'AssemblyStatus': {'Complete Genome': 2, 'Chromosome': 3, 'Scaffold': 4,
                                                       'Contig': 5, 'na': 6}})
        sorted_table = fact_table.sort_values(
            ['Refseq_category', 'AssemblyStatus', 'Contig_count', 'Release_date_Genbank'],
            ascending=[True, True, True, False])

        if self.rank_to_select != 'None':
            logging.info(f'Filtering according to {self.rank_to_select}, Refseq categories, assembly status, '
                         f'contig count and release date')
            select_index = []
            unique_list = list(set(sorted_table[self.rank_to_select]))
            if len(unique_list) > 1:
                for i in unique_list:
                    select_index.append(sorted_table[sorted_table[self.rank_to_select] == i].sample(1).index[0])
                    # randomly select one assembly ID for each unique selected rank (species for example)
                sorted_table = sorted_table.loc[select_index, :]
            if len(unique_list) == 1:
                logging.info(f'Same {self.rank_to_select} for all assemblies, Filtering according to Refseq '
                             f'categories, assembly status,contig count and release date')
            if len(unique_list) == 0:
                logging.error(f'{self.rank_to_select} is not a target rank')
        else:
            logging.info('No taxonomic rank specified, sorting according to Refseq category, '
                         'assembly status, contig count and release date')
        if len(sorted_table) >= self.nb:
            logging.info(f'Selecting {self.nb} sorted assemblies out of {len(sorted_table)}')
            sorted_table = sorted_table[0:self.nb]
        if len(sorted_table) < self.nb:
            logging.warning(f'Found less than {self.nb} assemblies in total, returning {len(sorted_table)} instead')
        sorted_table = sorted_table.replace({'Refseq_category': {0: 'reference genome', 1: 'representative genome',
                                                                 6: 'na'},
                                             'AssemblyStatus': {2: 'Complete Genome', 3: 'Chromosome', 4: 'Scaffold',
                                                                5: 'Contig', 6: 'na'}})
        return sorted_table

    def run(self):
        assemblies_found = self.search_assemblies()
        if len(assemblies_found) > self.nchunks:
            assemblies_chunks = self.chunks(assemblies_found,
                                            self.nchunks)  # Divide assembly lists by chunks of 10000
            logging.info(f'Parsing assemblies by chucks of {self.nchunks}')
            table_chunks = []
            for n, chunk in enumerate(assemblies_chunks):
                logging.info(f'chunk n°{n}')
                tb = self.generate_table(chunk)
                table_chunks.append(tb)
            non_filtered_tb = pd.concat(table_chunks, sort=False)
        else:
            non_filtered_tb = self.generate_table(assemblies_found)

        non_filtered_tb.to_csv(self.outnf, sep='\t', index=None)
        filtered_tb = self.select_assemblies(non_filtered_tb)
        filtered_tb.to_csv(self.outf, sep='\t', index=None)
        return filtered_tb


'''
Main
'''
Entrez.email = snakemake.params['ncbi_email']
Entrez.api_key = snakemake.params['ncbi_key']
comp = snakemake.params['comp']
ref = snakemake.params['ref']
rep = snakemake.params['rep']
met = snakemake.params['met']
gb = snakemake.params['gb']
rs = snakemake.params['rs']
entry = snakemake.wildcards.entry
intb = pd.read_csv(snakemake.input[0], sep='\t', dtype={'UserInputNames': 'str'})
intb.set_index('UserInputNames', inplace=True)
nb = intb.loc[entry]['nb_genomes']
rank = snakemake.params['rank_filter']
find_assemblies = AssemblyFinder(name=entry, genbank=gb, refseq=rs, representative=rep, reference=ref,
                                 complete=comp, exclude_metagenomes=met, nb=nb, rank_to_select=rank,
                                 outnf=snakemake.output.all, outf=snakemake.output.filtered)
find_assemblies.run()
