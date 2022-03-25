from Bio import Entrez
import pandas as pd
import warnings
import numpy as np
import logging
from ete3 import NCBITaxa
ncbi = NCBITaxa()


class AssemblyFinder:
    def __init__(self, name, isassembly=False, db='', source='latest[filter]', category=[], assembly_level = [], exclude=[], annotation=True,
                    nb=1, rank_to_select=False, outf='f.tsv', outnf='nf.tsv', n_by_rank=1):
        self.name = name
        self.db = db
        self.rcat = category
        self.alvl = assembly_level
        self.excl = exclude
        self.annot = annotation
        self.target_ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
        self.nchunks = 10000
        self.rank_to_select = rank_to_select
        self.nb = nb
        self.outf = outf
        self.outnf = outnf
        self.n_by_rank = n_by_rank
        logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', datefmt='%d %b %Y %H:%M:%S',
                            filename=snakemake.log[0], level=logging.DEBUG)

    # Static methods to apply functions on assembly summary table
    def get_stat(self, meta, stat):
        """
        function to extract assembly Meta stats (contig count, assembly length)
        """
        return meta.split(f' <Stat category="{stat}" sequence_tag="all">')[1].split('</Stat>')[0]

    def get_names(self, gbftp):
        """
        function to extract assembly file names
        """
        return gbftp.split('/')[-1]

    def get_lin_tax(self, lineages):
        """
        function to get lineages from a list of taxids
        """
        ranks = ncbi.get_rank(lineages).values()
        ranknames = ncbi.get_taxid_translator(lineages).values()
        return dict(zip(ranks, ranknames))

    def replace_nans(self, tb):
        """
        function to replace unknown taxonomic rank with placeholder names
        """
        tb = tb.replace(np.nan, 'unknown')
        for i in range(len(tb)):
            for n, col in enumerate(tb.columns):
                if tb.iloc[i, n] == 'unknown' and col != 'superkingdom':
                    tmpname = tb.iloc[i, n - 1] + '_' + col[0]
                    if col == 'species':
                        tmpname = tb.iloc[i, n - 1] + '_' + col[0:1]
                    tb.iloc[i, n] = tmpname
        return tb

    def chunks(self, ls, n):
        """
        function to split assembly list into chunks
        """
        return [ls[i:i + n] for i in range(0, len(ls), n)]

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
        if self.assembly:  # If the input is an assembly name or Gbuid use it as a search term
            search_term = f'{self.name}'
        else:  # If not, search ncbi taxonomy for the taxid
            taxid = self.taxid_find()
            if len(self.db)>0:
                self.source = f'"latest {db}"[filter]'
            if len(self.rcat)>1:
                refseq_category = f' AND ('
                for n, r in enumerate(self.rcat):
                    if n+1 == len(self.rcat):    
                        refseq_category += f'"{r} genome"[filter]'
                    else:
                        refseq_category += f'"{r} genome"[filter] OR '
                refseq_category += ')'
            elif len(self.rcat) == 1:
                refseq_category = f' AND "{rcat[0]} genome"[filter]'
            if len(self.alvl) > 1:
                assembly_level = ' AND ('
                for n, a in enumerate(self.alvl):
                    if n+1 == len(self.alvl):
                        assembly_level += f'"{a}"[filter]'
                    else:
                        assembly_level += f'"{a}"[filter] OR '
                assembly_level += ')'
            elif len(self.alvl) == 1:
                assembly_level = f' AND "{alvl[0]}"[filter]'
            if len(self.excl)>0:
                exclude = ' AND all[filter]'
                for e in self.excl:
                    exclude += f' NOT {e}[filter]'
            if self.annot and (len(db)>0):
                annotation = f' AND "{db} has annotation"[Properties]'
            elif self.annot and (len(db)<=0):
                annotation = f' AND "has annotation"[Properties]'
            search_term = f'txid{taxid}[Organism:exp] AND ({source}{refseq_category}{assembly_level}{exclude}{annotation})'
        assembly_ids = Entrez.read(Entrez.esearch(db='assembly', term=search_term, retmax=2000000))['IdList']
        logging.info(f'> Search term: {search_term}')
        logging.info(f'found {len(assembly_ids)} assemblies')
        if not assembly_ids:
            raise Warning('No assemblies found ! Change search term!')
        return assembly_ids

    def generate_assembly_table(self, assemblies):
        assembly_list = ','.join(assemblies)
        assembly_summary = Entrez.read(Entrez.esummary(db='assembly', id=assembly_list), validate=False)
        summaries = assembly_summary['DocumentSummarySet']['DocumentSummary']
        tb = pd.DataFrame.from_records(summaries)
        columns = ['GbUid', 'RefSeq_category', 'AssemblyStatus', 'FtpPath_GenBank', 'FtpPath_RefSeq', 'Meta',
                   'AsmReleaseDate_GenBank', 'ContigN50', 'ScaffoldN50', 'Coverage', 'Taxid']
        subset = tb[columns]
        lens = subset.apply(lambda x: self.get_stat(x['Meta'], stat='total_length'), axis=1)
        contigs = subset.apply(lambda x: self.get_stat(x['Meta'], stat='contig_count'), axis=1)
        subset.insert(loc=subset.shape[1] - 1, value=lens, column='Assembly_length')
        subset.insert(loc=subset.shape[1] - 1, value=contigs, column='Contig_count')
        subset.insert(loc=1, value=subset['FtpPath_GenBank'].apply(self.get_names), column='AssemblyNames')
        subset = subset.rename(columns={'Coverage': 'Assembly_coverage'})
        subset = subset.drop('Meta', axis=1)
        return subset

    def add_lineage(self, assembly_tb):
        unique_taxids = list(set(assembly_tb['Taxid']))
        taxid2lineage = ncbi.get_lineage_translator(unique_taxids)
        tax = {taxid: self.get_lin_tax(lineage) for taxid, lineage in taxid2lineage.items()}
        lineage_tb = pd.DataFrame.from_dict(tax, orient='index')
        lineage_tb.index.set_names('Taxid', inplace=True)
        lineage_tb.reset_index(inplace=True)
        ordered_ranks = self.target_ranks[::-1]
        ordered_ranks.append('Taxid')
        lin_cols = list(lineage_tb.columns)
        all_cols = list(set().union(lin_cols, ordered_ranks))
        lineage_tb = lineage_tb.reindex(columns=all_cols, fill_value=np.nan)
        lineage_tb = lineage_tb[ordered_ranks]
        lineage_tb = self.replace_nans(lineage_tb)
        lineage_tb = lineage_tb.astype({'Taxid': 'string'})
        merged_table = assembly_tb.merge(lineage_tb, on='Taxid')
        return merged_table

    def select_assemblies(self, table):
        fact_table = table.replace({'RefSeq_category': {'reference genome': 0, 'representative genome': 1, 'na': 6},
                                    'AssemblyStatus': {'Complete Genome': 2, 'Chromosome': 3, 'Scaffold': 4,
                                                       'Contig': 5, 'na': 6}})
        # make sure the path to an ftp is available
        if self.db == 'refseq':
            fact_table = fact_table[fact_table.FtpPath_Refseq != '']
        else:
            fact_table == fact_table[fact_table.FtpPath_GenBank != '']
        sorted_table = fact_table.sort_values(['RefSeq_category', 'AssemblyStatus', 'Contig_count',
                                               'ScaffoldN50', 'ContigN50', 'AsmReleaseDate_GenBank'],
                                              ascending=[True, True, True, False, False, False]).replace(
                                                  {'RefSeq_category': {0: 'reference genome', 1: 'representative genome',6: 'na'},
                                                  'AssemblyStatus': {2: 'Complete Genome', 3: 'Chromosome', 4: 'Scaffold',5: 'Contig', 6: 'na'}})
        
        if self.rank_to_select:
            logging.info(f'Filtering according to {self.rank_to_select}, Refseq categories, assembly status, '
                         f'contig count and release date')
            uniq_rank = set(sorted_table[f'{rank}'])
            sorted_table = pd.concat([sorted_table[f'{rank}' == ranks].iloc[:self.n_by_rank] for ranks in uniq_rank])
        if isinstance(self.nb, int):
            if len(sorted_table) >= self.nb:
                logging.info(f'Selecting {self.nb} from every {self.rank_to_select} out of {len(sorted_table)} assemblies')

            if len(sorted_table) < self.nb:
                logging.warning(f'Found less than {self.nb} assemblies in total, returning {len(sorted_table)} instead')
        else:
            logging.warning(f'No assembly number specified returning all assemblies found')
        return sorted_table

    def run(self):
        assemblies_found = self.search_assemblies()
        if len(assemblies_found) > self.nchunks:
            raise Warning(f'{len(assemblies_found)} assemblies found, restrict search term to find less assemblies')
            assemblies_chunks = self.chunks(assemblies_found,
                                            self.nchunks)  # Divide assembly lists by chunks of 10000
            logging.info(f'Parsing assemblies by chucks of {self.nchunks}')
            table_chunks = []
            for n, chunk in enumerate(assemblies_chunks):
                logging.info(f'chunk n°{n}')
                assembly_tb = self.generate_assembly_table(chunk)
                tb = self.add_lineage(assembly_tb)
                table_chunks.append(tb)
            non_filtered_tb = pd.concat(table_chunks, sort=False)
        else:
            assembly_tb = self.generate_assembly_table(assemblies_found)
            non_filtered_tb = self.add_lineage(assembly_tb)

        non_filtered_tb.to_csv(self.outnf, sep='\t', index=None)
        filtered_tb = self.select_assemblies(non_filtered_tb)
        filtered_tb.to_csv(self.outf, sep='\t', index=None)
        return filtered_tb


'''
Main
'''
Entrez.email = snakemake.params['ncbi_email']
Entrez.api_key = snakemake.params['ncbi_key']
entry = snakemake.wildcards.entry
assembly = snakemake.params['assembly']
db = snakemake.params['db']
cat = snakemake.params['rcat']
alvl = snakemake.params['alvl']
excl = snakemake.params['excl']
annot = snakemake.params['annot']
column = snakemake.params['column']
rank = snakemake.params['rank']
intb = pd.read_csv(snakemake.input[0], sep='\t', dtype={f'{column}': 'str'})
intb.set_index(f'{column}', inplace=True)
if assembly:
    nb = 1
else:
    try:
        nb = int(intb.loc[entry]['NbGenomes'])
    except KeyError:
        nb = None
n_by_rank = snakemake.params['n_by_rank']
find_assemblies = AssemblyFinder(name=entry, isassembly=assembly, db=db, category=cat, assembly_level=alvl, exclude=excl,
                                 nb=nb, rank_to_select=rank, annotation = annot, n_by_rank=n_by_rank,
                                 outnf=snakemake.output.all, outf=snakemake.output.filtered, )
find_assemblies.run()
