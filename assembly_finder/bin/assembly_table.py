from Bio import Entrez
import pandas as pd
import warnings
import numpy as np
import logging
from ete3 import NCBITaxa


class AssemblyFinder:
    def __init__(
        self,
        name,
        uid="False",
        db="refseq",
        source="latest[filter]",
        category="",
        assembly_level=[],
        exclude="",
        annotation="False",
        nb="all",
        rank_to_select=None,
        outf="f.tsv",
        outnf="nf.tsv",
        n_by_rank="none",
        release="AsmReleaseDate_RefSeq",
    ):
        self.n_by_rank = n_by_rank
        self.name = name.replace("-", " ")
        self.db = db
        self.nb = nb
        self.uid = eval(uid)
        self.release = release

        if self.db == "refseq":
            self.dbuid = "RsUid"
            self.ftp = "FtpPath_RefSeq"
            self.release = "AsmReleaseDate_RefSeq"

        else:
            self.db = "genbank"
            self.dbuid = "GbUid"
            self.ftp = "FtpPath_GenBank"
            self.release = "AsmReleaseDate_GenBank"

        self.columns = [
            "entry",
            self.dbuid,
            "AssemblyAccession",
            "AssemblyName",
            self.release,
            self.ftp,
            "AssemblyStatus",
            "RefSeq_category",
            "ContigCount",
            "ContigL50",
            "ContigN50",
            "Coverage",
            "TotalLength",
            "Taxid",
            "Organism",
            "Sub_type",
            "Sub_value",
        ]

        self.source = source
        self.rcat = category.split("_")
        if assembly_level == "all":
            self.alvl = "all"
        else:
            self.alvl = []
            for level in assembly_level.split("_"):
                if level == "complete":
                    self.alvl.append(f"{level} genome")
                else:
                    self.alvl.append(f"{level} level")

        self.excl = [e.replace("-", " ") for e in exclude.split("_")]
        self.annot = eval(annotation)
        self.target_ranks = [
            "species",
            "genus",
            "family",
            "order",
            "class",
            "phylum",
            "superkingdom",
        ]
        self.nchunks = 10000
        self.rank_to_select = rank_to_select

        try:
            self.nb = int(self.nb)
            self.n_by_rank = int(self.n_by_rank)
        except ValueError:
            self.nb = self.nb
            self.n_by_rank = self.n_by_rank
        self.outf = outf
        self.outnf = outnf

        logging.basicConfig(
            format="%(asctime)s %(levelname)s %(message)s",
            datefmt="%d %b %Y %H:%M:%S",
            filename=snakemake.log[0],
            level=logging.DEBUG,
        )

    # Static methods to apply functions on assembly summary table
    def get_stat(self, meta, stat):
        """
        function to extract assembly Meta stats (contig count, assembly length)
        """
        return meta.split(f' <Stat category="{stat}" sequence_tag="all">')[1].split(
            "</Stat>"
        )[0]

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
        tb = tb.replace(np.nan, "unknown")
        for i in range(len(tb)):
            for n, col in enumerate(tb.columns):
                if tb.iloc[i, n] == "unknown" and col != "superkingdom":
                    tmpname = tb.iloc[i, n - 1] + "_" + col[0]
                    if col == "species":
                        tmpname = tb.iloc[i, n - 1] + "_" + col[0:2]
                    tb.iloc[i, n] = tmpname
        return tb

    def chunks(self, ls, n):
        """
        function to split assembly list into chunks
        """
        return [ls[i : i + n] for i in range(0, len(ls), n)]

    def search_assemblies(self):
        taxid = None
        assembly_ids = []

        if self.name.isnumeric():
            logging.info(
                f"{self.name} is numeric, checking if it is a taxid or assembly uid"
            )
            if self.uid:
                logging.info(f"Checking if {self.name} is a uid")
                assembly_ids = Entrez.read(
                    Entrez.esearch(
                        db="assembly", term=f"{self.name}[uid]", retmax=2000000
                    )
                )["IdList"]
                if assembly_ids:
                    if len(assembly_ids) == 1:
                        logging.info(f"{self.name} is an assembly uid ")
                    elif len(assembly_ids) > 1:
                        raise Warning(
                            f"Found more than 1 uid with {self.name}, use a correct uid as input"
                        )
                    else:
                        raise Warning(f"{self.name} is not a uid")

            else:
                try:
                    taxid = list(ncbi.get_taxid_translator([self.name]).keys())[0]
                    logging.info(f"{self.name} is a taxid")
                except IndexError:
                    logging.error(f"{self.name} not a taxid, give a correct taxid !")

        elif self.name.isalpha():
            logging.info(f"{self.name} is alpha, trying to find a taxid")
            try:
                taxid = list(ncbi.get_name_translator([self.name]).values())[0][0]
                logging.info(f"found taxid {taxid} for {self.name}")
            except IndexError:
                logging.error(f"{self.name} not a correct taxonomic name !")

        elif (not self.name.isalpha()) and (
            self.name.isalnum() or "GCF" in self.name or "ASM" in self.name
        ):
            logging.info(
                f"trying to match {self.name} it to assembly name or accession"
            )
            assembly_ids = Entrez.read(
                Entrez.esearch(db="assembly", term=self.name, retmax=2000000)
            )["IdList"]
            if assembly_ids:
                if len(assembly_ids) == 1:
                    logging.info(f"{self.name} is an assembly name, found its uid")
                elif len(assembly_ids) > 1:
                    assembly_ids = []
                    raise Warning(
                        f"Found more than 1 uid with {self.name} in NCBI assembly, use a more specific name"
                    )
                else:
                    raise Warning(
                        f"{self.name} is alphanumeric but did not correspond to an assembly name, change it !"
                    )
        else:
            logging.info(
                f"{self.name} neither taxid, taxonomic name, uid or accession. Trying to find it in NCBI taxonomy"
            )
            taxid = Entrez.read(
                Entrez.esearch(db="taxonomy", term=self.name, retmax=2000000)
            )["IdList"]
            if len(taxid) == 1:
                taxid = taxid[0]
                logging.info(f"Found unique taxid: {taxid} for {self.name}")
            elif len(taxid) > 1:
                raise Warning(
                    f"Found more than 1 taxid for {self.name}, change it to a more precies entry"
                )
            else:
                raise Warning(f"Did not find a taxid for {self.name}")

        if taxid:
            annotation = ""
            refseq_category = ""
            exclude = ""
            assembly_level = ""
            if len(self.db) != "all":
                self.source = f'"latest {db}"[filter]'

            if len(self.rcat) > 1:
                refseq_category = f" AND ("
                for n, r in enumerate(self.rcat):
                    if n + 1 == len(self.rcat):
                        refseq_category += f'"{r} genome"[filter]'
                    else:
                        refseq_category += f'"{r} genome"[filter] OR '
                refseq_category += ")"

            elif len(self.rcat) == 1 and (self.rcat[0] != "all"):
                refseq_category = f' AND "{self.rcat[0]} genome"[filter]'

            if len(self.alvl) > 1 and (self.alvl != "all"):
                assembly_level = " AND ("
                for n, a in enumerate(self.alvl):
                    if n + 1 == len(self.alvl):
                        assembly_level += f'"{a}"[filter]'
                    else:
                        assembly_level += f'"{a}"[filter] OR '
                assembly_level += ")"
            elif len(self.alvl) == 1 and (self.alvl != "all"):
                assembly_level = f' AND "{self.alvl[0]}"[filter]'

            if len(self.excl) > 0:
                exclude = " AND all[filter]"
                for e in self.excl:
                    exclude += f" NOT {e}[filter]"
            if self.annot and (len(db) > 0):
                annotation = f' AND "{db} has annotation"[Properties]'
            elif self.annot and (len(db) <= 0):
                annotation = f' AND "has annotation"[Properties]'

            search_term = f"txid{taxid}[Organism:exp] AND ({self.source}{refseq_category}{assembly_level}{exclude}{annotation})"
            logging.info(f"> Search term: {search_term}")
            assembly_ids = Entrez.read(
                Entrez.esearch(db="assembly", term=search_term, retmax=2000000)
            )["IdList"]
        logging.info(f"found {len(assembly_ids)} assemblies")
        if not assembly_ids:
            raise Warning("No assemblies found ! Change search term!")
        return assembly_ids

    def generate_assembly_table(self, assemblies):
        assembly_list = ",".join(assemblies)
        assembly_summary = Entrez.read(
            Entrez.esummary(db="assembly", id=assembly_list), validate=False
        )
        tb = pd.DataFrame.from_records(
            assembly_summary["DocumentSummarySet"]["DocumentSummary"]
        )

        lens = tb.apply(lambda x: self.get_stat(x["Meta"], stat="total_length"), axis=1)
        contigs = tb.apply(
            lambda x: self.get_stat(x["Meta"], stat="contig_count"), axis=1
        )
        l50 = tb.apply(lambda x: self.get_stat(x["Meta"], stat="contig_l50"), axis=1)
        tb.insert(loc=tb.shape[1] - 1, value=lens, column="TotalLength")
        tb.insert(loc=tb.shape[1] - 1, value=contigs, column="ContigCount")
        tb.insert(loc=tb.shape[1] - 1, value=l50, column="ContigL50")
        sub_types = []
        sub_values = []
        for biosource in tb.Biosource:
            try:
                sub_types.append(biosource["InfraspeciesList"][0]["Sub_type"])
                sub_values.append(biosource["InfraspeciesList"][0]["Sub_value"])
            except IndexError:
                sub_types.append("na")
                sub_values.append("na")
        tb["Sub_type"] = sub_types
        tb["Sub_value"] = sub_values
        return tb

    def add_lineage(self, assembly_tb):
        unique_taxids = list(set(assembly_tb["Taxid"]))
        taxid2lineage = ncbi.get_lineage_translator(unique_taxids)
        tax = {
            taxid: self.get_lin_tax(lineage) for taxid, lineage in taxid2lineage.items()
        }
        lineage_tb = pd.DataFrame.from_dict(tax, orient="index")
        lineage_tb.index.set_names("Taxid", inplace=True)
        lineage_tb.reset_index(inplace=True)
        ordered_ranks = self.target_ranks[::-1]
        ordered_ranks.append("Taxid")
        lin_cols = list(lineage_tb.columns)
        all_cols = list(set().union(lin_cols, ordered_ranks))
        lineage_tb = lineage_tb.reindex(columns=all_cols, fill_value=np.nan)
        lineage_tb = lineage_tb[ordered_ranks]
        lineage_tb = self.replace_nans(lineage_tb)
        lineage_tb = lineage_tb.astype({"Taxid": "string"})
        merged_table = assembly_tb.merge(lineage_tb, on="Taxid")
        return merged_table

    def select_assemblies(self, table):
        table.drop("Meta", axis=1, inplace=True)
        table.insert(loc=0, value=[self.name] * len(table), column="entry")
        table["Coverage"] = pd.to_numeric(
            table["Coverage"], errors="coerce"
        )  # replace any non numeric values with NaN
        table["RefSeq_category"] = pd.Categorical(
            table["RefSeq_category"],
            ["reference genome", "representative genome", "na"],
            ordered=True,
        )
        table["AssemblyStatus"] = pd.Categorical(
            table["AssemblyStatus"],
            ["Complete Genome", "Chromosome", "Scaffold", "Contig", "na"],
            ordered=True,
        )
        sorted_table = table.sort_values(
            [
                "RefSeq_category",
                "AssemblyStatus",
                "ContigN50",
                "ContigL50",
                "Coverage",
                "ContigCount",
                "LastUpdateDate",
            ],
            ascending=[True, True, False, True, False, True, False],
        )
        if (self.rank_to_select != "none") and (self.n_by_rank != "none"):
            uniq_rank = set(sorted_table[f"{self.rank_to_select}"])
            sorted_table = pd.concat(
                [
                    sorted_table[sorted_table[f"{self.rank_to_select}"] == ranks].iloc[
                        : int(self.n_by_rank)
                    ]
                    for ranks in uniq_rank
                ]
            )
            logging.info(
                f"Selecting the top {self.n_by_rank} assemblies per {self.rank_to_select} per entry"
            )
        elif self.nb != "all":
            if len(sorted_table) >= int(self.nb):
                logging.info(
                    f"Selecting {self.nb} out of {len(sorted_table)} assemblies"
                )
                sorted_table = sorted_table.iloc[: self.nb]
            if len(sorted_table) < int(self.nb):
                logging.info(
                    f"Found less than {self.nb} assemblies in total, returning {len(sorted_table)} instead"
                )
        else:
            logging.info(f"Returning all assemblies found")

        sel_cols = self.columns + self.target_ranks[::-1]
        subset = sorted_table[sel_cols]
        renamed_cols = [
            "entry",
            "db_uid",
            "asm_accession",
            "asm_name",
            "asm_release_date",
            "ftp_path",
            "asm_status",
            "refseq_category",
            "contig_count",
            "contig_l50",
            "contig_n50",
            "coverage",
            "genome_size",
            "taxid",
            "organism",
            "sub_type",
            "sub_value",
        ] + self.target_ranks[::-1]
        subset.columns = renamed_cols
        subset.insert(loc=1, value=[self.db] * len(subset), column="database")
        return subset.sort_values(by=["entry", "refseq_category"])

    def run(self):
        assemblies_found = self.search_assemblies()
        if len(assemblies_found) > self.nchunks:
            warnings.warn(
                f"{len(assemblies_found)} assemblies found, restrict search term to find less assemblies"
            )
            assemblies_chunks = self.chunks(
                assemblies_found, self.nchunks
            )  # Divide assembly lists by chunks of 10000
            logging.info(f"Parsing assemblies by chucks of {self.nchunks}")
            table_chunks = []
            for n, chunk in enumerate(assemblies_chunks):
                logging.info(f"chunk n°{n}")
                assembly_tb = self.generate_assembly_table(chunk)
                tb = self.add_lineage(assembly_tb)
                table_chunks.append(tb)
            non_filtered_tb = pd.concat(table_chunks, sort=False)
        else:
            assembly_tb = self.generate_assembly_table(assemblies_found)
            non_filtered_tb = self.add_lineage(assembly_tb)

        non_filtered_tb.to_csv(self.outnf, sep="\t", index=None)
        filtered_tb = self.select_assemblies(non_filtered_tb)
        filtered_tb.to_csv(self.outf, sep="\t", index=None)
        return filtered_tb


"""
Main
"""
ncbi = NCBITaxa(dbfile=snakemake.input[0])
if snakemake.params["ncbi_email"] != "none":
    Entrez.email = snakemake.params["ncbi_email"]

if snakemake.params["ncbi_key"] != "none":
    Entrez.api_key = snakemake.params["ncbi_key"]

entry = snakemake.wildcards.entry
db = snakemake.params["db"]
uid = snakemake.params["uid"]
cat = snakemake.params["rcat"]
alvl = snakemake.params["alvl"]
excl = snakemake.params["excl"]
annot = snakemake.params["annot"]
rank = snakemake.params["rank"]
n_by_rank = snakemake.params["n_by_rank"]
nb = snakemake.params["nb"]

find_assemblies = AssemblyFinder(
    name=entry,
    db=db,
    uid=uid,
    category=cat,
    assembly_level=alvl,
    exclude=excl,
    annotation=annot,
    nb=nb,
    rank_to_select=rank,
    n_by_rank=n_by_rank,
    outnf=snakemake.output.all,
    outf=snakemake.output.filtered,
)
find_assemblies.run()
