import pandas as pd
import json
import os
import glob
import numpy as np


# Functions
def read_json(file):
    try:
        return pd.json_normalize(json.load(open(file)), record_path=["reports"])
    except KeyError:
        return pd.read_json(file)


def convert_query(wildcards):
    try:
        return int(wildcards.query)
    except ValueError:
        if ("_" in wildcards.query) and (
            ("GCF" not in wildcards.query) and ("GCA" not in wildcards.query)
        ):
            return wildcards.query.replace("_", " ")
        else:
            return wildcards.query


def get_limit(wildcards, nbs, dic):
    if nbs:
        return dic[wildcards.query]
    else:
        return ""


def get_abs_path(indir, accessions):
    return np.array(
        [
            os.path.abspath(glob.glob(os.path.join(indir, acc, f"{acc}*.fna*"))[0])
            for acc in accessions
        ]
    )


# ARGS
KEY = ""
if API_KEY:
    KEY = f"--api-key {API_KEY}"
ARGS = ""
if not TAXON:
    ARGS = ""
else:
    if ANNOTATED:
        ARGS += "--annotated "
    if ASM_LVL:
        ARGS += f"--assembly-level {ASM_LVL} "
    if SOURCE:
        ARGS += f"--assembly-source {SOURCE} "
    if ATYPICAL:
        ARGS += "--exclude-atypical "
    if MAG:
        ARGS += f"--mag {MAG} "
    if REFERENCE:
        ARGS += "--reference "

ARGS = ARGS.strip()

GZIP = ""
if COMPRESSED:
    GZIP = "--gzip"

QUERY2NB = {}

if TAXON:
    try:
        df = pd.read_csv(INPUT, sep="\t")
        QUERIES = list(df["taxon"])
        if "nb" in df.columns:
            LIMIT = list(df["nb"])
    except (FileNotFoundError, IsADirectoryError):
        QUERIES = INPUT.split(",")

    QUERIES = [str(query) for query in QUERIES]
    if LIMIT:
        if type(LIMIT) is not list:
            LIMIT = LIMIT.split(",")
        LIMIT = [f"--limit {nb}" for nb in LIMIT]
        if len(LIMIT) == 1:
            LIMIT = LIMIT * len(QUERIES)
        QUERY2NB = dict(zip(QUERIES, LIMIT))
