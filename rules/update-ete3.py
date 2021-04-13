#!/usr/bin/env python
import sys
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

from argparse import ArgumentParser
import os
from datetime import datetime
from ete3 import NCBITaxa
ncbi = NCBITaxa()

def update_ete3(env, out):
    sqldb = os.path.join(env, '.etetoolkit', 'taxa.sqlite')  # path to ete sql db
    db_modification_time = datetime.fromtimestamp(os.path.getctime(sqldb))
    database_age_days = abs((db_modification_time-datetime.now()).days)
    if database_age_days >= 5:
        ncbi.update_taxonomy_database()
        comment = f'taxa.sqlite is more than {database_age_days} days old, updating database'
    else:
        comment = 'taxa.sqlite is up to date'
    file = open(out, 'w')
    file.write(comment)
    file.close()


def update_ete3_args(parser):
    parser.add_argument('-i', '--home', required=True, help='path to home')
    parser.add_argument('-o', '--out', required=True, help='path to ouput file')
    return parser


def update_ete3_run(opts):
    update_ete3(opts.home, opts.out)


if __name__ == '__main__':
    parser = ArgumentParser(description='Initialise ete3 taxa.sqlite and check for updates')
    parser = update_ete3_args(parser)
    opts, unknown_args = parser.parse_known_args()
    update_ete3_run(opts)
