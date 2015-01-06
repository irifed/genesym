from pprint import pprint
import time
import re

import rpy2.robjects as robjects
import pandas.rpy.common as com
import pandas

import logging

logger = logging.getLogger('genesym')

data_dir = '/Users/irina/Projects/insilico medicine/datasets/'

def get_gpl_r(gpl_id):
    tic = time.time()
    robjects.r('''
        suppressMessages(library(GEOquery))
        destdir <- "{}"
        gpl <- getGEO("{}", destdir = destdir)
        gpldf <- Table(gpl)
    '''.format(data_dir, gpl_id))
    toc = time.time()
    logging.debug('getGEO took {} seconds'.format(toc - tic))

    tic = time.time()
    gpldf = com.load_data('gpldf')
    toc = time.time()
    logging.debug('load_data took {} seconds'.format(toc - tic))

    return gpldf

def get_gpl(gpl_id):
    tic = time.time()

    # find start and stop of platform table in SOFT file
    table_begin_ln = 0
    table_end_ln = 0
    soft_fn = data_dir + gpl_id + '.soft'
    with open(soft_fn) as soft_file:
        for num, line in enumerate(soft_file, 1):
            if '!platform_table_begin' in line:
                table_begin_ln = num
            elif '!platform_table_end' in line:
                table_end_ln = num

    # load platform table as dataframe
    skiprows = list(range(0, table_begin_ln)) + [table_end_ln-1]
    gpldf = pandas.read_csv(soft_fn, skiprows=skiprows, delimiter='\t')

    toc = time.time()
    logging.debug('loading GPL data table took {} seconds'.format(toc - tic))

    return gpldf

def get_hgnc_symbol(row):
    # get gene symbol value
    # if it does not exist, get entrez id and lookup in biomart
    # if it does not exist, get unigene id and lookup in biomart
    # ... genbank accession
    # ...


    gene_symbol = get_gene_symbol(row)
    if gene_symbol is not None:
        hgnc_symbol = hgnc_lookup(gene_symbol)
        return hgnc_symbol

    else:
        # lookup other ids
        return 'other ids'


def get_gene_symbol(row):
    # TODO check: 'symbol' or 'gene symbol' or other symbol name must be present or not?
    gene_symbol = None

    if 'Gene Symbol' in row:
        gene_symbol = row['Gene Symbol']
    elif 'Symbol' in row:
        gene_symbol = row['Symbol']

    return gene_symbol


def hgnc_lookup(gene_symbol):
    # TODO implement real HGNC lookup
    return gene_symbol + '_HGNC'


def get_symbol_from_entrez_id(entrez_id):
    return 'foo_entrez_hgnc'