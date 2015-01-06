from pprint import pprint
import time
import re

import rpy2.robjects as robjects
import pandas.rpy.common as com
import pandas
import numpy

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

def get_hgnc_id_symbol(row):
    tic = time.time()
    # get gene symbol value
    # if it does not exist, get entrez id and lookup in biomart
    # if it does not exist, get unigene id and lookup in biomart
    # ... genbank accession
    # ...

    gene_symbol = get_gene_symbol(row)
    if gene_symbol is not None:
        hgnc_id, hgnc_symbol = hgnc_lookup(gene_symbol)

    else:
        # lookup other ids
        # entrez
        # unigene
        # genbank id GI
        # genbank accession
        hgnc_id, hgnc_symbol = 'TODO_ID', 'TODO_SYMBOL'

    toc = time.time()

    logging.debug('lookup time = {}'.format(toc - tic))

    return hgnc_id, hgnc_symbol

def get_gene_symbol(row):
    # TODO check: 'symbol' or 'gene symbol' or other symbol name must be present or not?
    gene_symbol = None

    if 'Gene Symbol' in row:
        gene_symbol = row['Gene Symbol']
    elif 'Symbol' in row:
        gene_symbol = row['Symbol']
    else:
        logger.warning('Gene Symbol attribute is not found in row {}'.format(row))

    logger.debug('gene_symbol = {}'.format(gene_symbol))

    if type(gene_symbol) == float and numpy.isnan(gene_symbol):
        gene_symbol = None

    return gene_symbol


def hgnc_lookup(gene_symbol):
    hgnc_symbol = None
    hgnc_id = str(gene_symbol) + '_HGNC'

    hgnc_complete_set_fn = '/Users/irina/Projects/insilico medicine/datasets/hgnc_complete_set.txt'

    tic = time.time()
    hgnc_df = pandas.read_csv(hgnc_complete_set_fn, delimiter='\t')
    toc = time.time()
    logger.debug('HGNC file load took {} sec'.format(toc - tic))

    cols = ['Approved Symbol', 'Previous Symbols', 'Synonyms', 'HGNC ID']
    hgnc_df = hgnc_df[cols]

    found = False

    for col in cols:
        df = hgnc_df[hgnc_df[col] == gene_symbol]
        if df.shape[0] == 0:
            continue
        elif df.shape[0] == 1:
            hgnc_symbol = df.loc[df.index[0], 'Approved Symbol']
            hgnc_id = df.loc[df.index[0], 'HGNC ID']
            found = True
            break
        else:
            raise Exception('Non unique match for gene_symbol = {} and col = {}: shape = {}'.format(
                gene_symbol, col, df.shape
            ))

    if not found:
        logger.warning('HGNC match is not found for Symbol={}'.format(gene_symbol))

    logger.debug('gene_symbol = {}, hgnc_symbol = {}, hgnc_id = {}'.format(
        gene_symbol, hgnc_symbol, hgnc_id
    ))
    if gene_symbol != hgnc_symbol:
        logger.debug('gene_symbol {} was changed to {}!'.format(gene_symbol, hgnc_symbol))

    return hgnc_id, hgnc_symbol


def get_symbol_from_entrez_id(entrez_id):
    return 'foo_entrez_hgnc'