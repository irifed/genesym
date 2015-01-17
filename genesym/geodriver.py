import time

import rpy2.robjects as robjects
import pandas.rpy.common as com
import pandas
import numpy
import config

import logging
logger = logging.getLogger('genesym')

from genesym import hgnc, biomartfile

def get_gpl_r(gpl_id):
    tic = time.time()
    robjects.r('''
        suppressMessages(library(GEOquery))
        destdir <- "{}"
        gpl <- getGEO("{}", destdir = destdir)
        gpldf <- Table(gpl)
    '''.format(config.data_dir, gpl_id))
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
    soft_fn = config.data_dir + gpl_id + '.soft'
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
    """Main routine for finding HGNC id and symbol for row from GPL dataframe

    :param row: row from GPL dataframe
    :return: hgnc_id, hgnc_symbol for given row, None values if not found
    """
    tic = time.time()
    # get gene symbol value
    # if it does not exist, get entrez id and lookup in biomart
    # if it does not exist, get unigene id and lookup in biomart
    # ... genbank accession
    # ...
    hgnc_id, hgnc_symbol = None, None

    # TODO analyze row keys() to find attribute names of Symbol, Unigene, etc
    symbol_keynames = ['Symbol', 'Gene Symbol']
    unigene_keyname = ['Unigene_ID']

    gene_symbol = get_attribute(row, symbol_keynames)
    if gene_symbol is not None:
        hgnc_id, hgnc_symbol = hgnc.lookup_gene_symbol(gene_symbol)

    else:
        unigene_id = get_attribute(row, unigene_keyname)
        if unigene_id is not None:
            hgnc_id, hgnc_symbol = biomartfile.lookup_unigene_id(unigene_id)

        # TODO: lookup other ids
        # entrez
        # unigene
        # genbank id GI
        # genbank accession

    toc = time.time()

    logging.debug('lookup time = {}'.format(toc - tic))

    return hgnc_id, hgnc_symbol


def get_attribute(row, attr_list=list()):
    # find attribute in row, return first found attribute from attr_list

    result = None

    for attr in attr_list:
        if attr in row:
            result = row[attr]
        else:
            logger.debug('{} attribute is not found in row {}'.format(
                attr, row.keys()
            ))
        if result is not None:
            # protect from NA values
            if type(result) == float and numpy.isnan(result):
                result = None
            break

    return result
