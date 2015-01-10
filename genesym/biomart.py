import time
import pprint
import logging

import rpy2.robjects as robjects
import pandas.rpy.common as com

logger = logging.getLogger('genesym')


class BioMart:

    def __init__(self):
        logging.debug('Starting BioMart init...')
        tic = time.time()
        robjects.r('''
            suppressMessages(library(biomaRt))
            ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
        ''')

        toc = time.time()
        logging.debug('Starting BioMart done in {} sec'.format(toc - tic))

    def lookup_unigene_id(self, unigene_id):
        # return None if no match found
        hgnc_id = None
        hgnc_symbol = None

        tic = time.time()

        robjects.r('''
            hgnc_id <- getBM("hgnc_id", filters = "unigene", values = "{}", mart = ensembl)
        '''.format(unigene_id))
        hgnc_id_df = com.load_data('hgnc_id')
        if hgnc_id_df.shape == (1, 1):
            hgnc_id = hgnc_id_df.iloc[0, 0]
        else:
            logger.warning('Non-unique match found for Unigene ID = {}: hgnc_id = {}'.format(
                unigene_id, hgnc_id_df
            ))

        robjects.r('''
            hgnc_symbol <- getBM("hgnc_symbol", filters = "unigene", values = "{}", mart = ensembl)
        '''.format(unigene_id))
        hgnc_symbol_df = com.load_data('hgnc_symbol')
        if hgnc_symbol_df.shape == (1, 1):
            hgnc_symbol = hgnc_symbol_df.iloc[0, 0]
        else:
            logger.warning('Non-unique match found for Unigene ID = {}: hgnc_symbol = {}'.format(
                unigene_id, hgnc_symbol_df
            ))

        toc = time.time()
        logger.debug('BioMart Unigene ID lookup took {} sec'.format(toc-tic))
        logger.debug('BioMart Unigene ID lookup found ' +
                     'HGNC ID = {}, HGNC Symbol = {} for Unigene_ID={}'.format(
                         hgnc_id, hgnc_symbol, unigene_id
                     ))

        return hgnc_id, hgnc_symbol

biomart = BioMart()
