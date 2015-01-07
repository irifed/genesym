import time
import pandas
import logging
import pprint

logger = logging.getLogger('genesym')

class HGNCLookup:

    def __init__(self):
        logger.debug('Loading HGNC lookup file...')

        hgnc_complete_set_fn = '/Users/irina/Projects/insilico medicine/datasets/hgnc_complete_set.txt'

        tic = time.time()
        self.hgnc_df = pandas.read_csv(hgnc_complete_set_fn, delimiter='\t')
        toc = time.time()
        logger.debug('HGNC file load took {} sec'.format(toc - tic))

        self.cols = ['Approved Symbol', 'Previous Symbols', 'Synonyms', 'HGNC ID']
        self.hgnc_df = self.hgnc_df[self.cols]


    def lookup(self, gene_symbol):
        hgnc_symbol = None
        hgnc_id = str(gene_symbol) + '_HGNC'

        found = False

        for col in self.cols[0:3]:
            df = self.hgnc_df[self.hgnc_df[col] == gene_symbol]
            if df.shape[0] == 0:
                continue
            elif df.shape[0] == 1:
                hgnc_symbol = df.loc[df.index[0], 'Approved Symbol']
                hgnc_id = df.loc[df.index[0], 'HGNC ID']
                found = True
                break
            else:
                logger.warning('Non unique match for gene_symbol = {} and col = {}: shape = {}'.format(
                    gene_symbol, col, df.shape
                ))
                pass

        if not found:
            logger.warning('HGNC match is not found for Symbol={}'.format(gene_symbol))
            hgnc_id, hgnc_symbol = self.lookup_slow(gene_symbol)

        logger.debug('gene_symbol = {}, hgnc_symbol = {}, hgnc_id = {}'.format(
            gene_symbol, hgnc_symbol, hgnc_id
        ))
        if gene_symbol != hgnc_symbol:
            logger.debug('gene_symbol {} was changed to {}!'.format(gene_symbol, hgnc_symbol))

        return hgnc_id, hgnc_symbol

    def lookup_slow(self, gene_symbol):
        hgnc_symbol = None
        hgnc_id = str(gene_symbol) + '_HGNC'

        found = False
        for col in self.cols[1:3]:
            df = self.hgnc_df[self.hgnc_df[col].fillna('').str.contains(r'\b{}\b'.format(gene_symbol))]
            if df.shape[0] == 0:
                continue
            elif df.shape[0] == 1:
                hgnc_symbol = df.loc[df.index[0], 'Approved Symbol']
                hgnc_id = df.loc[df.index[0], 'HGNC ID']
                found = True
                break
            else:
                logger.debug('non-unique match: {}'.format(pprint.pprint(df)))
                logger.warning('Non unique match for gene_symbol = {} and col = {}: shape = {}'.format(
                    gene_symbol, col, df.shape
                ))
                pass

        if not found:
            logger.warning('HGNC match is not found for Symbol={}'.format(gene_symbol))

        return hgnc_id, hgnc_symbol
