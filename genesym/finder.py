import time
import pandas
import logging
import pprint


logger = logging.getLogger('genesym')


class FileFinder:

    def __init__(self, data_fn,
                 hgnc_symbol_colname=None,
                 hgnc_id_colname=None):
        logger.debug('Loading lookup file {}...'.format(data_fn))

        tic = time.time()
        self.hgnc_df = pandas.read_csv(data_fn, delimiter='\t')
        toc = time.time()
        logger.debug('File load took {} sec'.format(toc - tic))

        self.hgnc_symbol_colname = hgnc_symbol_colname
        self.hgnc_id_colname = hgnc_id_colname
        # TODO verify that these colnames exist

    def get_hgnc_symbol(self, df):
        if self.hgnc_symbol_colname in df.keys():
            return df.loc[df.index[0], self.hgnc_symbol_colname]
        return None

    def get_hgnc_id(self, df):
        if self.hgnc_id_colname in df.keys():
            return df.loc[df.index[0], self.hgnc_id_colname]
        return None

    def lookup(self, gene_symbol):
        return self.lookup_fast_low(
            gene_symbol,
            ['Approved Symbol', 'Previous Symbols', 'Synonyms']
        )

    def lookup_fast_low(self, gene_symbol, cols):
        # looks for exact match in given columns
        hgnc_symbol = None
        hgnc_id = None

        found = False

        for col in cols:
            df = self.hgnc_df[self.hgnc_df[col] == gene_symbol]
            if df.shape[0] == 0:
                continue
            elif df.shape[0] == 1:
                hgnc_symbol = self.get_hgnc_symbol(df)
                hgnc_id = self.get_hgnc_id(df)
                found = True
                break
            else:
                logger.warning(
                    'Non unique match for gene_symbol = {} and col = {}: '
                    'shape = {}'.format(gene_symbol, col, df.shape
                ))
                # skip symbol this for now,
                # this symbol will be verified using other available id's
                pass

        if not found:
            hgnc_id, hgnc_symbol = self.lookup_slow(gene_symbol)

        if hgnc_id is None or hgnc_symbol is None:
            logger.warning('HGNC match is not found for Symbol={}'.format(
                gene_symbol)
            )
        else:
            logger.debug('gene_symbol = {}, '
                         'hgnc_symbol = {}, hgnc_id = {}'.format(
                gene_symbol, hgnc_symbol, hgnc_id
            ))
            if gene_symbol != hgnc_symbol:
                logger.debug('gene_symbol {} was changed to {}!'.format(
                    gene_symbol, hgnc_symbol)
                )

        return hgnc_id, hgnc_symbol

    def lookup_slow(self, gene_symbol):
        return self.lookup_slow_low(
            gene_symbol,
            ['Previous Symbols', 'Synonyms']
        )

    def lookup_slow_low(self, gene_symbol, cols):
        # looks for substring (but whole word) match
        hgnc_symbol = None
        hgnc_id = None

        found = False
        for col in cols:

            df = self.hgnc_df[
                self.hgnc_df[col].fillna('').str.contains(
                    r'\b{}\b'.format(gene_symbol)
                )
            ]
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
                # skip this for now,
                # hope this symbol will be found using other methods
                pass

        if not found:
            logger.warning('HGNC match is not found for Symbol={}'.format(
                gene_symbol
            ))

        return hgnc_id, hgnc_symbol

    def lookup_entrez_id(self, entrez_id):
        return self.lookup_fast_low(entrez_id, ['Entrez_Gene_ID'])

    def lookup_unigene_id(self, unigene_id):
        # unigene lookup file has only hgnc id column, but no symbol
        hgnc_id, hgnc_symbol = self.lookup_fast_low(unigene_id, ['Unigene ID'])
        return hgnc_id, hgnc_symbol

