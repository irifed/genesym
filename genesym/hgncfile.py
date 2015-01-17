import logging
logger = logging.getLogger('genesym')

import config

from .finder import FileFinder


class HGNCFinder(FileFinder):
    def __init__(self):

        super().__init__(config.hgnc_complete_set_fn)

        self.hgnc_symbol_colname = 'Approved Symbol'
        self.hgnc_id_colname = 'HGNC ID'

    def get_hgnc_symbol(self, row):
        if self.hgnc_symbol_colname in row.keys():
            return row[self.hgnc_symbol_colname]
        return None

    def get_hgnc_id(self, row):
        if self.hgnc_id_colname in row.keys():
            return row[self.hgnc_id_colname]
        return None

    def lookup_gene_symbol(self, gene_symbol):
        hgnc_id = None
        hgnc_symbol = None

        df = self.lookup_fast(gene_symbol,
            ['Approved Symbol', 'Previous Symbols', 'Synonyms']
        )
        if df.shape[0] == 0:
            df = self.lookup_slow(gene_symbol, ['Previous Symbols', 'Synonyms'])

        if df.shape[0] >= 1:
            # return fist match
            # TODO smarter select in non-unique case
            row = df.iloc[0]
            hgnc_id = self.get_hgnc_id(row)
            hgnc_symbol = self.get_hgnc_symbol(row)
        else:
            logger.warning('Gene symbol {} is not found in file {}'.format(
                gene_symbol, self.data_fn)
            )

        return hgnc_id, hgnc_symbol

    def lookup_entrez_id(self, entrez_id):
        hgnc_id = None
        hgnc_symbol = None

        df = self.lookup_fast(entrez_id, ['Entrez_Gene_ID'])

        if df.shape[0] >= 1:
            # return fist match
            # TODO smarter select in non-unique case
            row = df.iloc[0]
            hgnc_id = self.get_hgnc_id(row)
            hgnc_symbol = self.get_hgnc_symbol(row)
        else:
            logger.warning('Entrez ID {} is not found in file {}'.format(
                entrez_id, self.data_fn)
            )

        return hgnc_id, hgnc_symbol

    def get_symbol_by_id(self, hgnc_id):
        hgnc_symbol = None
        df = self.lookup_fast(hgnc_id, [self.hgnc_id_colname])
        if df.shape[0] >= 1:
            # return fist match
            row = df.iloc[0]
            hgnc_symbol = self.get_hgnc_symbol(row)
        else:
            logger.warning('HGNC ID {} is not found in file {}'.format(
                hgnc_id, self.data_fn)
            )

        return hgnc_symbol
