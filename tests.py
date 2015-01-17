import unittest
import pandas
import numpy

from genesym import geodriver
from genesym import hgnc

from genesym import biomartfile as biomart


class MyTestCase(unittest.TestCase):

    def test_get_gpl(self):
        gpldf = geodriver.get_gpl('GPL570')

        columns = ['ID', 'GB_ACC', 'SPOT_ID', 'Species Scientific Name',
                   'Annotation Date', 'Sequence Type', 'Sequence Source',
                   'Target Description', 'Representative Public ID',
                   'Gene Title', 'Gene Symbol', 'ENTREZ_GENE_ID',
                   'RefSeq Transcript ID', 'Gene Ontology Biological Process',
                   'Gene Ontology Cellular Component',
                   'Gene Ontology Molecular Function']
        self.assertEqual(
            list(gpldf.columns),
            columns)

        nrows = 54675
        self.assertEqual(nrows, gpldf.shape[0])

        self.assertEqual('WRAP53', gpldf['Gene Symbol'][54400-1])

    def test_hgnc_lookup_1(self):
        hgnc_id, hgnc_sym = hgnc.lookup_gene_symbol('BTBD16')
        self.assertEqual('HGNC:26340', hgnc_id)
        self.assertEqual('BTBD16', hgnc_sym)

    def test_hgnc_lookup_2(self):
        hgnc_id, hgnc_sym = hgnc.lookup_gene_symbol('ACTGP5')
        self.assertEqual('HGNC:150', hgnc_id)
        self.assertEqual('ACTBP12', hgnc_sym)

    def test_hgnc_lookup_3(self):
        hgnc_id, hgnc_sym = hgnc.lookup_gene_symbol('XPMC2H')
        self.assertEqual('HGNC:12820', hgnc_id)
        self.assertEqual('REXO4', hgnc_sym)

    def test_hgnc_lookup_4(self):
        hgnc_id, hgnc_sym = hgnc.lookup_gene_symbol('FAM80B2')
        self.assertEqual('HGNC:34034', hgnc_id)
        self.assertEqual('RIMKLBP1', hgnc_sym)

    def test_hgnc_lookup_5(self):
        hgnc_id, hgnc_sym = hgnc.lookup_gene_symbol('TOB1')
        self.assertEqual('HGNC:11979', hgnc_id)
        self.assertEqual('TOB1', hgnc_sym)

    def test_hgnc_lookup_6(self):
        hgnc_id, hgnc_sym = hgnc.lookup_gene_symbol('ABC8')
        self.assertEqual('HGNC:73', hgnc_id)
        self.assertEqual('ABCG1', hgnc_sym)

    def test_hgnc_lookup_7(self):
        hgnc_id, hgnc_sym = hgnc.lookup_gene_symbol('L29')
        self.assertEqual('HGNC:10331', hgnc_id)
        self.assertEqual('RPL29', hgnc_sym)

    def test_get_hgnc_id_symbol_1(self):
        # trivial test: Gene Symbol equals to HGNC
        row = pandas.Series(data=pandas.Series(
            {'ID': '1007_s_at',
             'GB_ACC': 'U48705',
             'Gene Symbol': 'DDR1',
             'ENTREZ_GENE_ID': '780'}
        ))

        hgnc_id, hgnc_symbol = geodriver.get_hgnc_id_symbol(row)
        self.assertEqual('DDR1', hgnc_symbol)
        self.assertEqual('HGNC:2730', hgnc_id)

    def test_get_hgnc_id_symbol_2(self):
        # empty Symbol and Entrez Gene ID, but Unigene is present
        # example from GPL10558
        row = pandas.Series(data=pandas.Series(
            {'ID': 'ILMN_1913641',
             'GB_ACC': 'AL833463',
             'Symbol': numpy.nan,
             'Entrez_Gene_ID': numpy.nan,
             'Unigene_ID': 'Hs.87194'}
        ))

    # TODO test: non-empty symbol, but is outdated (hgnc approved is different)
    # TODO test: empty symbol, but entrez gene id is present

        hgnc_id, hgnc_symbol = geodriver.get_hgnc_id_symbol(row)
        self.assertEqual('PYGO1', hgnc_symbol)
        self.assertEqual('HGNC:30256', hgnc_id)

    def test_biomart_unigene_lookup_1(self):
        # normal case
        hgnc_id, hgnc_symbol = biomart.lookup_unigene_id('Hs.447469')

        self.assertEqual('HGNC:15210', hgnc_id)
        self.assertEqual('OR52B5P', hgnc_symbol)

    def test_biomart_unigene_lookup_2(self):
        # non-existent unigene id
        bad_id = 'Hs.000'
        hgnc_id, hgnc_symbol = biomart.lookup_unigene_id(bad_id)

        self.assertEqual(None, hgnc_id)
        self.assertEqual(None, hgnc_symbol)


if __name__ == '__main__':
    unittest.main()
