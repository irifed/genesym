import unittest

from genesym import geodriver, hgnc

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

    def test_hgnc_lookup(self):
        hl = hgnc.HGNCLookup()

        hgnc_id, hgnc_sym = hl.lookup('BTBD16')
        self.assertEqual('HGNC:26340', hgnc_id)
        self.assertEqual('BTBD16', hgnc_sym)

        hgnc_id, hgnc_sym = hl.lookup('ACTGP5')
        self.assertEqual('HGNC:150', hgnc_id)
        self.assertEqual('ACTBP12', hgnc_sym)

        hgnc_id, hgnc_sym = hl.lookup('XPMC2H')
        self.assertEqual('HGNC:12820', hgnc_id)
        self.assertEqual('REXO4', hgnc_sym)

        hgnc_id, hgnc_sym = hl.lookup('FAM80B2')
        self.assertEqual('HGNC:34034', hgnc_id)
        self.assertEqual('RIMKLBP1', hgnc_sym)

        hgnc_id, hgnc_sym = hl.lookup('TOB1')
        self.assertEqual('HGNC:11979', hgnc_id)
        self.assertEqual('TOB1', hgnc_sym)

        hgnc_id, hgnc_sym = hl.lookup('ABC8')
        self.assertEqual('HGNC:73', hgnc_id)
        self.assertEqual('ABCG1', hgnc_sym)

        hgnc_id, hgnc_sym = hl.lookup('L29')
        self.assertEqual('HGNC:10331', hgnc_id)
        self.assertEqual('RPL29', hgnc_sym)


if __name__ == '__main__':
    unittest.main()
