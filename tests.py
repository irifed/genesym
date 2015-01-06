import unittest

from genesym import geodriver

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
        self.assertEqual(gpldf.shape[0], nrows)

        self.assertEqual(gpldf['Gene Symbol'][54400-1], 'WRAP53')

if __name__ == '__main__':
    unittest.main()