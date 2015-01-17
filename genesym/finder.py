import time
import pandas
import logging


logger = logging.getLogger('genesym')


class FileFinder:

    def __init__(self, data_fn):
        self.data_fn = data_fn
        logger.debug('Loading lookup file {}...'.format(self.data_fn))
        tic = time.time()
        self.df = pandas.read_csv(self.data_fn, delimiter='\t')
        toc = time.time()
        logger.debug('File {} load took {} sec'.format(self.data_fn, toc - tic))

    def lookup_fast(self, value, cols):
        """
        Looks for exact match in given columns.
        Assume for now that given value may occur in only one of columns.
        TODO implement case when value may be found in different columns?

        :param value: exact value which is being searched for
        :param cols: list of dataframe columns where search is performed
        :return: dataframe with first found result, or empty df if nothing found
        """

        for col in cols:
            # exact match condition
            res_df = self.df[self.df[col] == value]
            if res_df.shape[0] == 0:
                # will look in next column
                continue
            elif res_df.shape[0] >= 1:
                if res_df.shape[0] > 1:
                    logger.warning('Non unique match for '
                                   'value = {} and col = {}: '
                                   'shape = {}'.format(value, col, res_df.shape)
                    )
                return res_df

        logger.warning('Exact match is not found for '
                       'value={} and columns={}'.format(value, cols)
        )
        return res_df

    def lookup_slow(self, value, cols):
        """
        Looks for substring match in given columns.
        Assume for now that given value may occur in only one of columns.
        TODO implement case when value may be found in different columns?

        :param value: exact value which is being searched for
        :param cols: list of dataframe columns where search is performed
        :return: dataframe with first found result, or empty df if nothing found
        """
        for col in cols:
            # whole word substring match condition
            res_df = self.df[
                self.df[col].fillna('').str.contains(
                    r'\b{}\b'.format(value)
                )
            ]
            if res_df.shape[0] == 0:
                # will look in next column
                continue
            elif res_df.shape[0] >= 1:
                if res_df.shape[0] > 1:
                    logger.warning('Non unique match for '
                                   'value = {} and col = {}: '
                                   'shape = {}'.format(value, col, res_df.shape)
                    )
                return res_df

        logger.warning('Whole word substring match is not found for '
                       'value={} and columns={}'.format(value, cols)
        )
        return res_df


