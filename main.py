#!/usr/bin/env python3

import pandas
import time
import csv
from genesym.geodriver import get_gpl, get_hgnc_id_symbol

import logging

logger = logging.getLogger('genesym')


def process_platform(gpl_id):
    logfh = logging.FileHandler(gpl_id + '.log')
    logger.addHandler(logfh)

    gpldf = get_gpl(gpl_id)
    print(gpldf.columns)

    # create new column in dataframe
    gpldf['HGNC_Symbol'] = pandas.Series(index=gpldf.index)
    gpldf['HGNC_ID'] = pandas.Series(index=gpldf.index)

    # DEBUG
    gpldf = gpldf[215:220]

    n_found_hgnc_ids = 0
    for idx, row in gpldf.iterrows():
        logger.info('Processing ID={}'.format(row['ID']))

        hgnc_id, hgnc_symbol = get_hgnc_id_symbol(row)
        if hgnc_id is not None:
            gpldf.loc[idx, 'HGNC_ID'] = hgnc_id
            n_found_hgnc_ids += 1
        if hgnc_symbol is not None:
            gpldf.loc[idx, 'HGNC_Symbol'] = hgnc_symbol

        logging.debug('hgnc_id = {}, hgnc_symbol = {}'.format(hgnc_id, hgnc_symbol))

    # print some stats
    print('Total number of rows: {}'.format(gpldf.shape[0]))
    print('Number or HGNC IDs found: {} ({}%)'.format(
        n_found_hgnc_ids, 100*n_found_hgnc_ids/gpldf.shape[0]))

    # DEBUG
    for idx, row in gpldf.iterrows():
        print(row)

    gpldf.to_csv(gpl_id + '.hgnc.txt', sep='\t', quoting=csv.QUOTE_NONNUMERIC)
    gpldf[['ID', 'HGNC_ID', 'HGNC_Symbol']].to_csv(gpl_id + '.only_hgnc.txt', sep='\t')


def main():
    logger.info('Starting processing data')
    # TODO load list of platforms
    tic = time.time()
    process_platform('GPL570')
    toc = time.time()

    logger.debug('Total time = {}'.format(toc - tic))

    logger.info('All done!')


if __name__ == '__main__':
    main()