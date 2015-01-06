#!/usr/bin/env python3

import pandas
import time
from genesym.geodriver import get_gpl, get_hgnc_id_symbol
from pprint import pprint

import logging

logger = logging.getLogger('genesym')


def setup_logging():
    format = '%(asctime)s: %(levelname)s: %(filename)s: %(funcName)s(): %(message)s'
    logging.basicConfig(
        format=format,
        level=logging.DEBUG)

    logger.setLevel(level=logging.DEBUG)

    logfh = logging.FileHandler('genesym.log')
    logfh.setFormatter(logging.Formatter(format))
    logger.addHandler(logfh)


def process_platform(gpl_id):
    gpldf = get_gpl(gpl_id)
    print(gpldf.columns)

    # create new column in dataframe
    gpldf['HGNC_Symbol'] = pandas.Series(index=gpldf.index)
    gpldf['HGNC_ID'] = pandas.Series(index=gpldf.index)

    # DEBUG
    gpldf = gpldf[215:220]


    for idx, row in gpldf.iterrows():
        logger.info('Processing ID={}'.format(row['ID']))

        hgnc_id, hgnc_symbol = get_hgnc_id_symbol(row)
        gpldf.loc[idx, 'HGNC_ID'] = hgnc_id
        gpldf.loc[idx, 'HGNC_Symbol'] = hgnc_symbol

        logging.debug('hgnc_id = {}, hgnc_symbol = {}'.format(hgnc_id, hgnc_symbol))

    # DEBUG
    for idx, row in gpldf.iterrows():
        print(row)

    result_fn = gpl_id + '.hgnc.txt'
    gpldf.to_csv(result_fn, sep='\t')


def main():
    logger.info('Starting processing data')
    # TODO load list of platforms
    tic = time.time()
    process_platform('GPL570')
    toc = time.time()

    logger.debug('Total time = {}'.format(toc - tic))

    logger.info('All done!')


if __name__ == '__main__':
    setup_logging()
    main()