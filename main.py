#!/usr/bin/env python3

import pandas
from genesym.geodriver import get_gpl, get_hgnc_symbol
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

    for idx, row in gpldf[0:10].iterrows():
        gpldf['HGNC_Symbol'] = get_hgnc_symbol(row)

    for idx, row in gpldf[0:10].iterrows():
        print(row)


def main():
    logger.info('Starting processing data')
    # load list of platforms
    # for each platform
    #   process platform
    process_platform('GPL570')

    logger.info('All done!')


if __name__ == '__main__':
    setup_logging()
    main()