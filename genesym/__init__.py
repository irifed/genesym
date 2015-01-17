import logging

logger = logging.getLogger('genesym')

# def setup_logging():
fmt = '%(asctime)s: %(levelname)s: %(filename)s: ' + \
      '%(funcName)s(): %(lineno)d: %(message)s'
logging.basicConfig(
    format=fmt,
    level=logging.DEBUG)

logger.setLevel(level=logging.DEBUG)

logfh = logging.FileHandler('genesym.log')
logfh.setFormatter(logging.Formatter(fmt))
logger.addHandler(logfh)

from .hgncfile import HGNCFinder
hgnc = HGNCFinder()

# from .biomartweb import BioMartWeb
# biomartweb = BioMartWeb()

from .biomartfile import BioMartFinder
biomartfile = BioMartFinder()

