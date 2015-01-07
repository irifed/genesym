import logging

logger = logging.getLogger('genesym')

# def setup_logging():
format = '%(asctime)s: %(levelname)s: %(filename)s: ' + \
         '%(funcName)s(): %(lineno)d: %(message)s'
logging.basicConfig(
    format=format,
    level=logging.DEBUG)

logger.setLevel(level=logging.DEBUG)

logfh = logging.FileHandler('genesym.log')
logfh.setFormatter(logging.Formatter(format))
logger.addHandler(logfh)
