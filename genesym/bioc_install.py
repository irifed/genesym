#!/usr/bin/env python3

from rpy2.robjects.packages import importr

# do the following _only the first time_, to install the package seqLogo
base = importr('base')
base.source('http://www.bioconductor.org/biocLite.R')
biocinstaller = importr('BiocInstaller')
biocinstaller.biocLite('seqLogo')

# load the installed package 'seqLogo'
seqlogo = importr('seqLogo')