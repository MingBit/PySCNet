from __future__ import absolute_import
import sys
sys.path.append('/Users/angelawu/GitHub/SCNetEnrich/BuildNet')

from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage
from rlib.snippets import FUNCTIONS

rsnippets = SignatureTranslatedAnonymousPackage(FUNCTIONS, 'rsnippets')

def base():
    return importr('base')


def wgcna():
    return importr('WGCNA')


def stats():
    return importr('stats')


def graphics():
    return importr('graphics')


def grdevices():
    return importr('grDevices')


def pheatmap():
    return importr('pheatmap')

def dynamicTreeCut():
    return importr('dynamicTreeCut')