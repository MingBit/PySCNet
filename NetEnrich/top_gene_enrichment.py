#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:50:03 2019

@author: mwu
"""

import gseapy
from scipy import stats as st
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage
from NetEnrich.signed_ks_test_snippets import FUNCTIONS
rsnippets = SignatureTranslatedAnonymousPackage(FUNCTIONS, 'rsnippets')



def gene_compare(gene_list, gene_set):
        """ it returns enrichment scores for given two ordered gene lists"""
        S = list(range(1, len(gene_set) + 1))
        overlap = [gene_set.index(x) + 1 for x in set(gene_list).intersection(set(gene_set))]
        S = [i for j, i in enumerate(S) if j not in overlap]
        return rsnippets.ks_test_2(overlap, S, maxCombSize = 10^10)
