#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 22:06:36 2019

@author: angelawu
"""

from __future__ import absolute_import
from BuildNet.pylib import genie3_raw as genie3
import numpy as np
import pandas as pd

# =============================================================================
# TO do list: call genie3 function
# =============================================================================

class PyGenie3():
    """ python funciton for genie3
    """
    def __init__(self, data, debug=False, params = None):
        self.Expr = data.T
        self.VIM = None
        self.Links = None
        self.params = params

        return None

    def run_genie3(self, filename = 'links.txt', gene_names = None, **kwargs1):

        self.VIM = genie3.GENIE3(np.array(self.Expr), gene_names = gene_names, regulators = gene_names, **kwargs1)
        genie3.get_link_list(self.VIM, file_name = filename, gene_names = gene_names, regulators = gene_names)
        self.Links = pd.read_csv(filename, sep= '\t', header = -1)

