#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 15:56:02 2019

@author: angelawu
"""

from __future__ import absolute_import
from BuildNet.pylib import bnn_raw as bnn

class Pybnn():
    """python function for bayesian neural network
    """

    def __init__(self, data, params = None):
        self.Expr = data
        self.params = params
        self.BNNExpr = None
        self.GenCor = None

        return None

    def run_bnn(self):

        subpara = {}
        subpara['output_path'] = self.params['output_path']
        subpara['filename'] = self.params['filename']
        subpara['hidden_layers']= [10, 5]

        self.BNNExpr = bnn.get_BNN_table(self.Expr, **subpara)
        links = self.BNNExpr.corr(method = 'spearman').unstack().sort_values().drop_duplicates().reset_index()
        links.columns = ['var1', 'var2','value']
        links = links.loc[(links['var1'] != links['var2'])]
        self.GenCor = links