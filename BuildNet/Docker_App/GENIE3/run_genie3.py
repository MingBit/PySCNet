from __future__ import absolute_import
import genie3_raw as genie3
import numpy as np
import pandas as pd


def run_genie3(Expr, filename = 'links.txt', gene_names = None, **kwargs1):

        VIM = genie3.GENIE3(np.array(Expr), gene_names = gene_names, regulators = gene_names, **kwargs1)
        genie3.get_link_list(VIM, file_name = filename, gene_names = gene_names, regulators = gene_names)
        Links = pd.read_csv(filename, sep= '\t', header = -1)

        return Links

Expr = pd.read_excel('run30_top100_K=6_cluster1_ExprData.xlsx', index_col=0)

run_genie3(Expr)