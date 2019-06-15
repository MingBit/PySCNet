from __future__ import absolute_import
import genie3_raw as genie3
import numpy as np
import pandas as pd
import multiprocessing


def run_genie3(Expr, filename = 'links.txt', gene_names = None, **kwargs1):
        numT = multiprocessing.cpu_count()
        VIM = genie3.GENIE3(np.array(Expr), gene_names = gene_names, regulators = gene_names, nthreads = int(numT/2) + 1, **kwargs1)
        genie3.get_link_list(VIM, file_name = filename, gene_names = gene_names, regulators = gene_names)

        return None

Expr = pd.read_csv('Expr.txt', sep = ',', header = 0, index_col = 0).T

run_genie3(Expr, gene_names=list(Expr.columns))