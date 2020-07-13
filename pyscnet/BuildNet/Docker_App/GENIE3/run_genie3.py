from __future__ import absolute_import
import os.path
from arboreto.algo import genie3
import numpy as np
import pandas as pd


def run_genie3(Expr, filename='links.txt', gene_names=None, **kwargs):
    links = genie3(np.asmatrix(Expr.T), gene_names=gene_names, **kwargs)
    links.to_csv(filename, sep='\t', index=False, header=False)


if __name__ == '__main__':
    tf_names = list(pd.read_csv('TF_Names.txt', sep='\t', header=0, index_col=0).index) if os.path.isfile('TF_Names.txt') else None
    Expr = pd.read_csv('Expr.txt', sep='\t', header=0, index_col=0)
    run_genie3(Expr, gene_names=list(Expr.index), tf_names = tf_names)
