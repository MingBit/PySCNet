from __future__ import absolute_import
from arboreto.algo import grnboost2
import numpy as np
import pandas as pd


def run_grnboost2(Expr, filename='links.txt', gene_names=None, **kwargs):
    links = grnboost2(np.asmatrix(Expr.T), gene_names=gene_names, **kwargs)
    links.to_csv(filename, sep='\t', index=False, header=False)


if __name__ == '__main__':
    Expr = pd.read_csv('Expr.txt', sep='\t', header=0, index_col=0)
    run_grnboost2(Expr, gene_names=Expr.index)
