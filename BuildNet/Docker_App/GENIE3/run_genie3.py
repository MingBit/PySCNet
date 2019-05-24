from __future__ import absolute_import
import numpy as np
from BuildNet.pylib import genie3_raw as genie3
import pandas as pd


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
