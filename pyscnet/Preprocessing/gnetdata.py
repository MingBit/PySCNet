#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed May 22 13:53:34 2019

@author: mwu
"""

import _pickle as pk
import copy
import pandas as pd
import scipy

class Gnetdata:
    """
    Gnetdata class includes dataset specifically
    ------------------------------------------------
    for building GRNs. It consists four sub-classes:
    1) Exp: A dict of sparse matrix, var_names and obj_names
    2) CellAttrs: A dict gives information about cells. eg. cluster_nr, annotation
    3) GeneAttrs: A dict gives information about genes. eg. module_nr, marker_annotation
    4) NetAttrs: A dict includes Networks attributes. eg. node centralities.
    """

    def __init__(self, Exp, CellAttrs=None, GeneAttrs=None):
        self.Exp = Exp
        self.CellAttrs = CellAttrs
        self.GeneAttrs = GeneAttrs
        self.NetAttrs = dict()
        self.NetAttrs['parameters'] = dict()

    def _add_cellattr(self, attr_name, attr):

        self.CellAttrs[attr_name] = attr

    def _add_geneattr(self, attr_name, attr):

        self.GeneAttrs[attr_name] = attr

    def _add_netattr(self, attr_name, attr):

        self.NetAttrs[attr_name] = attr

    def _add_netattr_para(self, attr_name, attr):

        self.NetAttrs['parameters'][attr_name] = attr

    @property
    def shape(self):
        """return the shape of Exp"""
        return self.Exp['matrix'].shape

    @property
    def deepcopy(self):
        """make a deepcopy of gnetData """
        return copy.deepcopy(self)

    @property
    def info(self):
        """ return the shape of Exp and the keys of CellAttrs, GeneAttrs, NetAttrs"""
        text = 'Gnetdata object with \nExp: {} x {}'.format(self.shape[0], self.shape[1])
        attr = ['CellAttrs', 'GeneAttrs', 'NetAttrs']
        for at in attr:
            keys = getattr(self, at).keys() if getattr(self, at) is not None else None
            text += '\n' + at + ':' + str(keys)

        print(text)

    def save_as(self, outpath):
        """save as pickle object """
        if outpath is not None:
            with open(outpath, 'wb') as outfile:
                pk.dump(self, outfile)
        else:
            raise Exception('filepath cannot be null!')


def load_Gnetdata_object(filepath):
    """
    load Gnetdata (pickle) from local
    ----------------------------------
    :param filepath: str, default None.
    :return: Gnetdata
    """
    with open(filepath, 'rb') as input:
        gnetdata = pk.load(input)

    return gnetdata


def load_from_anndata(anndata_obj=None):
    """
    load adata object
    --------------------------------------------
    :param anndata_obj: adata object, default None
    :return: Gnetdata
    """
    if anndata_obj is not None:
        if scipy.sparse.issparse(anndata_obj.X):
            gnetdata = Gnetdata(Exp=dict({'matrix': anndata_obj.X, 
                                          'feature': anndata_obj.var_names, 
                                          'cell': anndata_obj.obs_names}),
                                CellAttrs=dict({'CellInfo': anndata_obj.obs}),
                                GeneAttrs=dict({'GeneInfo': anndata_obj.var}))
        else:
            gnetdata = Gnetdata(Exp=dict({'matrix': scipy.sparse.csr_array(anndata_obj.X), 'feature': anndata_obj.var_names, 'cell': anndata_obj.obs_names}),
                                CellAttrs=dict({'CellInfo': anndata_obj.obs}),
                                GeneAttrs=dict({'GeneInfo': anndata_obj.var}))

    return gnetdata
