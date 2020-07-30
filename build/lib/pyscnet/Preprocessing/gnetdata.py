#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:53:34 2019

@author: mwu
"""

import _pickle as pk
import copy
import pandas as pd


class Gnetdata:
    """
    Gnetdata class includes dataset specifically
    ------------------------------------------------
    for building GRNs. It consists four sub-classes:
    1) ExpMatrix: A gene count matrix
    2) CellAttrs: A dict gives information about cells. eg. cluster_nr, annotation
    3) GeneAttrs: A dict gives information about genes. eg. module_nr, marker_annotation
    4) NetAttrs: A dict includes Networks attributes. eg. node centralities.
    """

    def __init__(self, ExpMatrix, CellAttrs=None, GeneAttrs=None):
        self.ExpMatrix = ExpMatrix
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
        """return the shape of ExpMatrix"""
        return self.ExpMatrix.shape

    @property
    def deepcopy(self):
        """make a deepcopy of gnetData """
        return copy.deepcopy(self)

    @property
    def info(self):
        """ return the shape of ExpMatrix and the keys of CellAttrs, GeneAttrs, NetAttrs"""
        text = 'Gnetdata object with \nExpMatrix: {} x {}'.format(self.shape[0], self.shape[1])
        attr = ['CellAttrs', 'GeneAttrs', 'NetAttrs']
        for at in attr:
            keys = getattr(self, at).keys()
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


def load_from_scanpy(scanpy_obj=None, filepath=None):
    """
    load adata object
    --------------------------------------------
    :param scanpy_obj: adata object, default None
    :param filepath: str, default None. filepath refers to adata object
    :return: Gnetdata
    """
    if scanpy_obj is not None:
        gnetdata = Gnetdata(pd.DataFrame(scanpy_obj.X.T, index=scanpy_obj.var_names, columns=scanpy_obj.obs_names),
                            CellAttrs=dict({'CellInfo': scanpy_obj.obs}),
                            GeneAttrs=dict({'GeneInfo': scanpy_obj.var}))

    elif filepath is not None:
        with open(filepath, 'rb') as input:
            scanpy_obj = pk.load(input)
            gnetdata = Gnetdata(pd.DataFrame(scanpy_obj.X.T, index=scanpy_obj.var_names, columns=scanpy_obj.obs_names),
                                CellAttrs=dict({'CellInfo': scanpy_obj.obs}),
                                GeneAttrs=dict({'GeneInfo': scanpy_obj.var}))

    return gnetdata
