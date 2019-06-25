#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:53:34 2019

@author: mwu
"""


import pickle as pk
from hdfs3 import HDFileSystem
import copy

class Gnetdata():
        """ As cell-oriented analysis can be done by scanpy, this class includes dataset specifically
        for building GRNs. It consists four sub-classes:
                1) GeneMatrix: A matrix for reconsctructing GRNs
                2) CellAttrs: A dict gives information about cells. eg. cluster_nr, annotation
                3) GeneAttrs: A dict gives information about genes. eg. module_nr, marker_annotation
                4) NetAttrs: A dict includes Networks attributes. eg. Node Centralities.
        """
        def __init__(self, GeneMatrix):

                self.GeneMatrix = GeneMatrix
                self.CellAttrs = dict()
                self.GeneAttrs = dict()
                self.NetAttrs = dict()
                self.NetAttrs['parameters'] = dict()

                return None

        def _add_cellattr(self, attr_name, attr):

                self.CellAttrs[attr_name] = attr

        def _add_geneattr(self, attr_name, attr):

                self.GeneAttrs[attr_name] = attr

        def _add_netattr(self, attr_name, attr):

               self.NetAttrs[attr_name] = attr

        def _add_netattr_para(self, attr_name, attr):

                self.NetAttrs['parameters'][attr_name] = attr

        def save_object(self, filepath):
                """ export Gnetdata via pickle protocols"""

                with open(filepath, 'wb') as output:
                        pk.dump(self, output, pk.HIGHEST_PROTOCOL)

                output.close()


def load_Gnetdata_object(gnetdata, filepath):

        hdfs = HDFileSystem(host = 'localhost', port = 8020)

        with hdfs.open(filepath, 'rb') as input:
             gnetdata = pk.load(input)

        input.close()

        return(gnetdata)


def copy_to(data):
        return(copy.deepcopy(data))