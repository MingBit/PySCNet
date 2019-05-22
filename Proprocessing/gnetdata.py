#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:53:34 2019

@author: mwu
"""


import pickle as pk

class Gnetdata():
        """ As cell-oriented analysis can be done by scanpy, this class includes dataset specifically
        for building GRNs. It consists four sub-classes:
                1) GeneMatrix: A matrix for reconsctructing GRNs
                2) CellAttrs: A dataframe gives information about cells. eg. cluster_nr, annotation
                3) GeneAttrs: A dataframe gives information about genes. eg. module_nr, marker_annotation
                4) NetAttrs: A dataframe includes Networks attributes. eg. Node Centralities.
        """
        def __init__(self, GeneMatrix, CellAttrs, GeneAttrs):

                self.GeneMatrix = GeneMatrix
                self.CellAttrs = CellAttrs
                self.GeneAttrs = GeneAttrs
                self.NetAttrs = None

                return None

        def save_Gnetdata_object(self, filepath):
                """ export Gnetdata via pickle protocols"""

                with open(filepath, 'wb') as output:
                        pk.dump(self, output, pk.HIGHEST_PROTOCOL)

                output.close()

        def load_Gnetdata_object(filepath):
                """ import Gnetdata via pickle protocols"""

                with open(filepath, 'rb') as input:
                        self = pk.load(input)

                input.close()

