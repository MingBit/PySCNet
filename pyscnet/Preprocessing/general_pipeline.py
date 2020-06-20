#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 11:59:54 2019
@author: mwu
"""
from __future__ import absolute_import
import _create_scanpy as scanpy

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri


def run_pipeline(gnetdata, pipeline='seurat'):
    """for given gnetdata, run scanpy/seurat pipeline and return cell clustering results
        """
    Expr = gnetdata.ExpMatrix
    Design = gnetdata.CellAttrs['Design']

    # TODO: check Design is a dataFrame

    if pipeline == 'seurat':
        print('Debug: Seurat')
    #                r_seurat = seurat.seurat_pipeline(Expr, Design)
    #                seurat_obj = ro.conversion.rpy2py(r_seurat)
    #                with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
    #                        gnetdata.CellAttrs['Design'] = ro.conversion.rpy2py(seurat_obj[0])
    #                gnetdata._add_cellattr('seurat_obj', seurat_obj[1])

    elif pipeline == 'scanpy':
        scanpy_obj = scanpy.scanpy_pipeline(Expr)
        scanpy_obj.obs.colunms = ['n_genes', 'n_counts', 'cluster_id']
        gnetdata.CellAttrs['Design'] = scanpy_obj.obs
        gnetdata._add_cellattr('scanpy_obj', scanpy_obj)

    else:
        raise Exception('Two methods are valid: seurat, scanpy')

    return gnetdata
