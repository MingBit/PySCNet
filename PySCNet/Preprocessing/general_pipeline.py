#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 11:59:54 2019

@author: mwu
"""

from PySCNet.Preprocessing import _create_scanpy as scanpy
from PySCNet.Preprocessing import _create_seurat as seurat



path = '/home/mwu/MING_V9T/PhD_Pro/Test/SS_Run30_D3-5/'

#scanpy pipeline
design = pd.read_csv(path + 'Seurat_Data/Cell_Info.tsv', delimiter = ' ', index_col = 0)
expr = pd.read_csv(path + 'Seurat_Data/Scale_Data.tsv', delimiter=' ', index_col = 0)


def run_pipeline():




















