#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 11:59:54 2019

@author: mwu
"""


import scanpy.api as sc
import pandas as pd
import numpy as np
import anndata
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import r
import rpy2.robjects as ro



path = '/home/mwu/MING_V9T/PhD_Pro/Test/SS_Run30_D3-5/'

#scanpy pipeline
design = pd.read_csv(path + 'Seurat_Data/Cell_Info.tsv', delimiter = ' ', index_col = 0)
expr = pd.read_csv(path + 'Seurat_Data/Scale_Data.tsv', delimiter=' ', index_col = 0)



def _create_scanpy(expr, feature, result_file, **kwargs):
        """ scanpy pipeline for QC, clustering, marker gene identification
        """
        adata = anndata.AnnData(expr.T)
        adata.var = feature

        #filtering
        sc.pp.filter_cells(adata, **kwargs)
        sc.pp.filter_genes(adata, **kwargs)

        mito_genes = adata.var_names.str.startswith('MT-')
        adata.obs['percent_mito'] = np.sum(
                        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1

        adata.obs['n_counts'] = adata.X.sum(axis=1).A1
        adata = adata[adata.obs['n_genes'] < 3000, :]
        adata = adata[adata.obs['percent_mito'] < 0.05, :]

        #normalize
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)

        #gene selection
        sc.pp.highly_variable_genes(adata,
                              min_mean = 0.025, max_mean = 2, min_disp = 0.4)

        adata_filter = adata [:,adata.var['highly_variable']]
        sc.pp.scale(adata_filter, max_value=10)

        #pca and clustering
        sc.tl.pca(adata)
        sc.pp.neighbors(adata_filter, n_neighbors=50, n_pcs=20)
        sc.tl.louvain(adata_filter)
        sc.tl.tsne(adata_filter)
        sc.tl.umap(adata_filter)

        #find marker gene
        sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
        adata.write(result_file)

        return(adata)


def _create_seurat(expr, design, result_file, **kwargs):
        """ seurat pipeline for QC, clustering and marker gene detection
        """
       #install & import required R packages
       utils = rpackages.importr('utils')
       utils.chooseCRANmirror(ind = 16)
       package_names = ('M3Drop', 'scater', 'tibble', 'ggpubr', 'dplyr', 'data.table',
                        'SingleCellExperiment', 'openxlsx', 'Seurat', 'base', 'Matrix')

       names_to_install = [x for x in package_names if not rpackages.isinstalled(x)]

       if len(names_to_install) > 0:
               utils.install_packages(StrVector(names_to_install))

       design = ro.pandas2ri.py2rpy.

       with localconverter(ro.default_converter + pandas2ri.converter):
               design = ro.conversion.py2rpy(design)

       #Seurat pipeline
       r('''
             print(dim(tmp))
             New_sc = SingleCellExperiment(assays = list(counts = as.matrix(as.data.frame(tmp)),
                                                         colData = as.data.frame(design))) %>%
                                                         calculateQCMetrics()

         ''')



     cell_counts <- as.data.frame(New_sce@colData)
     select_cells <- rownames(cell_counts[which(cell_counts$total_counts > 5000 & cell_counts$total_features > 500),])
     Filter_SC <- Expr_Df[, colnames(Expr_Df) %in% select_cells]
