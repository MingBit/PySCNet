#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:28:51 2019

@author: mwu
"""#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 14:09:33 2019

@author: mwu
"""

from __future__ import absolute_import
import sys
if('/home/mwu/MING_V9T/PhD_Pro/PySCNet/' not in sys.path): sys.path.append('/home/mwu/MING_V9T/PhD_Pro/PySCNet/')
from pyvis.network import Network
import scanpy.api as sc
import pandas as pd
import numpy as np
import scipy.io
import anndata
import matplotlib.pyplot as plt
from Preprocessing import gnetdata
from BuildNet import gne_dockercaller as gdocker
from NetEnrich import graph_toolkit as gt
from Plotting import dynamic_net as dn




path = '/home/mwu/MING_V9T/PhD_Pro/Test/KK_Run22_D0'

#scanpy pipeline
#cell_type = pd.read_csv('/home/mwu/V9T/LabDZ/Sequencing_old/SeqGeq/full_design.csv', delimiter = ',', index_col = 0)
#run22_run24 = pd.read_csv('/home/mwu/V9T/LabDZ/Sequencing_old/SeqGeq/raw_data_22_24.csv', delimiter=',', index_col = 0)
#run22_Expr = run22_run24[list(cell_type.loc[cell_type.Condition.isin(['Arm'])].index)].T

Expr = scipy.io.mmread('/home/mwu/Downloads/filtered_feature_bc_matrix/matrix.mtx').tocsr()
Feature = pd.read_csv('/home/mwu/Downloads/filtered_feature_bc_matrix/features.tsv', sep = '\t', header = -1)
Feature.columns = ['GeneId', 'GeneName', 'GeneType']

adata = anndata.AnnData(X = Expr.T)
adata.var = Feature
adata.var_names = Feature.GeneName
adata.var_names_make_unique()

#adata = sc.datasets.pbmc3k()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_counts=500)

mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1


sc.pl.scatter(adata, x='n_counts', y='percent_mito')
sc.pl.scatter(adata, x='n_counts', y='n_genes')

adata = adata[adata.obs['n_genes'] < 5000, :]
adata = adata[adata.obs['percent_mito'] < 0.25, :]


sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata,
                              min_mean = 0.025, max_mean = 2, min_disp = 0.4)
sc.pl.highly_variable_genes(adata)
adata_filter = adata [:,adata.var['highly_variable']]
adata_filter.shape

sc.pp.regress_out(adata_filter, ['n_counts', 'percent_mito'])

sc.pp.scale(adata_filter, max_value=10)
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata_filter, n_neighbors=50, n_pcs=20)
sc.tl.umap(adata_filter)
sc.tl.louvain(adata_filter)
sc.tl.tsne(adata_filter)
sc.tl.umap(adata_filter)
sc.pl.umap(adata_filter, color = ['louvain'], save = True)
sc.pl.tsne(adata_filter, color = ['louvain'], save = True)

sc.tl.rank_genes_groups(adata_filter, 'louvain', method='wilcoxon')
sc.pl.rank_genes_groups(adata_filter, n_genes=25, sharey = False, save=True)

marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'CD14',
                'LGALS3', 'GNLY', 'NKG7', 'KLRB1', 'MZB1', 'IL32', 'CD28',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']


tmp = pd.DataFrame(adata_filter.uns['rank_genes_groups']['names']).head(25)


new_cluster_names = [
    'Dendritic', 'Plasmacytoid DC',
    'Naive CD4', 'Mem CD4',
    'NK', 'B', 'CD8',
    'FCGR3A+ Monocytes', 'Megakaryocytes']
adata_filter.rename_categories('louvain', new_cluster_names)

sc.pl.tsne(adata_filter, color='louvain', legend_loc='on data', title='', frameon=False, save='.pdf')

ax = sc.pl.dotplot(adata_filter, marker_genes, groupby='louvain', save='.pdf')


# create gedata object
Mms_fator = pd.read_csv('/home/mwu/MING_V9T/PhD_Pro/PySCNet/Mus_TF_and_TFcofactor/Mus_musculus_TF.txt', sep = '\t')
Expr = pd.DataFrame(adata.X.T, columns=list(adata.obs.index))
Expr.index = [x.lower() for x in list(adata.var_names)]

comm = list(set([x.lower() for x in list(Mms_fator.Symbol)]) & set(Expr.index))
Expr = Expr.loc[comm, :]

run22_gne = gnetdata.Gnetdata(Expr)
run22_gne.CellAttrs['CellInfo'] = adata.obs
run22_gne_PIDC = gdocker.rundocker(run22_gne, method = 'PIDC')
run22_gne_GENIE3 = gdocker.rundocker(run22_gne, method = 'GENIE3')

run22_gne_PIDC = gdocker.buildnet(run22_gne_PIDC, top = 500)
run22_gne_GENIE3 = gdocker.buildnet(run22_gne_GENIE3, top = 500)

run22_gne_GENIE3 = gt.get_centrality(run22_gne_GENIE3)
run22_gne_GENIE3 = gt.community_detect(run22_gne_GENIE3)

run22_gne_PIDC = gt.get_centrality(run22_gne_PIDC)
run22_gne_PIDC = gt.community_detect(run22_gne_PIDC)

dn.dynamic_netShow(run22_gne_PIDC, filepath = '/home/mwu/test_pidc.html')
dn.dynamic_netShow(run22_gne_GENIE3, filepath = '/home/mwu/test_genie3.html')

dn.Genes_modules_Cell_Clusters(run22_gne_GENIE3)

filter_genie3 = run22_gne_GENIE3.NetAttrs['links'].loc[run22_gne_GENIE3.NetAttrs['links'].weight > 0.1]
filter_pidc = run22_gne_PIDC.NetAttrs['links'].loc[run22_gne_PIDC.NetAttrs['links'].weight > 2.0]
merge_network = gt.graph_merge(filter_genie3, filter_pidc, method = 'intersection')

dn.static_netShow(merge_network, filepath = '/home/mwu/test.html')





