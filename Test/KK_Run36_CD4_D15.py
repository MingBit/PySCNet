
#!/usr/bin/env python3
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
import copy
import numpy as np
from PySCNet.Preprocessing import gnetdata
from PySCNet.BuildNet import gne_dockercaller as gdocker
from PySCNet.NetEnrich import graph_toolkit as gt
from PySCNet.Plotting import show_net as sn


path = '/home/mwu/MING_V9T/PhD_Pro/Test/KK_Run36_D41_CD4/'

#scanpy pipeline
cell_info = pd.read_csv(path + 'Seurat_Data/Cell_Info.tsv', delimiter = ' ', index_col = 0)
Expr = pd.read_csv(path + 'Seurat_Data/top50_marker_raw_data.csv', delimiter=',', index_col = 0)
#Expr = Expr[list(cell_info.loc[cell_info.cluster_id.isin(['1'])].index)]

# create gedata object
Mms_fator = pd.read_csv('/home/mwu/MING_V9T/PhD_Pro/PySCNet/Mus_TF_and_TFcofactor/Mus_musculus_TF.txt', sep = '\t')
Mms_fator['Symbol'] = [x.upper() for x in list(Mms_fator['Symbol'])]
#Expr = pd.DataFrame(adata.X.T, columns=list(adata.obs.index))

Expr.index = [x.upper() for x in list(Expr.index)]
#comm = list(Expr.mean(1).sort_values(ascending = False).head(300).index)

#comm = list(set([x.upper() for x in list(Mms_fator.Symbol)]) & set(Expr.index))
#Expr = Expr.loc[comm, :]

gene_info = pd.DataFrame({'Gene': Expr.index, 'TF_Gene': ['TF' if x in list(Mms_fator['Symbol']) else 'Gene' for x in Expr.index]})

#Expr.to_csv('/home/mwu/MING_V9T/PhD_Pro/Test/Expr.txt', sep = '\t' )
run36_gne = gnetdata.Gnetdata(ExpMatrix = Expr)
run36_gne.CellAttrs = cell_info
run36_gne.GeneAttrs = gene_info
run36_gne_GENIE3 = gdocker.rundocker(run36_gne.deepcopy, method = 'GENIE3')
run36_Links = run36_gne_GENIE3.NetAttrs['links']
run36_Links['cell_clusterid'] = ['All' for i in range(run36_Links.shape[0])]

for i in range(5):
        print(i+1)
        Tmp_Expr = Expr[list(cell_info.loc[cell_info.cluster_id.isin([i+1])].index)]
        Tmp_gne = gnetdata.Gnetdata(ExpMatrix = Tmp_Expr)
        Tmp_gne_GENIE3 = gdocker.rundocker(Tmp_gne.deepcopy, method = 'GENIE3')
        Tmp_Links = Tmp_gne_GENIE3.NetAttrs['links']
        Tmp_Links['cell_clusterid'] = [str(i+1) for j in range(Tmp_Links.shape[0])]
        run36_Links = run36_Links.append(Tmp_Links)



run36_gne_GENIE3.NetAttrs['links'] = copy.deepcopy(run36_Links)
run36_gne_GENIE3.save_as(path + 'Run36_Top50_MarkerGenes.pk')

run36_gne_PIDC = gdocker.rundocker(run36_gne.deepcopy, method = 'PIDC')
run36_gne_GENIE3 = gdocker.rundocker(run36_gne.deepcopy, method = 'GENIE3')
run36_gne_CORR = gdocker.rundocker(run36_gne.deepcopy, method = 'CORR')
time_point = pd.DataFrame({'cell_id': list(Expr.columns), 'time_point': [1 for i in range(Expr.shape[1])]})
run36_gne_SCODE = gdocker.rundocker(run36_gne.deepcopy, method = 'SCODE', time_point = time_point)

run36_gne_PIDC = gdocker.buildnet(run36_gne_PIDC, top = 200)
run36_gne_GENIE3 = gdocker.buildnet(run36_gne_GENIE3, top = 200)

run36_gne_GENIE3 = gt.get_centrality(run36_gne_GENIE3)
run36_gne_GENIE3 = gt.community_detect(run36_gne_GENIE3)

run36_gne_PIDC = gt.get_centrality(run36_gne_PIDC)
run36_gne_PIDC = gt.community_detect(run36_gne_PIDC)

#save as pickle object
run36_gne_GENIE3.save_as(path + 'cluster_4_GENIE3.pk')
run36_gne_PIDC.save_as(path + 'cluster_4_PIDC.pk')

#dn.dynamic_netShow(run36_gne_PIDC, filepath = path + 'cluster_5_wo_CD4UniqueEff_pidc.html')
#dn.dynamic_netShow(run36_gne_GENIE3, filepath = path + 'cluster_5_wo_CD4UniqueEff_genie3.html')


gt_rw_dt = pd.DataFrame(columns=['step_' + str(x) for x in range(10)])
nodes = list(run36_gne_GENIE3.NetAttrs['graph'].node)

for i in range(len(nodes)):
        rw_path = gt.random_walk(gnetdata = run36_gne_GENIE3, start = nodes[i], supervisedby = 'pageRank', steps = 9)
        if len(rw_path) < 10:

                for i in range(10 - len(rw_path)):
                        rw_path.append('NA')

        gt_rw_dt.loc[i] = rw_path

gt_rw_dt.to_csv(path + 'cluster_5_wo_CD4UniqueEff_RandomPath_GENIE3.csv', sep=' ')

#merge test
merge_graph = gt.graph_merge(run36_gne_GENIE3.NetAttrs['links'], run36_gne_PIDC.NetAttrs['links'], method = 'knn')
merge_graph.edges



