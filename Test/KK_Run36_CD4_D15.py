# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 14:09:33 2019
@author: mwu
"""

from __future__ import absolute_import
import sys
import os
import itertools

sys.path.append(os.getenv('HOME') + '/MING_V9T/PhD_Pro/pyscnet')
from pyvis.network import Network
import pandas as pd
import copy
import numpy as np
from pyscnet.Preprocessing import gnetdata
from pyscnet.BuildNet import gne_dockercaller as gdocker
from pyscnet.NetEnrich import graph_toolkit as gt
from pyscnet.Plotting import show_net as sn

import importlib
importlib.reload(gdocker)

path = os.getenv('HOME') + '/MING_V9T/KK_Run36/GEO_Submission/mingwu/GEO_mingwu_NGS_scRNA_Run36/ProcessedData/'
sub_path = os.getenv('HOME') + '//MING_V9T/PhD_Pro/Test/Run36_GRN/'

##############Run36_D41################
genesets = ['Tcf7', 'Gzma', 'Gzmb', 'Gzmk', 'Fasl', 'Id2', 'Tbx21', 'Pdcd1',
            'Cd244', 'Cd160', 'Entpd1', 'Klf2', 'S1pr5', 'S1pr1', 'Cx3cr1',
            'Tox', 'Id3', 'Top2a', 'Birc5', 'Cxcr5', 'Emoes', 'Havcr2',
            'Bcl2l11', 'Cd200', 'Siah2']

selected_genes = ['Tbx21', 'Cx3cr1', 'Klf2', 'S1pr1', 'S1pr4', 'Jak', 'Il18r1', 'Il18rap', 'Klrk1',
                  'Ctsd', 'Lef1', 'Klrc1', 'Klrd1', 'Lgals1', 'Lgals3', 'Vim', 'Pdcd1', 'Cd160', 'Cd244',
                  'Entpd1', 'Cd37', 'Tnfrsf9', 'Tox', 'Nr4a2', 'Psmb8', 'Ctsb', 'Lax1', 'Gzmb', 'Srgn',
                  'Tnfrsf9', 'Tnfrsf4', 'Tnfrsf18', 'Cd200r1', 'Ifngr1', 'Il2rb', 'Sh2d1a', 'Ccl3', 'Ccl4',
                  'Lsp1', 'Vasp', 'Sin3b', 'Ctsb', 'Ctsc', 'Calm1', 'Stat3', 'Id2', 'Ly6a', 'Ly6e', 'Trim30a',
                  'Il10rb', 'Fyb', 'Tnfsf8', 'Ptpn6', 'Cd84', 'Emb', 'Traf1', 'Cd9', 'Icos', 'Cd72', 'Ptpn11',
                  'Cd69', 'Sh2d1a', 'Penk', 'Slamf6', 'Dgka', 'Inpp4b', 'Tespa1', 'Ifnar2', 'Cxcr3', 'Itgb1',
                  'Prf1', 'Fasl', 'Cd244', 'Chn2', 'Fgl2']

selected_genes = list(set(genesets + selected_genes))
# markergene = pd.read_excel(sub_path + 'SelectedGene_MergedState_Avg.xlsx')[1:]
# selected_genes = markergene.geneName
# import raw data and cell clusterID
cell_info = pd.read_csv(path + 'Run36_KK_Cell_ClusterID.txt', delimiter=' ')
Expr = pd.read_csv(path + 'Run36_KK_UMI_Top2000_Raw.txt', delimiter=' ', index_col=0)
Expr.columns = [x.replace('.', '-') for x in Expr.columns]

selected_Expr = Expr[Expr.index.isin(selected_genes)]
# Expr = Expr[list(cell_info.loc[cell_info.cluster_id.isin(['1'])].index)]

# create gedata object
Mms_fator = pd.read_csv(os.getenv('HOME') + '/MING_V9T/PhD_Pro/pyscnet/Mus_TF_and_TFcofactor/Mus_musculus_TF.txt',
                        sep='\t')
Mms_fator['Symbol'] = [x.upper() for x in list(Mms_fator['Symbol'])]
gene_info = pd.DataFrame({'Gene': selected_Expr.index,
                          'TF_Gene': ['TF' if x.upper() in list(Mms_fator['Symbol']) else 'Gene' for x in
                                      selected_Expr.index]})
# prepare reference from stringppi
gne_Ref_raw = pd.read_csv(sub_path + 'string_interaction_selectedGenes.tsv', sep='\t')[['node1', 'node2']]
gne_Ref_raw.columns = ['source', 'target']
gne_Ref_raw['weight'] = 1

Full_Ref = pd.DataFrame(itertools.permutations(selected_Expr.index, 2), columns=['source', 'target'])
gne_Ref = pd.merge(Full_Ref, gne_Ref_raw, how='outer').fillna(0)

run36_gne = gnetdata.Gnetdata(ExpMatrix=selected_Expr)
run36_gne.CellAttrs = cell_info
run36_gne.GeneAttrs = gene_info

param_grid = {
    'n_estimators': [50, 100],
    'max_depth': [10, 20, 50],
    'max_features': [3, 5],
    'min_samples_leaf': [3, 4, 5],
    'criterion': ['gini', 'entropy']
}

parameters = {'reference_links': gne_Ref,
              'p': 0.1, 'q': 1.0,
              'size': 5, 'walk_len': 5,
              'num_walks': 1000,
              'workers': 10, 'n_pc': 5,
              'param_grid': param_grid}

run36_gne_PIDC = gdocker.rundocker(run36_gne.deepcopy, method='PIDC')
run36_gne_GENIE3 = gdocker.rundocker(run36_gne.deepcopy, method='GENIE3')
run36_gne_GRNBOOST2 = gdocker.rundocker(run36_gne.deepcopy, method='GRNBOOST2')
run36_gne_NODE2VEC = gdocker.rundocker(run36_gne.deepcopy, method='SCNODE2VEC', parameters=parameters)

Genie3_Links = run36_gne_GENIE3.NetAttrs['links'].reindex(
    run36_gne_GENIE3.NetAttrs['links'].weight.abs().sort_values(ascending=False).index).head(300)

# Node2Vec_Links = run36_gne_node2vec.NetAttrs['links'].reindex(
#     run36_gne_node2vec.NetAttrs['links'].weight.abs().sort_values(ascending=False).index).head(300)

# Genie3_Links.to_excel(sub_path + 'Run36_Genie3_Top300.xlsx')
# Node2Vec_Links.to_excel(sub_path + 'Run36_Node2Vec_Top300_NoRef.xlsx')
#
# Genie3_Links['cell_clusterid'] = ['All' for i in range(Genie3_Links.shape[0])]

# for i in range(5):
#     print(i + 1)
#     Tmp_Expr = selected_Expr.iloc[:, cell_info.loc[cell_info.cluster_id.isin([i + 1])].index]
#     Tmp_gne = gnetdata.Gnetdata(ExpMatrix=Tmp_Expr)
#     Tmp_gne_GENIE3 = gdocker.rundocker(Tmp_gne.deepcopy, method='GENIE3')
#     Tmp_Links = Tmp_gne_GENIE3.NetAttrs['links']
#     Tmp_Links['cell_clusterid'] = [str(i + 1) for j in range(Tmp_Links.shape[0])]
#     Genie3_Links = Genie3_Links.append(Tmp_Links)
#
# run36_gne_GENIE3.NetAttrs['links'] = Genie3_Links
#
# run36_gne_GENIE3 = gdocker.buildnet(run36_gne_GENIE3, top=300)
# run36_gne_GENIE3 = gt.get_centrality(run36_gne_GENIE3)
# run36_gne_GENIE3 = gt.community_detect(run36_gne_GENIE3)
run36_gne_GENIE3.save_as(os.getenv('HOME') + '/run36_gne_GENIE3_test.pk')
run36_gne_GRNBOOST2.save_as(os.getenv('HOME') + '/run36_gne_GRNBOOST2_test.pk')
# run36_gne_node2vec = gdocker.buildnet(run36_gne_node2vec, top=500)
# run36_gne_node2vec = gt.get_centrality(run36_gne_node2vec)
# run36_gne_node2vec = gt.community_detect(run36_gne_node2vec)
#
# run36_gne_node2vec.save_as(sub_path + 'run36_gne_node2vec_NoRef.pk')

# importlib.reload(sn)
# sn.dynamic_netShow(run36_gne_GENIE3, filename=path + 'run36_genie3_top300.html')
# sn.dynamic_netShow(run36_gne_node2vec, filename=sub_path + 'run36_node2vec_top300_NoRef.html')
# sn.static_netShow(run36_gne_GENIE3, filename=path + 'run36_genie3_top300.pdf')
# sn.static_netShow(run36_gne_node2vec, filename=path + 'run36_node2vec_top300.pdf')


# gt_rw_dt = pd.DataFrame(columns=['step_' + str(x) for x in range(10)])
# nodes = list(run36_gne_GENIE3.NetAttrs['graph'].node)
#
# for i in range(len(nodes)):
#     rw_path = gt.random_walk(gnetdata=run36_gne_GENIE3, start=nodes[i], supervisedby='pageRank', steps=9)
#     if len(rw_path) < 10:
#
#         for i in range(10 - len(rw_path)):
#             rw_path.append('NA')
#
#     gt_rw_dt.loc[i] = rw_path
#
# gt_rw_dt.to_csv(path + 'cluster_5_wo_CD4UniqueEff_RandomPath_GENIE3.csv', sep=' ')
