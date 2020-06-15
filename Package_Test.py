#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 14:35:57 2019

@author: angelawu
"""

from __future__ import absolute_import
import pandas as pd
import numpy as np
import sys
import os
import seaborn as sns

sys.path.append('/home/mwu/MING_V9T/PhD_Pro/PySCNet/')
import _pickle as pk
from PySCNet.Preprocessing import gnetdata
from PySCNet.Preprocessing import general_pipeline as pipeline
from PySCNet.BuildNet import gne_dockercaller as gdocker
from PySCNet.NetEnrich import graph_toolkit as gt
from PySCNet.Plotting import show_net as sn
import matplotlib.pyplot as plt
import networkx as nx
from PySCNet.BuildNet import gne_modelcaller as gmodel

path = '/home/mwu/MING_V9T/PhD_Pro/TODOLIST/BN_Test/SC_Published_Data/'
# HSC_Data = pd.read_excel(path + 'HSC_DATA/Moignard_HSC_Data.xlsx', index_col = 0, sheetname = 0)
# HSC_Data[HSC_Data == 25] = 0
# HSC_Ref = pd.read_csv(path + 'HSC_DATA/Moignard_HSC_ReferenceEdges.tsv', sep = '\t')

# path = '/home/mwu/MING_V9T/PhD_Pro/Test/Simulation/GNW_Data'

Sim_1 = pd.read_csv(path + 'simulated_datasets/100_yeast3_large_dropouts_low.txt', sep='\t', index_col=0)
Sim_1_Ref = pd.read_csv(path + 'goldstandards/100_yeast3.tsv', sep='\t', header=None)
Sim_1_Ref.columns = ['node1', 'node2', 'value']
Sim_1_Ref = Sim_1_Ref.loc[Sim_1_Ref['value'] > 0]

Sim_2 = pd.read_csv(path + 'simulated_datasets/100_yeast3_medium.txt', sep='\t', index_col=0)
Sim_2_Ref = pd.read_csv(path + 'goldstandards/100_yeast3.tsv', sep='\t', header=None)
Sim_2_Ref.columns = ['node1', 'node2', 'value']
Sim_2_Ref = Sim_2_Ref.loc[Sim_2_Ref['value'] > 0]

HSC_Data = pd.read_excel(path + 'HSC_DATA/Moignard_HSC_Data.xlsx', index_col=0, sheet_name=0).T
HSC_Data[HSC_Data == 25] = 0
HSC_Data = HSC_Data[HSC_Data.sum(1) > 0]
HSC_Ref = pd.read_csv(path + 'HSC_DATA/Moignard_HSC_ReferenceEdges.tsv', sep='\t')

ESC_Data = pd.read_excel(path + 'ESC_DATA/GSE59892_ESC_Dataset.xlsx', index_col=0, sheet_name=0)
ESC_Gene_Symbol = pd.read_csv(path + 'ESC_DATA/GSE59892_ESC_Genesymobol.txt', sep='\t')
ESC_Gene_Symbol_dict = dict(zip(ESC_Gene_Symbol.Ref_ID, ESC_Gene_Symbol.GeneSymbol))
ESC_Ref = pd.read_csv(path + 'ESC_DATA/GSE59892_ESC_ReferenceEdges.tsv', sep='\t')

ESC_Data.index = list(ESC_Gene_Symbol_dict.get(i) for i in list(ESC_Data.index))
ESC_Data = ESC_Data[ESC_Data.sum(1) > 0]

# random set 0 for GNW simulation
path_2 = '/home/mwu/MING_V9T/PhD_Pro/Test/Simulation/GNW_Data/'
Ecoli_G50_C200 = pd.read_csv(path_2 + 'G50_C200/Ecoli-1_dream4_timeseries.tsv', sep='\t', index_col=0)

Ecoli_G50_C200_Ref = pd.read_csv(path_2 + 'G50_C200/')

# Simlated data from BOOLODE
BODE_Expr = pd.read_csv('/home/mwu/MING_V9T/PhD_Pro/Test/Simulation/BoolODE_Data/G300_C1000/ExpressionData.csv',
                        sep=',', index_col=0)
BODE_Ref = pd.read_csv('/home/mwu/MING_V9T/PhD_Pro/Test/Simulation/BoolODE_Data/G300_C1000/refNetwork.csv', sep=',',
                       index_col=None)
BODE_Ref['value'] = 1
BODE_Ref = BODE_Ref.drop(columns='Type')
BODE_Ref.columns = ['node1', 'node2', 'value']
# =============================================================================
# Test with Preprocessing and docker run
# =============================================================================


Sim1_gne = gnetdata.Gnetdata(Sim_1)
Sim2_gne = gnetdata.Gnetdata(Sim_2)
ESC_gne = gnetdata.Gnetdata(ESC_Data)
HSC_gne = gnetdata.Gnetdata(HSC_Data)
BODE_gne = gnetdata.Gnetdata(BODE_Expr)

# gne_tmp_GENIE3 = gdocker.rundocker(tmp_gne.deepcopy, method='GENIE3')
gne_Sim1_PIDC = gdocker.rundocker(Sim1_gne.deepcopy, method='PIDC')
gne_Sim1_GENIE3 = gdocker.rundocker(Sim1_gne.deepcopy, method='GENIE3')
gne_Sim1_CORR = gdocker.rundocker(Sim1_gne.deepcopy, method='CORR')
time_point = pd.DataFrame({'cell_id': list(Sim_1.columns), 'time_point': [1 for i in range(Sim1_gne.shape[1])]})
gne_Sim1_SCODE = gdocker.rundocker(Sim1_gne.deepcopy, method='SCODE', time_point=time_point)

gne_Sim2_PIDC = gdocker.rundocker(Sim2_gne.deepcopy, method='PIDC')
gne_Sim2_GENIE3 = gdocker.rundocker(Sim2_gne.deepcopy, method='GENIE3')
gne_Sim2_CORR = gdocker.rundocker(Sim2_gne.deepcopy, method='CORR')
time_point = pd.DataFrame({'cell_id': list(Sim_2.columns), 'time_point': [1 for i in range(Sim2_gne.shape[1])]})
gne_Sim2_SCODE = gdocker.rundocker(Sim2_gne.deepcopy, method='SCODE', time_point=time_point)

gne_ESC_PIDC = gdocker.rundocker(ESC_gne.deepcopy, method='PIDC')
gne_ESC_GENIE3 = gdocker.rundocker(ESC_gne.deepcopy, method='GENIE3')
gne_ESC_CORR = gdocker.rundocker(ESC_gne.deepcopy, method='CORR')
time_point = pd.DataFrame({'cell_id': list(ESC_Data.columns), 'time_point': [1 for i in range(ESC_gne.shape[1])]})
gne_ESC_SCODE = gdocker.rundocker(ESC_gne.deepcopy, method='SCODE', time_point=time_point)

gne_HSC_PIDC = gdocker.rundocker(HSC_gne.deepcopy, method='PIDC')
gne_HSC_GENIE3 = gdocker.rundocker(HSC_gne.deepcopy, method='GENIE3')
gne_HSC_CORR = gdocker.rundocker(HSC_gne.deepcopy, method='CORR')
time_point = pd.DataFrame({'cell_id': list(HSC_Data.columns), 'time_point': [1 for i in range(HSC_gne.shape[1])]})
gne_HSC_SCODE = gdocker.rundocker(HSC_gne.deepcopy, method='SCODE', time_point=time_point)

gne_BODE_PIDC = gdocker.rundocker(BODE_gne.deepcopy, method='PIDC')
gne_BODE_GENIE3 = gdocker.rundocker(BODE_gne.deepcopy, method='GENIE3')
gne_BODE_CORR = gdocker.rundocker(BODE_gne.deepcopy, method='CORR')
time_point = pd.DataFrame({'cell_id': list(BODE_Expr.columns), 'time_point': [1 for i in range(BODE_Expr.shape[1])]})
gne_BODE_SCODE = gdocker.rundocker(BODE_gne.deepcopy, method='SCODE', time_point=time_point)

gne_BODE_PIDC.save_as(path + 'gne_BODE_PIDC.pk')
gne_BODE_GENIE3.save_as(path + 'gne_BODE_GENIE3.pk')
gne_BODE_CORR.save_as(path + 'gne_BODE_CORR.pk')
# gne_HSC_SCODE.save_as(path + 'gnetData/HSC_SCODE.pk')


# Add node centralities
sim1_list = [gne_Sim1_CORR, gne_Sim1_GENIE3, gne_Sim1_PIDC, gne_Sim1_SCODE]
sim2_list = [gne_Sim2_CORR, gne_Sim2_GENIE3, gne_Sim2_PIDC, gne_Sim2_SCODE]

for obj in sim1_list:
    obj = gdocker.buildnet(obj, top=50)
    obj = gt.community_detect(obj)
    obj = gt.get_centrality(obj)


def build_table(gnedata):
    links = gnedata.NetAttrs['links']
    centralities = gnedata.NetAttrs['centralities']

    centralities['avg'] = centralities.mean(axis=1)
    centrality_avg = list()
    for source, target in zip(links.source, links.target):

        if source in list(centralities.node):
            source_val = centralities[centralities.node == source].avg.values[0]
        else:
            source_val = 0
        if target in list(centralities.node):
            target_val = centralities[centralities.node == target].avg.values[0]
        else:
            target_val = 0

        centrality_avg.append(source_val + target_val)

    links['centrality_avg'] = centrality_avg
    return (links)


gne_Sim1_PIDC = gnetdata.load_Gnetdata_object(path + 'gnetData/Sim1_PIDC.pk')
gne_Sim1_GENIE3 = gnetdata.load_Gnetdata_object(path + 'gnetData/Sim1_GENIE3.pk')
gne_Sim1_CORR = gnetdata.load_Gnetdata_object(path + 'gnetData/Sim1_CORR.pk')
gne_Sim1_SCODE = gnetdata.load_Gnetdata_object(path + 'gnetData/Sim1_SCODE.pk')

gne_Sim2_PIDC = gnetdata.load_Gnetdata_object(path + 'gnetData/Sim2_PIDC.pk')
gne_Sim2_GENIE3 = gnetdata.load_Gnetdata_object(path + 'gnetData/Sim2_GENIE3.pk')
gne_Sim2_CORR = gnetdata.load_Gnetdata_object(path + 'gnetData/Sim2_CORR.pk')
gne_Sim2_SCODE = gnetdata.load_Gnetdata_object(path + 'gnetData/Sim2_SCODE.pk')

gne_HSC_CORR = gnetdata.load_Gnetdata_object(path + 'gnetData/ESC_CORR.pk')
gne_HSC_PIDC = gnetdata.load_Gnetdata_object(path + 'gnetData/ESC_PIDC.pk')
gne_HSC_GENIE3 = gnetdata.load_Gnetdata_object(path + 'gnetData/ESC_GENIE3.pk')
gne_HSC_SCODE = gnetdata.load_Gnetdata_object(path + 'gnetData/ESC_SCODE.pk')

# Sim1_links_dict = {'genie3': build_table(gne_Sim1_GENIE3), 'corr': build_table(gne_Sim1_CORR),
#              'pidc': build_table(gne_Sim1_PIDC), 'scode': build_table(gne_Sim1_SCODE)}
#
# Sim2_links_dict = {'genie3': build_table(gne_Sim2_GENIE3), 'corr': build_table(gne_Sim2_CORR),
#              'pidc': build_table(gne_Sim2_PIDC), 'scode': build_table(gne_Sim2_SCODE)}
#
# HSC_links_dict = {'genie3': build_table(gne_HSC_GENIE3), 'corr': build_table(gne_HSC_CORR),
#              'pidc': build_table(gne_HSC_PIDC), 'scode': build_table(gne_HSC_SCODE)}
#
# ESC_links_dict = {'genie3': build_table(gne_ESC_GENIE3), 'corr': build_table(gne_ESC_CORR),
#              'pidc': build_table(gne_ESC_PIDC), 'scode': build_table(gne_ESC_SCODE)}


Sim1_links_dict = {'genie3': gne_Sim1_GENIE3.NetAttrs['links'], 'corr': gne_Sim1_CORR.NetAttrs['links'],
                   'pidc': gne_Sim1_PIDC.NetAttrs['links'], 'scode': gne_Sim1_SCODE.NetAttrs['links']}

Sim2_links_dict = {'genie3': gne_Sim2_GENIE3.NetAttrs['links'], 'corr': gne_Sim2_CORR.NetAttrs['links'],
                   'pidc': gne_Sim2_PIDC.NetAttrs['links'], 'scode': gne_Sim2_SCODE.NetAttrs['links']}
#
#
ESC_links_dict = {'genie3': gne_ESC_GENIE3.NetAttrs['links'], 'corr': gne_ESC_CORR.NetAttrs['links'],
                  'pidc': gne_ESC_PIDC.NetAttrs['links'], 'scode': gne_ESC_SCODE.NetAttrs['links']}
#
HSC_links_dict = {'genie3': gne_HSC_GENIE3.NetAttrs['links'], 'corr': gne_HSC_CORR.NetAttrs['links'],
                  'pidc': gne_HSC_PIDC.NetAttrs['links'], 'scode': gne_HSC_SCODE.NetAttrs['links']}

BODE_links_dict = {'genie3': gne_BODE_GENIE3.NetAttrs['links'], 'corr': gne_BODE_CORR.NetAttrs['links'],
                   'pidc': gne_BODE_PIDC.NetAttrs['links']}

Eclass = gt.ensemble_classifier(Sim2_links_dict, threshold=0.5, model='autoSK', seed=5)
Eclass_pos = Eclass.sort_values('weight', ascending=False).head(500)

Eclass_RF = gt.ensemble_classifier(Sim2_links_dict, threshold=0.5, num_trees=1000, seed=5, model='RF')
Eclass_pos_RF = Eclass_RF.sort_values('weight', ascending=False).head(500)

RF_fpr, RF_tpr, RF_auc = Eva_Test.test_run(links=Eclass_pos_RF,
                                           Ref_links=Sim_2_Ref, input_dataset=Sim_2)

Auto_RF_fpr, Auto_RF_tpr, Auto_RF_auc = Eva_Test.test_run(links=Eclass_pos,
                                                          Ref_links=Sim_2_Ref, input_dataset=Sim_2)

PIDC_fpr, PIDC_tpr, PIDC_auc = Eva_Test.test_run(
    links=gne_Sim2_PIDC.NetAttrs['links'].sort_values('weight', ascending=False).head(500),
    Ref_links=Sim_2_Ref, input_dataset=Sim_2)

GENIE3_fpr, GENIE3_tpr, GENIE3_auc = Eva_Test.test_run(
    links=gne_Sim2_GENIE3.NetAttrs['links'].sort_values('weight', ascending=False).head(500),
    Ref_links=Sim_2_Ref, input_dataset=Sim_2)

CORR_fpr, CORR_tpr, CORR_auc = Eva_Test.test_run(
    links=gne_Sim2_CORR.NetAttrs['links'].sort_values('weight', ascending=False).head(500),
    Ref_links=Sim_2_Ref, input_dataset=Sim_2)

SCODE_fpr, SCODE_tpr, SCODE_auc = Eva_Test.test_run(
    links=gne_HSC_SCODE.NetAttrs['links'].sort_values('weight', ascending=False).head(500),
    Ref_links=HSC_Ref, input_dataset=HSC_Data)

#######################test for node2vec################################
tmp = model_node2vec._build_coexp_graph(Sim_2)
tmp['weight'] = tmp['weight'].abs()

DL_fpr = dict()
DL_tpr = dict()
DL_auc = dict()

for dimensions in [32]:
    for walk_length in [10, 20]:
        for num_walks in [3000]:
            tmp1 = my_node2vec.my_node2vec(Sim_2, p=1, q=0.25, walk_len=walk_length,
                                           num_walks=num_walks, size=dimensions)
            #            tmp1 = model_node2vec._build_node2vec_model(tmp, p = 1, q = 4, walk_length=walk_length,
            #                                        num_walks=num_walks, dimensions=dimensions, method='cosine', workers = 12)
            #            tmp1['weight'] = tmp1['weight'].abs()
            tmp1 = gdocker._remove_duplicate(tmp1)
            # tmp1[tmp1['weight'] < np.quantile(tmp1['weight'], 0.75)]['weight']

            # linkage_table = pd.DataFrame(list(nx.from_pandas_adjacency(tmp1).edges), columns = ['source', 'target'])
            # linkage_table['weight'] = [tmp1[linkage_table['source'][i]][linkage_table['target'][i]] for i in range(linkage_table.shape[0])]

            tmp_fpr, tmp_tpr, tmp_auc = Eva_Test.test_run(
                links=tmp1[tmp1['weight'] > np.quantile(tmp1['weight'], 0.90)],
                Ref_links=Sim_2_Ref, input_dataset=Sim_2)
            DL_fpr['dims_' + str(dimensions) + 'walk_len_' + str(walk_length) + 'num_walks_' + str(num_walks)] = tmp_fpr
            DL_tpr['dims_' + str(dimensions) + 'walk_len_' + str(walk_length) + 'num_walks_' + str(num_walks)] = tmp_tpr
            DL_auc['dims_' + str(dimensions) + 'walk_len_' + str(walk_length) + 'num_walks_' + str(num_walks)] = tmp_auc

#########################test for node2vec###############################

fpr_dict = {'GENIE3': GENIE3_fpr, 'PIDC': PIDC_fpr, 'CORR': CORR_fpr,
            'Node2Vec': DL_fpr['dims_32walk_len_20num_walks_3000']}
tpr_dict = {'GENIE3': GENIE3_tpr, 'PIDC': PIDC_tpr, 'CORR': CORR_tpr,
            'Node2Vec': DL_tpr['dims_32walk_len_20num_walks_3000']}
auc_dict = {'GENIE3': GENIE3_auc, 'PIDC': PIDC_auc, 'CORR': CORR_auc,
            'Node2Vec': DL_auc['dims_32walk_len_20num_walks_3000']}

import seaborn as sns

colors = sns.color_palette().as_hex()[1:8]
plt.figure(figsize=(8, 8))
for i, color in zip(list(auc_dict.keys()), colors):
    plt.plot(fpr_dict[i], tpr_dict[i], color=color,
             label='ROC curve of method {0} (area = {1:0.2f})'.format(i, auc_dict[i]))

# for i, color in zip(list(DL_auc.keys()), colors):
#        plt.plot(DL_fpr[i], DL_tpr[i], color = color, label =  'ROC of method {0} (area = {1:0.2f})'.format(i, DL_auc[i]))

# plt.plot(RF_fpr, RF_tpr)
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.xlabel('False Positive Rate', fontsize=20)
plt.ylabel('True Positive Rate', fontsize=20)
plt.title('Sim2: Algorithms performance', fontsize=25)
plt.legend(loc="lower right")

plt.savefig('/home/mwu/Sim2.pdf')
plt.close()

# obj = gne_HSC_PIDC.deepcopy
obj = gnetdata.Gnetdata(HSC_Data)
obj.NetAttrs['links'] = Eclass_RF
obj = gdocker.buildnet(obj, top=200)
obj = gt.community_detect(obj)
obj = gt.get_centrality(obj)
random_path = gt.random_walk(obj, start='Gata1', supervisedby='pageRank', steps=10)

tmp1 = ['Gata1']
tmp1.extend(list(obj.NetAttrs['graph']['Gata1'].keys()))
sn.static_netShow(obj, '/home/mwu/Gata_Path.pdf', figure_size=[5, 5],
                  #                  start='Gata1', neighhours_highlight=True,
                  random_path=random_path, path_highlight=True,
                  scale=2, node_size=100, font_size=4, edge_color='grey', width=1)

sn.static_netShow(obj, '/home/mwu/test_7.pdf', figure_size=[5, 2], scale=2, node_size=100, width=0.8,
                  edge_color='grey', font_size=5, random_path=path, path_highlight=True)

sn.dynamic_netShow(obj, filename='/home/mwu/test.html')
# simulate linkage table
tmp = pd.DataFrame()

# make network with node encoded by continues colors
import seaborn as sns

path = os.getenv('HOME') + '/MING_V9T/PhD_Pro/Test/Run36_GRN/MergedState_MarkersGRN/'
run36_gne_GENIE3 = gnetdata.load_Gnetdata_object(path + 'GENIE3_GRN/run36_gne_GENIE3_allclusters.pk')
gene_avg = pd.read_excel(os.getenv('HOME') + '/MING_V9T/PhD_Pro/Test/Run36_GRN/SelectedGene_MergedState_Avg.xlsx',
                         index_col=0)

run36_gne_GENIE3 = gdocker.buildnet(run36_gne_GENIE3, top=300)
run36_gne_GENIE3 = gt.get_centrality(run36_gne_GENIE3)
run36_gne_GENIE3 = gt.community_detect(run36_gne_GENIE3)

G = run36_gne_GENIE3.NetAttrs['graph']
gene_avg = gene_avg.reindex(G.nodes())
pos = nx.random_layout(G)
colors = sns.color_palette().as_hex() + sns.color_palette('Paired', 100).as_hex()
highlight_node = ['Tcf7', 'Pdcd1']

for i in gene_avg.columns:
    plt.figure(figsize=(15, 12))
    nc = nx.draw_networkx_nodes(G, pos=pos, node_size=1500, node_color=gene_avg[i], scale=3)
    ec = nx.draw_networkx_edges(G, pos=pos, scale=3, edge_color='grey')
    labels = nx.draw_networkx_labels(G, pos=pos, font_color='white', font_size=10)
    # highligh nodes edges
    n = 0
    for node in highlight_node:
        n = n + 1
        neighbors = [n for n in G.neighbors(node)]
        for j in range(len(neighbors)):
            nx.draw_networkx_edges(G, pos, {(node, neighbors[j]): str(j)},
                                   edge_color=colors[n], width=4)
    plt.colorbar(nc)
    plt.axis('off')
    plt.show()
    plt.title(str(i), fontsize=20)
    plt.savefig(path + 'GENIE3_GRN/Genie3_' + str(i) + '_randm.pdf')
