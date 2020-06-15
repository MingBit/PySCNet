#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 14:16:47 2019

@author: mwu
"""
from __future__ import absolute_import
import os
import sys

sys.path.append(os.getenv("HOME") + '/MING_V9T/PhD_Pro/PySCNet/')

from importlib import reload
from PySCNet.Preprocessing import gnetdata
from PySCNet.BuildNet import gne_dockercaller as gdocker
from AlgEval import Eva_Test
from AlgEval import basic_functions as bf
import matplotlib.pyplot as plt
from PySCNet.BuildNet import gne_modelcaller as gmodel
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2, genie3

# path = os.getenv("HOME") + '/MING_V9T/PhD_Pro/Test/Simulation/BEELINE-data/inputs/Curated/HSC/'
# filename = 'G30_C1000_Clu_4_Tra_3'

p = 1.0
q = 0.1
dropRate = 0.75

####################HSC data##########################
path = os.getenv('HOME') + '/MING_V9T/PhD_Pro/TODOLIST/BN_Test/SC_Published_Data/'
import pandas as pd
import itertools
import copy
import numpy as np

HSC_Data = pd.read_excel(path + 'HSC_DATA/Moignard_HSC_Data.xlsx', index_col=0, sheet_name=0).T
HSC_Ref_raw = pd.read_csv(path + 'HSC_DATA/Moignard_HSC_ReferenceEdges.tsv', sep='\t')[['node1', 'node2']]
HSC_Ref_raw.columns = ['source', 'target']
HSC_Ref_raw['weight'] = 1

Full_Ref = pd.DataFrame(itertools.permutations(HSC_Data.index, 2), columns=['source', 'target'])
HSC_Ref = pd.merge(Full_Ref, HSC_Ref_raw, how='outer').fillna(0)

HSC_Ref_dropout = copy.deepcopy(HSC_Ref)
pos = list(HSC_Ref_dropout[abs(HSC_Ref_dropout.weight) == 1].index)
replaceNum = np.random.choice(pos, int(len(pos) * dropRate), replace=False)
HSC_Ref_dropout.loc[replaceNum, 'weight'] = 0

Sim = HSC_Data
Sim_Ref = HSC_Ref
Sim_Ref_dropout = HSC_Ref_dropout
top_edges = HSC_Ref_raw.shape[0]

#######################ESC data######################################
ESC_Data = pd.read_excel(path + 'ESC_DATA/GSE59892_ESC_Dataset.xlsx', index_col=0, sheet_name=0)
ESC_Gene_Symbol = pd.read_csv(path + 'ESC_DATA/GSE59892_ESC_Genesymobol.txt', sep='\t')

ESC_Gene_Symbol_dict = dict(zip(ESC_Gene_Symbol.Ref_ID, ESC_Gene_Symbol.GeneSymbol))
ESC_Ref_raw = pd.read_csv(path + 'ESC_DATA/GSE59892_ESC_ReferenceEdges.tsv', sep='\t')[['node1', 'node2']]
ESC_Data.index = list(ESC_Gene_Symbol_dict.get(i) for i in list(ESC_Data.index))
ESC_Data = ESC_Data[ESC_Data.sum(1) > 0]

ESC_Ref_raw.columns = ['source', 'target']
ESC_Ref_raw['weight'] = 1

Full_Ref = pd.DataFrame(itertools.permutations(ESC_Data.index, 2), columns=['source', 'target'])
ESC_Ref = pd.merge(Full_Ref, ESC_Ref_raw, how='outer').fillna(0)

ESC_Ref_dropout = copy.deepcopy(ESC_Ref)
pos = list(ESC_Ref_dropout[abs(ESC_Ref_dropout.weight) == 1].index)
replaceNum = np.random.choice(pos, int(len(pos) * dropRate), replace=False)
ESC_Ref_dropout.loc[replaceNum, 'weight'] = 0

Sim = ESC_Data
Sim_Ref = ESC_Ref
Sim_Ref_dropout = ESC_Ref_dropout
top_edges = ESC_Ref_raw.shape[0]
#######################ESC data######################################

filename = sys.argv[1]
p = sys.argv[2]
q = sys.argv[3]
dropRate = sys.argv[4]

# Sim, Sim_Ref, Sim_Ref_dropout, top_edges = bf.read_data(path + filename, float(dropRate))
Sim_gne = gnetdata.Gnetdata(Sim)

gne_Sim_PIDC = gdocker.rundocker(Sim_gne.deepcopy, method='PIDC')
gne_Sim_GENIE3 = gdocker.rundocker(Sim_gne.deepcopy, method='GENIE3')
gne_Sim_CORR = gdocker.rundocker(Sim_gne.deepcopy, method='CORR')

# reload(gmodel.model_node2vec)
gne_Sim_node2vec = gmodel.call_node2vec(Sim_gne.deepcopy, Sim_Ref_dropout, method='pearson',
                                        p=1.0, q=0.1, dim_list=[5], walk_list=[5],
                                        num_walks_list=[1000], n_pc=5, workers=4)

PIDC_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'],
                     Eva_Test.test_run(links=gne_Sim_PIDC.NetAttrs['links'].reindex(
                         gne_Sim_PIDC.NetAttrs['links'].weight.abs().sort_values(ascending=False).index).head(
                         top_edges),
                         Ref_links=Sim_Ref, input_dataset=Sim)))

GENIE3_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'],
                       Eva_Test.test_run(links=gne_Sim_GENIE3.NetAttrs['links'].reindex(
                           gne_Sim_GENIE3.NetAttrs['links'].weight.abs().sort_values(ascending=False).index).head(
                           top_edges),
                           Ref_links=Sim_Ref, input_dataset=Sim)))

CORR_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'],
                     Eva_Test.test_run(links=gne_Sim_CORR.NetAttrs['links'].reindex(
                         gne_Sim_CORR.NetAttrs['links'].weight.abs().sort_values(ascending=False).index).head(
                         top_edges),
                         Ref_links=Sim_Ref, input_dataset=Sim)))

Node2Vec_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'],
                         Eva_Test.test_run(links=gne_Sim_node2vec.NetAttrs['links'].reindex(
                             gne_Sim_node2vec.NetAttrs['links'].weight.abs().sort_values(ascending=False).index).head(
                             top_edges), Ref_links=Sim_Ref, input_dataset=Sim)))

input_path = os.getenv("HOME") + '/MING_V9T/PhD_Pro/Test/Simulation/Node2Vec_BEELINE_Data_Res/'
output_path = os.getenv("HOME") + '/MING_V9T/PhD_Pro/Test/Simulation/Node2Vec_BEELINE_Data_Res/PCA_ROC_Res/'

node_dict_list = list()
for i in range(10):
    gne_Sim_node2vec = gnetdata.load_Gnetdata_object(input_path +
                                                     filename + '_PCA/' + filename + '_node2vec_p_' +
                                                     str(p) + '_q_' + str(q) + '_dim_5_walk_5_dropRate_' +
                                                     str(dropRate) + '_repeat_' + str(i + 1) + '.pk')
    top_links = gne_Sim_node2vec.NetAttrs['links'].reindex(
        gne_Sim_node2vec.NetAttrs['links'].weight.abs().sort_values(ascending=False).index).head(top_edges)

    node_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'],
                         Eva_Test.test_run(links=top_links, Ref_links=Sim_Ref, input_dataset=Sim)))

    node_dict_list.append(node_dict)

fig = plt.figure(figsize=(20, 15))
grid = plt.GridSpec(2, 2, wspace=0.2, hspace=0.2)
ax1 = fig.add_subplot(grid[0, 0])
ax2 = fig.add_subplot(grid[0, 1])
ax3 = fig.add_subplot(grid[1, 0:])
fig.suptitle(filename + '_p_' + str(p) + '_q_' + str(q) + ':Algorithms performance', fontsize=30)
bf.build_curves(ax1, node_dict_list, GENIE3_dict, PIDC_dict, CORR_dict, 'ROC', filename, p, q, dropRate)
bf.build_curves(ax2, node_dict_list, GENIE3_dict, PIDC_dict, CORR_dict, 'PCR', filename, p, q, dropRate)
bf.build_plot(ax3, node_dict_list, GENIE3_dict, PIDC_dict, CORR_dict)

# fig.savefig(output_path + filename + '_PCA_p_' + str(p) + '_q_' + str(q) + '_dropRate_' +
#             str(dropRate) + '_topTP.pdf', bbox_inches='tight')

fig.savefig('/Users/angelawu/Desktop/' + filename + '_p_' + str(p) + '_q_' + str(q) + '_dropRate_' +
            str(dropRate) + '_topTP.pdf', bbox_inches='tight')

gne_Sim_node2vec.save_as('/Users/angelawu/Desktop/HomeOffice/' + filename + '_node2vec_p_' + str(p) + '_q_' + str(q) +
                         '_dim_5_walk_5' + '_dropRate_' + str(dropRate) + '.pk')

gne_dict = {'GENIE3': gne_Sim_GENIE3, 'PIDC': gne_Sim_PIDC, 'CORR': gne_Sim_CORR, 'Node2Vec': gne_Sim_node2vec}
link_dict = dict()
for key in gne_dict.keys():
    gne = gne_dict[key]
    links = gne.NetAttrs['links'].reindex(
        gne.NetAttrs['links'].weight.abs().sort_values(ascending=False).index).head(top_edges)
    edge1 = links[['source', 'target']].apply(lambda x: '_'.join(x), axis=1)
    edge2 = links[['target', 'source']].apply(lambda x: '_'.join(x), axis=1)

    link_dict[key] = edge1.append(edge2)

merge_links = np.unique(link_dict['GENIE3'].append(link_dict['PIDC']).append(link_dict['CORR']))
unique_links = list(set(link_dict['Node2Vec']) - set(merge_links))

Sim_Ref_Raw = Sim_Ref[Sim_Ref.weight != 0][['source', 'target']]
Sim_Ref_Raw_edge = list(Sim_Ref_Raw[['source', 'target']].apply(lambda x: '_'.join(x), axis=1).append(
    Sim_Ref_Raw[['target', 'source']].apply(lambda x: '_'.join(x), axis=1)))

unique_falseNeg_links = set(unique_links) - set(Sim_Ref_Raw_edge)
unique_falseNeg_links_df = pd.DataFrame({'source': [x.split('_')[0] for x in unique_falseNeg_links],
                                         'target': [x.split('_')[1] for x in unique_falseNeg_links]})

unique_falseNeg_links_df = gdocker._remove_duplicate(unique_falseNeg_links_df)
unique_falseNeg_links_df.to_excel('/Users/angelawu/Desktop/HomeOffice/HSC_unique_falseNeg.xlsx')

Sim_Ref[Sim_Ref.weight != 0][['source', 'target']].to_excel('/Users/angelawu/Desktop/HomeOffice/HSC_Ref.xlsx')

# test arboreto
network = grnboost2(expression_data=np.matrix(Sim.T),
                    gene_names=list(Sim.index))

network.to_csv('/Users/angelawu/Desktop/output2.tsv', sep='\t', index=False, header=False)
