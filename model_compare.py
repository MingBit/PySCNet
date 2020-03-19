#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 14:16:47 2019

@author: mwu
"""
from __future__ import absolute_import

import pandas as pd
import numpy as np
import sys
import itertools
import os
sys.path.append(os.getenv("HOME") + '/MING_V9T/PhD_Pro/PySCNet/')
import _pickle as pk
from PySCNet.Preprocessing import gnetdata
from PySCNet.Preprocessing import general_pipeline as pipeline
from PySCNet.BuildNet import gne_dockercaller as gdocker
#from PySCNet.NetEnrich import graph_toolkit as gt
from PySCNet.Plotting import show_net as sn
import Eva_Test
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import pdb
from PySCNet.BuildNet import gne_modelcaller as gmodel


path = os.getenv("HOME") + '/MING_V9T/PhD_Pro/Test/Simulation/BoolODE_Data/'

filename = 'G30_C3000_Clu_6_Tra_4'
p = 0.1
q = 1
topology = 'random'

Sim = pd.read_csv(path + filename + '/ExpressionData.csv', sep = ',', index_col = 0)
Sim_Ref_raw = pd.read_csv(path + filename + '/refNetwork.csv', sep = ',', index_col = None)
print(Sim_Ref_raw.shape)

Sim_Ref_raw['weight'] = Sim_Ref_raw['Type'].map({'+': 1, '-': -1})
Sim_Ref_raw = Sim_Ref_raw.drop(columns = 'Type')
Sim_Ref_raw.columns = ['source', 'target', 'weight']

#Full_Ref = pd.DataFrame([p for p in itertools.product(Sim.index, repeat=2)], columns = ['source', 'target'])
Full_Ref = pd.DataFrame(itertools.permutations(Sim.index, 2), columns = ['source', 'target'])
Sim_Ref = pd.merge(Full_Ref, Sim_Ref_raw, how = 'outer').fillna(0)

# Sim_Ref = Sim_Ref[Sim_Ref.node1 != Sim_Ref.node2].reset_index()

print(Sim_Ref.shape)

#Sim_Top = pd.read_csv(path + filename.split('_')[0] + '.txt', sep = '\t', index_col = None)
#Sim_Top = Sim_Top[Sim_Top.Topology != topology]['Gene']
#
#Sim_Ref = Sim_Ref[Sim_Ref.node1.isin(Sim_Top)]
#Sim = Sim[Sim.index.isin(Sim_Top)]
#print(Sim.shape)

Sim_gne = gnetdata.Gnetdata(Sim)

gne_Sim_PIDC = gdocker.rundocker(Sim_gne.deepcopy, method = 'PIDC')
gne_Sim_GENIE3 = gdocker.rundocker(Sim_gne.deepcopy, method = 'GENIE3')
gne_Sim_CORR = gdocker.rundocker(Sim_gne.deepcopy, method = 'CORR')

#gne_Sim_node2vec = gmodel.call_node2vec(Sim_gne.deepcopy, 
#                                         p = 1, q = 0.1, dim_list = [5], walk_list = [5],
#                                         num_walks_list = [1000], workers = 10)

#choose top head(Sim_Ref.shape[0])
top_edges = Sim_Ref_raw.shape[0]

PIDC_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'], 
                     Eva_Test.test_run(links = gne_Sim_PIDC.NetAttrs['links'].reindex(gne_Sim_PIDC.NetAttrs['links'].weight.abs().sort_values(ascending = False).index).head(top_edges),
                                       Ref_links=Sim_Ref, input_dataset = Sim)))

GENIE3_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'], 
                       Eva_Test.test_run(links = gne_Sim_GENIE3.NetAttrs['links'].reindex(gne_Sim_GENIE3.NetAttrs['links'].weight.abs().sort_values(ascending = False).index).head(top_edges),
                                         Ref_links=Sim_Ref, input_dataset = Sim)))

CORR_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'], 
                     Eva_Test.test_run(links = gne_Sim_CORR.NetAttrs['links'].reindex(gne_Sim_CORR.NetAttrs['links'].weight.abs().sort_values(ascending = False).index).head(top_edges),
                                       Ref_links=Sim_Ref, input_dataset = Sim)))

def build_curves(ax, node_dict_list, GENIE3_dict, PIDC_dict, CORR_dict, 
                 curve, p, q, filename, output_path):

    keywords = ['fpr', 'tpr'] if curve == 'ROC' else ['recall_list', 'pre']
    node_fpr_pre_df = pd.DataFrame(node_dict_list[i][keywords[0]] for i in range(len(node_dict_list))).fillna(1)
    node_tpr_recall_df = pd.DataFrame(node_dict_list[i][keywords[1]] for i in range(len(node_dict_list))).fillna(1)
    node_auc_list = list(node_dict_list[i]['auc'] for i in range(len(node_dict_list)))
    node_avgpre_list = list(node_dict_list[i]['avg_pre'] for i in range(len(node_dict_list)))
    
    colors = sns.color_palette().as_hex() + sns.color_palette('hls', 8).as_hex()
    
    for i in range(len(node_dict_list)):
        
        ax.plot(node_dict_list[i][keywords[0]], node_dict_list[i][keywords[1]], color = 'grey', alpha = 0.2)    
    
    fpr_dict = {'GENIE3': GENIE3_dict[keywords[0]], 'PIDC': PIDC_dict[keywords[0]], 'CORR': CORR_dict[keywords[0]], 'Node2Vec': node_fpr_pre_df.mean(0)}
    tpr_dict = {'GENIE3': GENIE3_dict[keywords[1]], 'PIDC': PIDC_dict[keywords[1]], 'CORR': CORR_dict[keywords[1]], 'Node2Vec': node_tpr_recall_df.mean(0)}
    auc_dict = {'GENIE3': GENIE3_dict['auc'], 'PIDC': PIDC_dict['auc'], 'CORR': CORR_dict['auc'], 'Node2Vec': np.mean(node_auc_list)}
    avgpre_dict = {'GENIE3': GENIE3_dict['avg_pre'], 'PIDC': PIDC_dict['avg_pre'], 'CORR': CORR_dict['avg_pre'], 'Node2Vec': np.mean(node_avgpre_list)}
    
    for i, color in zip(list(auc_dict.keys()), colors):
            ax.plot(fpr_dict[i], tpr_dict[i], 
                     label='ROC curve {0} (area = {1:0.2f})'.format(i, auc_dict[i]) if curve == 'ROC' else 
                     'PR curve {0} (area = {1:0.2f})'.format(i, avgpre_dict[i]),
                     color = color)
    if curve == 'ROC':
        
        ax.plot([0, 1], [0, 1], 'k--') 
        
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.0])
    ax.set_xlabel('False Positive Rate', fontsize = 20)
    ax.set_ylabel('True Positive Rate', fontsize = 20)
    # ax.set_title(filename + '_p_' + str(p) + '_q_' + str(q) + ':Algorithms performance', fontsize = 12)
    ax.legend(loc="lower right")




def build_boxplot(ax, node_dict_list, GENIE3_dict, PIDC_dict, CORR_dict):
    node_precision_list = list(node_dict_list[i]['precision'] for i in range(len(node_dict_list)))
    node_recall_list = list(node_dict_list[i]['recall'] for i in range(len(node_dict_list)))
    node_f1score_list = list(node_dict_list[i]['f1_score'] for i in range(len(node_dict_list)))
    
    Algorithms = ['GENIE3', 'PIDC', 'CORR', 'Node2Vec']
    Eva_Methods = ['Precision', 'Recall', 'F1_Score']
    
    data = [[GENIE3_dict['precision'], PIDC_dict['precision'], CORR_dict['precision'], np.mean(node_precision_list)],
            [GENIE3_dict['recall'], PIDC_dict['recall'], CORR_dict['recall'], np.mean(node_recall_list)],
            [GENIE3_dict['f1_score'], PIDC_dict['f1_score'], CORR_dict['f1_score'], np.mean(node_f1score_list)]]
    
    colors = sns.color_palette().as_hex() + sns.color_palette('hls', 8).as_hex()
    # fig, ax = plt.subplots(figsize = (12,8))    
    X = np.arange(4)
    # ax = fig.add_axes([0,0,1,1])
    for i in range(len(data)):
        ax.bar(X + 0.25*i, data[i], color = colors[i], width = 0.25)    
    
    ax.set_xlabel('Algorithms', fontsize=25)
    ax.set_ylabel('Performance', fontsize = 25)
    ax.set_xticks(X + 0.25)
    ax.set_xticklabels(Algorithms, fontsize = 20)
    ax.legend(Eva_Methods, loc = 'upper left')

    
    

input_path =  os.getenv("HOME") + '/MING_V9T/PhD_Pro/Test/Simulation/Node2Vec_BoolODE_Res/'
output_path = os.getenv("HOME") + '/MING_V9T/PhD_Pro/Test/Simulation/Node2Vec_BoolODE_Res/ROC_Res/'
    
node_dict_list = list()
for i in range(10):
    gne_Sim_node2vec = gnetdata.load_Gnetdata_object(input_path + 
                                                    filename + '/'+ filename + '_node2vec_p_' + str(p) + '_q_' + str(q) + '_dim_5_walk_5_repeat_' + str(i+1) + '.pk')
    top_links = gne_Sim_node2vec.NetAttrs['links'].reindex(gne_Sim_node2vec.NetAttrs['links'].weight.abs().sort_values(ascending = False).index).head(top_edges)
       
    node_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'], 
                         Eva_Test.test_run(links = top_links, Ref_links=Sim_Ref, input_dataset = Sim)))
    
    node_dict_list.append(node_dict)

    

fig = plt.figure(figsize = (20,20))
grid = plt.GridSpec(2, 2, wspace=0.2, hspace=0.2)
ax1 = fig.add_subplot(grid[0,0])
ax2 = fig.add_subplot(grid[0,1])
ax3 = fig.add_subplot(grid[1,0:])
fig.suptitle(filename + '_p_' + str(p) + '_q_' + str(q) + ':Algorithms performance', fontsize=30)
build_curves(ax1, node_dict_list, GENIE3_dict, PIDC_dict, CORR_dict,'ROC', p, q, filename, output_path)
build_curves(ax2, node_dict_list, GENIE3_dict, PIDC_dict, CORR_dict,'PCR', p, q, filename, output_path)
build_boxplot(ax3, node_dict_list, GENIE3_dict, PIDC_dict, CORR_dict)

fig.savefig(output_path + filename +'_p_' + str(p) + '_q_' + str(q) + '_' + '_topTP.pdf')