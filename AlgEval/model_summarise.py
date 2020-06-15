#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 14:16:47 2019

@author: mwu
"""
from __future__ import absolute_import
import sys
import pandas as pd
import numpy as np
import os
import glob
import ntpath
import seaborn as sns
import matplotlib.pyplot as plt

from AlgEval import basic_functions as bf
from AlgEval import Eva_Test
sys.path.append(os.getenv("HOME") + '/MING_V9T/PhD_Pro/PySCNet/')
from PySCNet.Preprocessing import gnetdata


input_path = os.getenv('HOME') + '/MING_V9T/PhD_Pro/Test/Simulation/Node2Vec_BEELINE_Data_Res/'
path = os.getenv("HOME") + '/MING_V9T/PhD_Pro/Test/Simulation//BEELINE-data/inputs/Curated/HSC/'
filename_dict = {'HSC-2000-10': '0', 'HSC-2000-10-50': '0.5',
                 'HSC-2000-10-70': '0.7'}

dataSum = pd.DataFrame({'Data': np.repeat(['HSC-2000-10', 'HSC-2000-10-50',
                                           'HSC-2000-10-70'], 3),
                        'SC_drop': np.repeat(['0', '0.5', '0.7'], 3),
                        'Algorithms': ['GENIE3', 'PIDC', 'CORR'] * 3,
                        # G30_C1000
                        # 'AUROC': [0.56, 0.55, 0.52, 0.55, 0.56, 0.52, 0.54, 0.55, 0.51],
                        # 'AUPRC': [0.29, 0.30, 0.25, 0.27, 0.31, 0.24, 0.26, 0.28, 0.23]})
                        # #G30_C3000
                        # 'AUROC': [0.58, 0.57, 0.50, 0.55, 0.56, 0.51, 0.53, 0.55, 0.51],
                        # 'AUPRC': [0.30, 0.32, 0.23, 0.27, 0.32, 0.23, 0.25, 0.28, 0.23]})
                        # GSD-2000-1
                        # 'AUROC': [0.57, 0.57, 0.51, 0.56, 0.57, 0.51, 0.55, 0.52, 0.49],
                        # 'AUPRC': [0.29, 0.30, 0.25, 0.29, 0.31, 0.25, 0.29, 0.26, 0.23]})
                        # GSD-2000-10
                        # 'AUROC': [0.59, 0.57, 0.53, 0.56, 0.57, 0.51, 0.55, 0.52, 0.49],
                        # 'AUPRC': [0.30, 0.30, 0.25, 0.29, 0.31, 0.25, 0.29, 0.26, 0.23]})
                        # HSC-2000-1
                        # 'AUROC': [0.64, 0.64, 0.56, 0.62, 0.63, 0.55, 0.63, 0.62, 0.55],
                        # 'AUPRC': [0.45, 0.47, 0.37, 0.41, 0.44, 0.38, 0.39, 0.41, 0.35]})
                        # HSC-2000-10
                        'AUROC': [0.66, 0.63, 0.57, 0.65, 0.63, 0.55, 0.59, 0.57, 0.54],
                        'AUPRC': [0.44, 0.44, 0.35, 0.43, 0.43, 0.34, 0.37, 0.36, 0.32]})



for file in filename_dict.keys():

    print(file)
    Sim, Sim_Ref, Sim_Ref_dropout, top_edges = bf.read_data(path + file)
    SC_drop = filename_dict[file]
    file_list = glob.glob(input_path + file + '_PCA/*dropRate_0.75*')

    node_dict_list = list()
    for sub_file in file_list:
        gne_Sim_node2vec = gnetdata.load_Gnetdata_object(sub_file)
        top_links = gne_Sim_node2vec.NetAttrs['links'].reindex(
            gne_Sim_node2vec.NetAttrs['links'].weight.abs().sort_values(ascending=False).index).head(top_edges)

        node_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'],
                             Eva_Test.test_run(links=top_links, Ref_links=Sim_Ref, input_dataset=Sim)))

        node_dict_list.append(node_dict)

    node_auc_avg = np.mean(list(node_dict_list[i]['auc'] for i in range(len(node_dict_list))))
    node_pre_avg = np.mean(list(node_dict_list[i]['avg_pre'] for i in range(len(node_dict_list))))
    dataSum = dataSum.append(pd.Series([file, SC_drop, 'Node2Vec', node_auc_avg, node_pre_avg], index=dataSum.columns),
                             ignore_index=True)

    # node_precision_list = np.mean(list(node_dict_list[i]['precision'] for i in range(len(node_dict_list))))
dataSum.to_csv(input_path + 'PCA_ROC_Res/HSC-2000-10_Algorithms_Vs_Node2vec_0.75.csv')
# colors = ['#C25E6E', '#B4E3DA', '#ECAD7C', '#5AAA9A', '#7779A8', '#A0BBCB', '#E6BEC5', '#6C6B74']
fig, axs = plt.subplots(ncols=2, figsize=(12, 6))
sns.barplot(y='AUROC', x='Algorithms', hue='SC_drop', data=dataSum, ax=axs[0], palette="Set2").set(title='AUROC values')
sns.barplot(y='AUPRC', x='Algorithms', hue='SC_drop', data=dataSum, ax=axs[1], palette="Set2").set(title='AUPRC values')

fig.savefig(input_path + 'PCA_ROC_Res/HSC-2000-10_Algorithms_Vs_Node2vec_0.75.pdf')

##########parameters boxplot###############
file_stat = list()
for file in filename_dict.keys():
    print(input_path + file + '_PCA/')
    Sim, Sim_Ref, Sim_Ref_dropout, top_edges = bf.read_data(path + file)
    SC_drop = filename_dict[file]
    file_list = glob.glob(input_path + file + '_PCA/*')

    for sub_file in file_list:
        filename = ntpath.basename(sub_file)
        data_name = filename.split('_node2vec_')[0]
        pq_state = filename.split('_node2vec_')[1][:11]
        dropState = filename.split('_walk_5_')[1].split('_repeat_')[0]

        gne_Sim_node2vec = gnetdata.load_Gnetdata_object(sub_file)
        top_links = gne_Sim_node2vec.NetAttrs['links'].reindex(
            gne_Sim_node2vec.NetAttrs['links'].weight.abs().sort_values(ascending=False).index).head(top_edges)

        node_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'],
                             Eva_Test.test_run(links=top_links, Ref_links=Sim_Ref, input_dataset=Sim)))

        file_stat.append([data_name, pq_state, dropState, node_dict['auc'], node_dict['avg_pre']])

file_stat_df = pd.DataFrame(file_stat)
file_stat_df.columns = ['data_name', 'pq_value', 'dropRate', 'AUROC', 'AUPRC']
file_stat_df.to_csv(input_path + 'PCA_ROC_Res/HSC-2000-10_Node2vec_Parameters_Summary.csv')

one_file = list(filename_dict.keys())[0]
order = ["dropRate_0.25", "dropRate_0.5", "dropRate_0.75"]

fig, axs = plt.subplots(ncols=2, figsize=(12, 6))
ax1 = sns.stripplot(y='AUPRC', x='dropRate',
                    data=file_stat_df[file_stat_df.data_name == one_file],
                    marker='o',
                    alpha=0.5,
                    hue='pq_value',
                    color='grey',
                    order=order,
                    ax=axs[0])
sns.boxplot(y='AUPRC', x='dropRate', hue='pq_value',
            data=file_stat_df[file_stat_df.data_name == one_file],
            ax=axs[0], order=order, palette='Set1').set(title=one_file)

handles, labels = ax1.get_legend_handles_labels()
# # specify just one legend
l = ax1.legend(handles[0:2], labels[0:2])

ax2 = sns.stripplot(y='AUROC', x='dropRate',
                    data=file_stat_df[file_stat_df.data_name == one_file],
                    marker='o',
                    alpha=0.5,
                    hue='pq_value',
                    color='grey',
                    order=order,
                    ax=axs[1])

sns.boxplot(y='AUROC', x='dropRate', hue='pq_value',
            data=file_stat_df[file_stat_df.data_name == one_file],
            ax=axs[1], order=order, palette='Set1').set(title=one_file)
# # get legend information from the plot object
l = ax2.legend(handles[0:2], labels[0:2])
# ax1.get_legend().remove()
fig.savefig(input_path + 'PCA_ROC_Res/' + one_file + '_Node2vec_Parameters.pdf')

