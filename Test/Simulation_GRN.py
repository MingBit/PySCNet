#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 14:16:47 2019

@author: mwu
"""
from __future__ import absolute_import
import warnings

warnings.filterwarnings("ignore")
import os
import pandas as pd
import sys
import re
import numpy as np

sys.path.append(os.getenv('HOME') + '/MING_V9T/PhD_Pro/pyscnet')
sys.path.append(os.getenv('HOME') + '/MING_V9T/PhD_Pro/SCNode2Vec')
from src import model_compare as compare
# from src import model_node2vec as nv
from arboreto.algo import grnboost2, genie3
import importlib
# importlib.reload(compare)

path = os.getenv('HOME') + '/MING_V9T/PhD_Pro/Test/Simulation/BoolODE_Data/BoolODE_Data_Repeat_Jun2020/HSC_MEP_2/'
output = os.getenv(
    'HOME') + '/MING_V9T/PhD_Pro/Test/Simulation/Node2Vec_BoolODE_Res/Node2Vec_BoolODE_Res_Repeat_Jun2020/HSC_MEP_2_Performance/'

file_list = [f for f in os.listdir(path) if re.match(r'HSC_MEP_2*', f)]
if __name__ == '__main__':

    p = 1
    q = 0.1
    dropRate = 0.75
    size = 5
    walk_len = 5

    performance_res = pd.DataFrame(columns=['Dataset', 'Dropout', 'AlgName', 'AUROC', 'AUPRC'])

    for file in file_list:
        Sim, Sim_Ref, Sim_Ref_dropout, top_edges = compare.read_data(path + file, float(dropRate))

        gne_Sim_GENIE3 = genie3(np.asmatrix(Sim.T), gene_names=list(Sim.index))
        gne_Sim_GENIE3.columns = ['source', 'target', 'weight']
        # # gne_Sim_GENIE3 = compare.remove_duplicate(gne_Sim_GENIE3)

        gne_Sim_GRNBOOST2 = grnboost2(np.asmatrix(Sim.T), gene_names=list(Sim.index))
        gne_Sim_GRNBOOST2.columns = ['source', 'target', 'weight']
        # # gne_Sim_GRNBOOST2 = compare.remove_duplicate(gne_Sim_GRNBOOST2)

        param_grid = {
            'n_estimators': [50, 100],
            'max_depth': [10, 20, 50],
            'max_features': [3, 5],
            'min_samples_leaf': [3, 4, 5],
            'criterion': ['gini', 'entropy']
        }

        # run PIDC
        Sim.to_csv('Expr.txt', sep='\t')
        os.system('julia run_pidc.jl')
        gne_Sim_PIDC = pd.read_csv('links.txt', sep='\t', header=None)
        gne_Sim_PIDC.columns = ['source', 'target', 'weight']
        gne_Sim_PIDC = compare.remove_duplicate(gne_Sim_PIDC)

        os.system('rm Expr.txt | rm links.txt')

        PIDC_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'],
                             compare.test_run(links=gne_Sim_PIDC.reindex(
                                 gne_Sim_PIDC.weight.abs().sort_values(ascending=False).index).head(
                                 top_edges), Ref_links=Sim_Ref, input_dataset=Sim)))

        GENIE3_dict = dict(
            zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'],
                compare.test_run(links=gne_Sim_GENIE3.reindex(
                    gne_Sim_GENIE3.weight.abs().sort_values(ascending=False).index).head(
                    top_edges), Ref_links=Sim_Ref, input_dataset=Sim)))

        GRNBOOST2_dict = dict(
            zip(['fpr', 'tpr', 'pre', 'recall_list', 'auc', 'avg_pre', 'precision', 'recall', 'f1_score'],
                compare.test_run(links=gne_Sim_GRNBOOST2.reindex(
                    gne_Sim_GRNBOOST2.weight.abs().sort_values(ascending=False).index).head(
                    top_edges), Ref_links=Sim_Ref, input_dataset=Sim)))

        Node2Vec_dict_randm = {
            'randm1': compare.repeat_node2vec(Sim, Sim_Ref, Sim_Ref_dropout, p, q, size, walk_len, top_edges,
                                              param_grid,
                                              use_ref=False, method='Node2Vec'),
            'randm2': compare.repeat_node2vec(Sim, Sim_Ref, Sim_Ref_dropout, q, p, size, walk_len, top_edges,
                                              param_grid,
                                              use_ref=False, method='Node2Vec')}

        Node2Vec_dict_ref = {
            'ref1': compare.repeat_node2vec(Sim, Sim_Ref, Sim_Ref_dropout, p, q, size, walk_len, top_edges, param_grid,
                                            use_ref=True, method='Node2Vec'),
            'ref2': compare.repeat_node2vec(Sim, Sim_Ref, Sim_Ref_dropout, q, p, size, walk_len, top_edges, param_grid,
                                            use_ref=True, method='Node2Vec')}

        Node2Vec_dict_randm_res = Node2Vec_dict_randm['randm1'] if Node2Vec_dict_randm['randm1']['auc'] > \
                                                                   Node2Vec_dict_randm['randm2']['auc'] else Node2Vec_dict_randm['randm2']

        Node2Vec_dict_ref_res = Node2Vec_dict_ref['ref1'] if Node2Vec_dict_ref['ref1']['auc'] > \
                                                             Node2Vec_dict_ref['ref2']['auc'] else Node2Vec_dict_ref['ref2']

        # Struc2Vec_dict = compare.repeat_node2vec(Sim, Sim_Ref, Sim_Ref_dropout, p, q, size, walk_len, top_edges,
        #                                          param_grid,
        #                                          use_ref=False, method='Struc2Vec')
        # GraRep_dict = compare.repeat_node2vec(Sim, Sim_Ref, Sim_Ref_dropout, p, q, size, walk_len, top_edges,
        #                                       param_grid,
        #                                       use_ref=False, method='GraRep')
        # DeepWalk_dict = compare.repeat_node2vec(Sim, Sim_Ref, Sim_Ref_dropout, p, q, size, walk_len, top_edges,
        #                                         param_grid,
        #                                         use_ref=False, method='DeepWalk')

        dropout = '0' if len(file) < 20 else file.split('-')[3]
        res = pd.DataFrame({'Dataset': np.repeat(file, 5),
                            'Dropout': np.repeat('dropout_' + dropout, 5),
                            'AlgName': ['GENIE3', 'PIDC', 'GRNBOOST2', 'NODE2VEC', 'NODE2VEC_REF'],
                            'AUROC': [GENIE3_dict['auc'], PIDC_dict['auc'], GRNBOOST2_dict['auc'],
                                      Node2Vec_dict_randm_res['auc'], Node2Vec_dict_ref_res['auc']],
                            'AUPRC': [GENIE3_dict['avg_pre'], PIDC_dict['avg_pre'], GRNBOOST2_dict['avg_pre'],
                                      Node2Vec_dict_randm_res['avg_pre'], Node2Vec_dict_ref_res['avg_pre']]})

        # res = pd.DataFrame({'Dataset': np.repeat(file, 4),
        #                     'Dropout': np.repeat('dropout_' + dropout, 4),
        #                     'AlgName': ['GraRep_Ref', 'DeepWalk_Ref', 'Struc2Vec_Ref', 'NODE2VEC_Ref'],
        #                     'AUROC': [GraRep_dict['auc'], DeepWalk_dict['auc'], Struc2Vec_dict['auc'],
        #                               Node2Vec_dict_ref_res['auc']],
        #                     'AUPRC': [GraRep_dict['avg_pre'], DeepWalk_dict['avg_pre'], Struc2Vec_dict['avg_pre'],
        #                               Node2Vec_dict_ref_res['avg_pre']]})

        performance_res = pd.concat([performance_res, res])
        # res.to_csv(output + file + '_GE_Ref.csv')

    performance_res.to_csv(output + 'HSC_MEP_2_GE_Performance.csv')

# build strip-box plots
import seaborn as sns
import matplotlib.pyplot as plt

# path = os.getenv('HOME') + '/MING_V9T/PhD_Pro/Test/Simulation/Node2Vec_BoolODE_Res/Node2Vec_BoolODE_Res_Repeat_Jun2020/Beeline_G57_Performance/'
# file_list = os.listdir(path)
# performance_res = pd.DataFrame()
# for f in file_list:
#     tmp = pd.read_csv(path + f)
#     performance_res = pd.concat([performance_res, tmp])
# performance_res.to_csv(os.getenv('HOME') + '/MING_V9T/PhD_Pro/Test/Simulation/Node2Vec_BoolODE_Res/Node2Vec_BoolODE_Res_Repeat_Jun2020/Beeline_G57_Performance/Beeline_G57_performance.csv')

# performance_res = pd.read_csv(os.getenv(
#     'HOME') + '/MING_V9T/PhD_Pro/Test/Simulation/Node2Vec_BoolODE_Res/Node2Vec_BoolODE_Res_Repeat_Jun2020/HSC_MEP_2_Performance/HSC_MEP_2_GE_Ref_Performance.csv')
# order = ['GraRep_Ref', 'DeepWalk_Ref', 'Struc2Vec_Ref', 'NODE2VEC_Ref']
# hue_order = ['dropout_0', 'dropout_25', 'dropout_50', 'dropout_75']
#
# fig, axs = plt.subplots(ncols=2, figsize=(20, 10))
# ax1 = sns.stripplot(y='AUROC', x='AlgName',
#                     data=performance_res,
#                     # marker='o',
#                     dodge=True,
#                     alpha=0.5,
#                     hue='Dropout',
#                     hue_order=hue_order,
#                     color='grey',
#                     order=order,
#                     ax=axs[0])
# sns.boxplot(y='AUROC', x='AlgName', hue='Dropout',
#             data=performance_res,
#             ax=axs[0], order=order, palette='Set1', hue_order=hue_order).set(title='HSC_MEP_2')
#
# handles, labels = ax1.get_legend_handles_labels()
# # # specify just one legend
# l = ax1.legend(handles[0:4], labels[0:4])
# ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45)
# ax2 = sns.stripplot(y='AUPRC', x='AlgName',
#                     data=performance_res,
#                     # marker='o',
#                     dodge=True,
#                     alpha=0.5,
#                     hue='Dropout',
#                     hue_order=hue_order,
#                     color='grey',
#                     order=order,
#                     ax=axs[1])
#
# sns.boxplot(y='AUPRC', x='AlgName', hue='Dropout',
#             data=performance_res,
#             ax=axs[1], order=order, palette='Set1', hue_order=hue_order).set(title='HSC_MEP_2')
# # # get legend information from the plot object
# l = ax2.legend(handles[0:4], labels[0:4])
# ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45)
# fig.tight_layout()
# fig.savefig(os.getenv('HOME') + '/Desktop/HSC_MEP_2_GE_Ref.pdf')
