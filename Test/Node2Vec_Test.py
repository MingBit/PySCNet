#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:51:43 2019

@author: mwu
"""

from __future__ import absolute_import

import pandas as pd
import numpy as np
import sys
sys.path.append('/home/mwu/MING_V9T/PhD_Pro/PySCNet/')
import _pickle as pk
from PySCNet.Preprocessing import gnetdata
from PySCNet.Preprocessing import general_pipeline as pipeline
from PySCNet.BuildNet import gne_dockercaller as gdocker
from PySCNet.NetEnrich import graph_toolkit as gt
from PySCNet.Plotting import show_net as sn
import Eva_Test
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from PySCNet.BuildNet.Models import my_node2vec
from PySCNet.BuildNet.Models import model_node2vec


def test_node2vec(Expr, Ref, p, q, dim_list, walk_list, num_walks_list):

    links_list = dict()
    
    for dimensions in dim_list:
        for walk_length in walk_list:
            for num_walks in num_walks_list:
                
                links = my_node2vec.my_node2vec(Expr, p = p, q = q, walk_len=walk_length, 
                                            num_walks=num_walks, size=dimensions)
    
                links = gdocker._remove_duplicate(links)
                links_list['dims_' + str(dimensions) + 'walk_len_' + str(walk_length) + 'num_walks_' + str(num_walks)] = links
    
                
    return(links_list)
                
    
def test_thresholds(links_list, Expr, Ref, threshold, filepath):
    
    for key in links_list.keys():
        links = links_list[key]
        DL_fpr = dict()
        DL_tpr = dict()
        DL_auc = dict()
        
        tmp_fpr, tmp_tpr, tmp_auc = Eva_Test.test_run(links = links[links['weight'] > np.quantile(links['weight'], threshold)],
                                           Ref_links = Ref, input_dataset = Expr)
                
        DL_fpr[key] = tmp_fpr
        DL_tpr[key] = tmp_tpr
        DL_auc[key] = tmp_auc
                
    colors = sns.color_palette().as_hex()[1:len(links_list)]
    plt.figure(figsize = (8,8))
    
    for i, color in zip(list(DL_auc.keys()), colors):
            plt.plot(DL_fpr[i], DL_tpr[i], color = color, label = 'ROC curve of method {0} (area = {1:0.2f})'.format(i, DL_auc[i]))
    
    
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate', fontsize = 20)
    plt.ylabel('True Positive Rate', fontsize = 20)
    plt.title('Algorithms performance', fontsize = 25)
    plt.legend(loc="lower right")
    
    plt.savefig(filepath)
    plt.close()
                


                
                
                
                
                