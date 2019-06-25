#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 14:35:57 2019

@author: angelawu
"""



from __future__ import absolute_import

import pandas as pd
import numpy as np
from rpy2.robjects import pandas2ri

pandas2ri.activate()

#path = '/Users/angelawu/Desktop/PhD_Pro/TODOLIST/BN_Test/SC_Published_Data/'
#HSC_Data = pd.read_excel(path + 'HSC_DATA/Moignard_HSC_Data.xlsx', index_col = 0, sheetname = 0)
#HSC_Data[HSC_Data == 25] = 0
#HSC_Ref = pd.read_csv(path + 'HSC_DATA/Moignard_HSC_ReferenceEdges.tsv', sep = '\t')

#
#path = '/home/mwu/MING_V9T/PhD_Pro/network_inference_tutorials/'
path = '/Users/angelawu/GitHub/network_inference_tutorials/'
Sim_1 = pd.read_csv(path + 'simulated_datasets/100_yeast3_large_dropouts_low.txt', sep = '\t', index_col = 0)
Sim_1_Ref = pd.read_csv(path + 'goldstandards/100_yeast3.tsv', sep = '\t', header = -1)
Sim_1_Ref.columns = ['node1', 'node2', 'value']
Sim_1_Ref = Sim_1_Ref.loc[Sim_1_Ref['value'] > 0]
#

Sim_2 = pd.read_csv(path + 'simulated_datasets/100_yeast3_medium.txt', sep = '\t', index_col = 0)
Sim_2_Ref = pd.read_csv(path + 'goldstandards/100_yeast3.tsv', sep = '\t', header = -1)
Sim_2_Ref.columns = ['node1', 'node2', 'value']
Sim_2_Ref = Sim_2_Ref.loc[Sim_2_Ref['value'] > 0]



# =============================================================================
# Test with Preprocessing and docker run
# =============================================================================

import sys
sys.path.append('/home/mwu/MING_V9T/PhD_Pro/PySCNet/')

from Preprocessing import gnetdata
from BuildNet import gne_dockercaller as gdocker
import pandas as pd
from NetEnrich import graph_toolkit as gt
from Plotting import dynamic_net as dn


path =  '/home/mwu/MING_V9T/PhD_Pro/PySCNet/BuildNet/Docker_App/'
Expr = pd.read_csv(path + "PIDC/100_yeast2_medium.txt", sep = '\t', header = 0, index_col = 0)
pd.DataFrame.to_csv(Expr, path + 'tmp.txt', sep='\t')
gne_exp = gnetdata.Gnetdata(Expr)

gne_exp_1 = gdocker.rundocker(gne_exp, method='PIDC')
gne_exp_1 = gdocker.buildnet(gne_exp, threshold=0.002)

gne_exp_2 = gdocker.rundocker(gne_exp, method='GENIE3')
gne_exp_2 = gdocker.buildnet(gne_exp, threshold=0.001)

gne_exp_1 = gt.get_centrality(gne_exp_1)
gne_exp_2 = gt.get_centrality(gne_exp_2)

gne_exp_1 = gt.community_detect(gne_exp_1)
gne_exp_2 = gt.community_detect(gne_exp_2)

merge_gne = gt.graph_merge(gne_exp_1.NetAttrs['links'], gne_exp_2.NetAttrs['links'], method='knn')
dfs_path_1 = gt.graph_traveral(graph = gne_exp_1.NetAttrs['graph'], start='G53', threshold=4, method='dfs')
dfs_path_2 = gt.graph_traveral(graph = gne_exp_2.NetAttrs['graph'], start='G53', threshold=4, method='dfs')

random_walk_1 = gt.random_walk(gnetdata=gne_exp_1, start='G10', supervisedby='degree', steps=10)
random_walk_2 = gt.random_walk(gnetdata=gne_exp_2, start='G10', supervisedby='degree', steps=10)

path_merge = gt.path_merge(random_walk_1, random_walk_2, 5)

dn.dynamic_netShow(gnetdata=gne_exp_2, filepath= path + 'test.html')

tmp = [['A','B', 0.1], ['A', 'C', 0.2], ['A', 'D', 0.3]]
tmp_2 = [['B','A', 0.1], ['A', 'C', 0.2], ['A', 'E', 0.3]]

link_1 = pd.DataFrame(tmp, columns=['source', 'target', 'weight'])
link_2 = pd.DataFrame(tmp_2, columns=['source', 'target', 'weight'])

g = gt.graph_merge(link_1, link_2, method='intersection')