#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 14:35:57 2019

@author: angelawu
"""



from __future__ import absolute_import

import pandas as pd
import numpy as np
import copy

#path = '/Users/angelawu/Desktop/PhD_Pro/TODOLIST/BN_Test/SC_Published_Data/'
#HSC_Data = pd.read_excel(path + 'HSC_DATA/Moignard_HSC_Data.xlsx', index_col = 0, sheetname = 0)
#HSC_Data[HSC_Data == 25] = 0
#HSC_Ref = pd.read_csv(path + 'HSC_DATA/Moignard_HSC_ReferenceEdges.tsv', sep = '\t')

#
path = '/home/mwu/MING_V9T/PhD_Pro/TODOLIST/BN_Test/SC_Published_Data/'
#path = '/Users/angelawu/GitHub/network_inference_tutorials/'
#Sim_1 = pd.read_csv(path + 'simulated_datasets/100_yeast3_large_dropouts_low.txt', sep = '\t', index_col = 0)
#Sim_1_Ref = pd.read_csv(path + 'goldstandards/100_yeast3.tsv', sep = '\t', header = -1)
#Sim_1_Ref.columns = ['node1', 'node2', 'value']
#Sim_1_Ref = Sim_1_Ref.loc[Sim_1_Ref['value'] > 0]
#
#
#Sim_2 = pd.read_csv(path + 'simulated_datasets/100_yeast3_medium.txt', sep = '\t', index_col = 0)
#Sim_2_Ref = pd.read_csv(path + 'goldstandards/100_yeast3.tsv', sep = '\t', header = -1)
#Sim_2_Ref.columns = ['node1', 'node2', 'value']
#Sim_2_Ref = Sim_2_Ref.loc[Sim_2_Ref['value'] > 0]

HSC_Data = pd.read_excel(path + 'HSC_DATA/Moignard_HSC_Data.xlsx', index_col = 0, sheetname = 0)
HSC_Data[HSC_Data == 25] = 0
HSC_Ref = pd.read_csv(path + 'HSC_DATA/Moignard_HSC_ReferenceEdges.tsv', sep = '\t')

ESC_Data = pd.read_excel(path + 'ESC_DATA/GSE59892_ESC_Dataset.xlsx', index_col = 0).T
ESC_Gene_Symbol = pd.read_csv(path + 'ESC_DATA/GSE59892_ESC_Genesymobol.txt', sep = '\t')
ESC_Gene_Symbol_dict = dict(zip(ESC_Gene_Symbol.Ref_ID, ESC_Gene_Symbol.GeneSymbol))
ESC_Ref = pd.read_csv(path + 'ESC_DATA/GSE59892_ESC_ReferenceEdges.tsv', sep = '\t')

ESC_Data.columns = list(ESC_Gene_Symbol_dict.get(i) for i in list(ESC_Data.columns))
ESC_Data.to_csv(path + 'Expr.txt', sep = '\t')
# =============================================================================
# Test with Preprocessing and docker run
# =============================================================================

import sys
sys.path.append('/home/mwu/MING_V9T/PhD_Pro/PySCNet/')
import _pickle as pk
from PySCNet.Preprocessing import gnetdata
from PySCNet.BuildNet import gne_dockercaller as gdocker
from PySCNet.NetEnrich import graph_toolkit as gt
from PySCNet.Plotting import dynamic_net as dn
import Eva_Test

#path =  '/home/mwu/MING_V9T/PhD_Pro/PySCNet/BuildNet/Docker_App/'
#Expr = pd.read_csv(path + "PIDC/100_yeast2_medium.txt", sep = '\t', header = 0, index_col = 0)
#pd.DataFrame.to_csv(Expr, path + 'tmp.txt', sep='\t')
#gne_exp = gnetdata.Gnetdata(Sim_1)
ESC_gne = gnetdata.Gnetdata(ESC_Data)
HSC_gne = gnetdata.Gnetdata(HSC_Data)

#gne_exp_PIDC = gdocker.rundocker(gne_exp.deepcopy, method = 'PIDC')
#gne_exp_GENIE3 = gdocker.rundocker(gne_exp.deepcopy, method = 'GENIE3')
#gne_exp_CORR = gdocker.rundocker(gne_exp.deepcopy, method = 'CORR')
#time_point = pd.DataFrame({'cell_id': list(Sim_1.columns), 'time_point': [1 for i in range(Sim_1.shape[1])]})
#gne_exp_SCODE = gdocker.rundocker(gne_exp.deepcopy, method = 'SCODE', time_point = time_point)

gne_ESC_PIDC = gdocker.rundocker(ESC_gne.deepcopy, method = 'PIDC')
gne_ESC_GENIE3 = gdocker.rundocker(ESC_gne.deepcopy, method = 'GENIE3')
gne_ESC_CORR = gdocker.rundocker(ESC_gne.deepcopy, method = 'CORR')
time_point = pd.DataFrame({'cell_id': list(ESC_Data.columns), 'time_point': [1 for i in range(ESC_gne.shape[1])]})
gne_ESC_SCODE = gdocker.rundocker(ESC_gne.deepcopy, method = 'SCODE', time_point = time_point)

gne_HSC_PIDC = gdocker.rundocker(HSC_gne.deepcopy, method = 'PIDC')
gne_HSC_GENIE3 = gdocker.rundocker(HSC_gne.deepcopy, method = 'GENIE3')
gne_HSC_CORR = gdocker.rundocker(HSC_gne.deepcopy, method = 'CORR')
time_point = pd.DataFrame({'cell_id': list(HSC_Data.columns), 'time_point': [1 for i in range(HSC_gne.shape[1])]})
gne_HSC_SCODE = gdocker.rundocker(HSC_gne.deepcopy, method = 'SCODE', time_point = time_point)



#tmp = [['A','B', 0.1], ['A', 'C', 0.2], ['A', 'D', 0.3]]
#tmp_2 = [['B','A', 0.1], ['A', 'C', 0.2], ['A', 'E', 0.3]]
#
#link_1 = pd.DataFrame(tmp, columns=['source', 'target', 'weight'])
#link_2 = pd.DataFrame(tmp_2, columns=['source', 'target', 'weight'])
#
#g = gt.graph_merge(link_1, link_2, method='intersection')

#ensemble test
#path = '/home/mwu/MING_V9T/PhD_Pro/Test/KK_Run36_D15_CD4/'
#gne_exp_CORR = gnetdata.load_Gnetdata_object(path + 'cluster_3_CORR.pk')
#gne_exp_GENIE3 = gnetdata.load_Gnetdata_object(path + 'cluster_3_GENIE3.pk')
#gne_exp_PIDC = gnetdata.load_Gnetdata_object(path + 'cluster_3_PIDC.pk')
#gne_exp_SCODE = gnetdata.load_Gnetdata_object(path + 'cluster_3_SCODE.pk')

links_dict = {'genie3': gne_exp_GENIE3.NetAttrs['links'], 'corr': gne_exp_CORR.NetAttrs['links'],
              'pidc': gne_exp_PIDC.NetAttrs['links'], 'scode': gne_exp_SCODE.NetAttrs['links']}

ESC_links_dict = {'genie3': gne_ESC_GENIE3.NetAttrs['links'], 'corr': gne_ESC_CORR.NetAttrs['links'],
              'pidc': gne_ESC_PIDC.NetAttrs['links'], 'scode': gne_ESC_SCODE.NetAttrs['links']}

HSC_links_dict = {'genie3': gne_HSC_GENIE3.NetAttrs['links'], 'corr': gne_HSC_CORR.NetAttrs['links'],
              'pidc': gne_HSC_PIDC.NetAttrs['links'], 'scode': gne_HSC_SCODE.NetAttrs['links']}




Eclass = gt.ensemble_classifier(links_dict, threshold=0.5, num_trees=3000)
Eclass_pos = Eclass.sort_values('connected', ascending = False).head(500)
Eva_Test.test_run(links = Eclass, Ref_links=Sim_1_Ref, input_dataset = Sim_1, filename='ensemble')


Eclass = gt.ensemble_classifier(HSC_links_dict, threshold=0.5, num_trees=1000)
Eclass_pos = Eclass.sort_values('connected', ascending = False).head(500)
Eva_Test.test_run(links = Eclass_pos, Ref_links=HSC_Ref, input_dataset = HSC_Data, filename='ensemble')








links = copy.deepcopy(gne_exp_GENIE3.NetAttrs['links'])
links.columns = ['source', 'target', 'connected']
Eva_Test.test_run(links = links,
                  Ref_links=Sim_1_Ref, input_dataset = Sim_1, filename='Genie3')

Bclass = gt.bnn_classifier(links_dict, threshold=0.5)
Bclass_pos = Bclass.sort_values('connected', ascending = False).head(500)
Eva_Test.test_run(links = Bclass_pos, Ref_links=Sim_1_Ref, input_dataset = Sim_1, filename='test')





