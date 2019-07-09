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
import sys
sys.path.append('/home/mwu/MING_V9T/PhD_Pro/PySCNet/')
import _pickle as pk
from PySCNet.Preprocessing import gnetdata
from PySCNet.BuildNet import gne_dockercaller as gdocker
from PySCNet.NetEnrich import graph_toolkit as gt
from PySCNet.Plotting import show_net as sn
import Eva_Test
import copy
import matplotlib.pyplot as plt

#path = '/Users/angelawu/Desktop/PhD_Pro/TODOLIST/BN_Test/SC_Published_Data/'
#HSC_Data = pd.read_excel(path + 'HSC_DATA/Moignard_HSC_Data.xlsx', index_col = 0, sheetname = 0)
#HSC_Data[HSC_Data == 25] = 0
#HSC_Ref = pd.read_csv(path + 'HSC_DATA/Moignard_HSC_ReferenceEdges.tsv', sep = '\t')

#
path = '/home/mwu/MING_V9T/PhD_Pro/TODOLIST/BN_Test/SC_Published_Data/'

Sim_1 = pd.read_csv(path + 'simulated_datasets/100_yeast3_large_dropouts_low.txt', sep = '\t', index_col = 0)
Sim_1_Ref = pd.read_csv(path + 'goldstandards/100_yeast3.tsv', sep = '\t', header = -1)
Sim_1_Ref.columns = ['node1', 'node2', 'value']
Sim_1_Ref = Sim_1_Ref.loc[Sim_1_Ref['value'] > 0]

Sim_2 = pd.read_csv(path + 'simulated_datasets/100_yeast3_medium.txt', sep = '\t', index_col = 0)
Sim_2_Ref = pd.read_csv(path + 'goldstandards/100_yeast3.tsv', sep = '\t', header = -1)
Sim_2_Ref.columns = ['node1', 'node2', 'value']
Sim_2_Ref = Sim_2_Ref.loc[Sim_2_Ref['value'] > 0]

HSC_Data = pd.read_excel(path + 'HSC_DATA/Moignard_HSC_Data.xlsx', index_col = 0, sheet_name = 0).T
HSC_Data[HSC_Data == 25] = 0
HSC_Data = HSC_Data[HSC_Data.sum(1) > 0]
HSC_Ref = pd.read_csv(path + 'HSC_DATA/Moignard_HSC_ReferenceEdges.tsv', sep = '\t')

ESC_Data = pd.read_excel(path + 'ESC_DATA/GSE59892_ESC_Dataset.xlsx', index_col = 0, sheet_name=0)
ESC_Gene_Symbol = pd.read_csv(path + 'ESC_DATA/GSE59892_ESC_Genesymobol.txt', sep = '\t')
ESC_Gene_Symbol_dict = dict(zip(ESC_Gene_Symbol.Ref_ID, ESC_Gene_Symbol.GeneSymbol))
ESC_Ref = pd.read_csv(path + 'ESC_DATA/GSE59892_ESC_ReferenceEdges.tsv', sep = '\t')

ESC_Data.index = list(ESC_Gene_Symbol_dict.get(i) for i in list(ESC_Data.index))
ESC_Data = ESC_Data[ESC_Data.sum(1) > 0]
# =============================================================================
# Test with Preprocessing and docker run
# =============================================================================



Sim1_gne = gnetdata.Gnetdata(Sim_1)
Sim2_gne = gnetdata.Gnetdata(Sim_2)
ESC_gne = gnetdata.Gnetdata(ESC_Data)
HSC_gne = gnetdata.Gnetdata(HSC_Data)

gne_Sim1_PIDC = gdocker.rundocker(Sim1_gne.deepcopy, method = 'PIDC')
gne_Sim1_GENIE3 = gdocker.rundocker(Sim1_gne.deepcopy, method = 'GENIE3')
gne_Sim1_CORR = gdocker.rundocker(Sim1_gne.deepcopy, method = 'CORR')
time_point = pd.DataFrame({'cell_id': list(Sim_1.columns), 'time_point': [1 for i in range(Sim1_gne.shape[1])]})
gne_Sim1_SCODE = gdocker.rundocker(Sim1_gne.deepcopy, method = 'SCODE', time_point = time_point)

gne_Sim2_PIDC = gdocker.rundocker(Sim2_gne.deepcopy, method = 'PIDC')
gne_Sim2_GENIE3 = gdocker.rundocker(Sim2_gne.deepcopy, method = 'GENIE3')
gne_Sim2_CORR = gdocker.rundocker(Sim2_gne.deepcopy, method = 'CORR')
time_point = pd.DataFrame({'cell_id': list(Sim_2.columns), 'time_point': [1 for i in range(Sim2_gne.shape[1])]})
gne_Sim2_SCODE = gdocker.rundocker(Sim2_gne.deepcopy, method = 'SCODE', time_point = time_point)

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

gne_HSC_PIDC.save_as(path + 'gnetData/HSC_PIDC.pk')
gne_HSC_GENIE3.save_as(path + 'gnetData/HSC_GENIE3.pk')
gne_HSC_CORR.save_as(path + 'gnetData/HSC_CORR.pk')
gne_HSC_SCODE.save_as(path + 'gnetData/HSC_SCODE.pk')


#Add node centralities
sim1_list = [gne_Sim1_CORR, gne_Sim1_GENIE3, gne_Sim1_PIDC, gne_Sim1_SCODE]
sim2_list = [gne_Sim2_CORR, gne_Sim2_GENIE3, gne_Sim2_PIDC, gne_Sim2_SCODE]

for obj in sim1_list:
        obj = gdocker.buildnet(obj, top = 500)
        obj = gt.community_detect(obj)
        obj = gt.get_centrality(obj)


def build_table(gnedata):

        links = gnedata.NetAttrs['links']
        centralities = gnedata.NetAttrs['centralities']

        centralities['avg'] = centralities.mean(axis = 1)
        centrality_avg = list()
        for source, target in zip(links.source, links.target):

                if source in list(centralities.node):
                        source_val = centralities[centralities.node == source].avg.values[0]
                else:
                        source_val =0
                if target in list(centralities.node):
                        target_val = centralities[centralities.node == target].avg.values[0]
                else:
                        target_val = 0

                centrality_avg.append(source_val + target_val)

        links['centrality_avg'] = centrality_avg
        return(links)


gne_HSC_CORR = gnetdata.load_Gnetdata_object(path + 'gnetData/HSC_CORR.pk')
gne_HSC_PIDC = gnetdata.load_Gnetdata_object(path + 'gnetData/HSC_PIDC.pk')
gne_HSC_GENIE3 = gnetdata.load_Gnetdata_object(path + 'gnetData/HSC_GENIE3.pk')
gne_HSC_SCODE = gnetdata.load_Gnetdata_object(path + 'gnetData/HSC_SCODE.pk')

#Sim1_links_dict = {'genie3': build_table(gne_Sim1_GENIE3), 'corr': build_table(gne_Sim1_CORR),
#              'pidc': build_table(gne_Sim1_PIDC), 'scode': build_table(gne_Sim1_SCODE)}
#
#Sim2_links_dict = {'genie3': build_table(gne_Sim2_GENIE3), 'corr': build_table(gne_Sim2_CORR),
#              'pidc': build_table(gne_Sim2_PIDC), 'scode': build_table(gne_Sim2_SCODE)}
#
#HSC_links_dict = {'genie3': build_table(gne_HSC_GENIE3), 'corr': build_table(gne_HSC_CORR),
#              'pidc': build_table(gne_HSC_PIDC), 'scode': build_table(gne_HSC_SCODE)}
#
#ESC_links_dict = {'genie3': build_table(gne_ESC_GENIE3), 'corr': build_table(gne_ESC_CORR),
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


Eclass = gt.ensemble_classifier(ESC_links_dict, threshold=0.45, num_trees=3000, seed=7, model = 'RF')
Eclass_pos = Eclass.sort_values('weight', ascending = False).head(500)

RF_fpr, RF_tpr, RF_auc = Eva_Test.test_run(links = Eclass_pos,
                                           Ref_links=ESC_Ref, input_dataset = ESC_Data)

PIDC_fpr, PIDC_tpr, PIDC_auc = Eva_Test.test_run(links = gne_ESC_PIDC.NetAttrs['links'].sort_values('weight', ascending = False).head(500),
                                           Ref_links=ESC_Ref, input_dataset = ESC_Data)

GENIE3_fpr, GENIE3_tpr, GENIE3_auc = Eva_Test.test_run(links = gne_ESC_GENIE3.NetAttrs['links'].sort_values('weight', ascending = False).head(500),
                                                       Ref_links=ESC_Ref, input_dataset = ESC_Data)

CORR_fpr, CORR_tpr, CORR_auc = Eva_Test.test_run(links = gne_ESC_CORR.NetAttrs['links'].sort_values('weight', ascending = False).head(500),
                                           Ref_links=ESC_Ref, input_dataset = ESC_Data)

SCODE_fpr, SCODE_tpr, SCODE_auc = Eva_Test.test_run(links = gne_ESC_SCODE.NetAttrs['links'].sort_values('weight', ascending = False).head(500),
                                           Ref_links=ESC_Ref, input_dataset = ESC_Data)

fpr_dict = {'RF': RF_fpr, 'GENIE3': GENIE3_fpr, 'PIDC': PIDC_fpr, 'CORR': CORR_fpr, 'SCODE': SCODE_fpr}
tpr_dict = {'RF': RF_tpr, 'GENIE3': GENIE3_tpr, 'PIDC': PIDC_tpr, 'CORR': CORR_tpr, 'SCODE': SCODE_tpr}
auc_dict = {'RF': RF_auc, 'GENIE3': GENIE3_auc, 'PIDC': PIDC_auc, 'CORR': CORR_auc, 'SCODE': SCODE_auc}

import seaborn as sns
colors = sns.color_palette().as_hex()[1:6]
plt.figure()
for i, color in zip(['RF', 'GENIE3', 'PIDC', 'CORR', 'SCODE'], colors):
        plt.plot(fpr_dict[i], tpr_dict[i], color = color, label = 'ROC curve of method {0} (area = {1:0.2f})'.format(i, auc_dict[i]))

plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Algorithms performance')
plt.legend(loc="lower right")
#plt.show()

plt.savefig('/home/mwu/test.png')
plt.close()

obj = gne_HSC_PIDC.deepcopy
obj = gdocker.buildnet(obj, top = 200)
obj = gt.community_detect(obj)
obj = gt.get_centrality(obj)
path = gt.random_walk(obj, start='Gata1', supervisedby='pageRank', steps=10)

tmp1 = ['Gata1']
tmp1.extend(list(obj.NetAttrs['graph']['Gata1'].keys()))
sn.static_netShow(obj, '/home/mwu/test_1.pdf', figure_size=[10, 10],
                  start = 'Gata1', neighhours_highlight=True,
                  scale=2, node_size = 200, font_size = 10, edge_color = 'grey')
