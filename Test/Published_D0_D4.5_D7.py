#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 09:21:08 2019

@author: mwu
"""


from __future__ import absolute_import
import sys
if('/home/mwu/MING_V9T/PhD_Pro/PySCNet/' not in sys.path): sys.path.append('/home/mwu/MING_V9T/PhD_Pro/PySCNet/')
import seaborn as sns
import networkx as nx
import scanpy.api as sc
import pandas as pd
import copy
import numpy as np
from PySCNet.Preprocessing import gnetdata
from PySCNet.BuildNet import gne_dockercaller as gdocker
from PySCNet.BuildNet import gne_modelcaller as gmodel
from PySCNet.NetEnrich import graph_toolkit as gt
from PySCNet.Plotting import show_net as sn
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

path = '/home/mwu/MING_V9T/PhD_Pro/SC_GRN_Pro/Merge_All_SCData/'

#read seurat object via loom
D0_D4_loom = ad.read_loom(path + 'D0_Naive_D4.5_C13/Loom_K=6_PC=30_K.para=500.rds')
Expr = pd.DataFrame(D0_D4_loom.X.T.todense(), columns = list(D0_D4_loom.obs.index), index = list(D0_D4_loom.var.index))
Expr.index = [x.upper() for x in list(Expr.index)]
cell_info = pd.DataFrame(D0_D4_loom.obs)
cell_info.columns = ['cluster_id'] + list(cell_info.columns[1:10])

# create signatrue list and tf list
Signatures = pd.read_csv(path + 'Signatures/FastProject_Input/Merge_Signatures.txt', sep = '\t', header = None)
Signatures.columns = ['Signatures_type', 'Symbol']
Signatures['Symbol'] = [x.upper() for x in list(Signatures['Symbol'])]

Mms_fator = pd.read_csv('/home/mwu/MING_V9T/PhD_Pro/PySCNet/Mus_TF_and_TFcofactor/Mus_musculus_TF.txt', sep = '\t')
Mms_fator['Symbol'] = [x.upper() for x in list(Mms_fator['Symbol'])]

#comm = list(Expr.mean(1).sort_values(ascending = False).head(300).index)

comm = list(set([x.upper() for x in list(Signatures.Symbol)]) & set(Expr.index))
Expr = Expr.loc[comm, :]

gene_info = pd.DataFrame({'Gene': Expr.index, 'TF_Gene': ['TF' if x in list(Mms_fator['Symbol']) else 'Gene' for x in Expr.index]})

#create gnedata
top_gene = list(Expr.mean(1).sort_values(ascending = False).head(500).index)
test_gne = gnetdata.Gnetdata(ExpMatrix = Expr.loc[top_gene,])
test_gne.CellAttrs = cell_info
test_gne.GeneAttrs = gene_info
test_gne_GENIE3 = gdocker.rundocker(test_gne.deepcopy, method = 'GENIE3')
test_Links = test_gne_GENIE3.NetAttrs['links']
test_Links['cell_clusterid'] = ['All' for i in range(test_Links.shape[0])]

for i in range(6):
        print(i+1)
        Tmp_Expr = Expr[list(cell_info.loc[cell_info.cluster_id.isin([i+1])].index)]
        top_Genes = list(Tmp_Expr.mean(1).sort_values(ascending = False).head(500).index)
        Tmp_Expr = Tmp_Expr.loc[top_Genes,]
        Tmp_gne = gnetdata.Gnetdata(ExpMatrix = Tmp_Expr)
        Tmp_gne_GENIE3 = gdocker.rundocker(Tmp_gne.deepcopy, method = 'GENIE3')
        Tmp_Links = Tmp_gne_GENIE3.NetAttrs['links']
        Tmp_Links['cell_clusterid'] = [str(i+1) for j in range(Tmp_Links.shape[0])]
        test_Links = test_Links.append(Tmp_Links)



test_gne_GENIE3.NetAttrs['links'] = copy.deepcopy(test_Links)
test_gne_GENIE3.save_as(path + 'D0_Naive_D4.5_C13/GRN/D0_Naive_D4.5_C13_Top500.pk')

test_gne_PIDC = gdocker.rundocker(test_gne.deepcopy, method = 'PIDC')
test_gne_GENIE3 = gdocker.rundocker(test_gne.deepcopy, method = 'GENIE3')
test_gne_CORR = gdocker.rundocker(test_gne.deepcopy, method = 'CORR')
time_point = pd.DataFrame({'cell_id': list(Expr.columns), 'time_point': [1 for i in range(Expr.shape[1])]})
test_gne_SCODE = gdocker.rundocker(test_gne.deepcopy, method = 'SCODE', time_point = time_point)

test_gne_PIDC = gdocker.buildnet(test_gne_PIDC, top = 200)
test_gne_GENIE3 = gdocker.buildnet(test_gne_GENIE3, top = 200)

test_gne_GENIE3 = gt.get_centrality(test_gne_GENIE3)
test_gne_GENIE3 = gt.community_detect(test_gne_GENIE3)

test_gne_PIDC = gt.get_centrality(test_gne_PIDC)
test_gne_PIDC = gt.community_detect(test_gne_PIDC)

#save as pickle object
test_gne_GENIE3.save_as(path + 'cluster_4_GENIE3.pk')
test_gne_PIDC.save_as(path + 'cluster_4_PIDC.pk')

#dn.dynamic_netShow(test_gne_PIDC, filepath = path + 'cluster_5_wo_CD4UniqueEff_pidc.html')
#dn.dynamic_netShow(test_gne_GENIE3, filepath = path + 'cluster_5_wo_CD4UniqueEff_genie3.html')


gt_rw_dt = pd.DataFrame(columns=['step_' + str(x) for x in range(10)])
nodes = list(test_gne_GENIE3.NetAttrs['graph'].node)

for i in range(len(nodes)):
        rw_path = gt.random_walk(gnetdata = test_gne_GENIE3, start = nodes[i], supervisedby = 'pageRank', steps = 9)
        if len(rw_path) < 10:

                for i in range(10 - len(rw_path)):
                        rw_path.append('NA')

        gt_rw_dt.loc[i] = rw_path

gt_rw_dt.to_csv(path + 'cluster_5_wo_CD4UniqueEff_RandomPath_GENIE3.csv', sep=' ')

#merge test
merge_graph = gt.graph_merge(test_gne_GENIE3.NetAttrs['links'], test_gne_PIDC.NetAttrs['links'], method = 'knn')
merge_graph.edges


tmp = gnetdata.load_Gnetdata_object(path + 'test.pk')
tmp.CellAttrs.columns = ['cluster_id'] + list(tmp.CellAttrs.columns[1:9])
tmp.save_as(path + 'test.pk')


#conbine networks from each clusters
pickle_data = gnetdata.load_Gnetdata_object('/home/mwu/MING_V9T/PhD_Pro/SC_GRN_Pro/Merge_All_SCData/D0_Naive_D4.5_C13/GRN/D0_Naive_D4.5_C13_Top500.pk')
links = copy.deepcopy(pickle_data.NetAttrs['links'])


sub_links = pd.DataFrame(columns = links.columns, index = None)

for i in range(6):
    tmp_links = links[links.cell_clusterid.isin([str(i + 1)])]
#    top_links = pd.DataFrame(tmp_links.sort_values('weight', ascending=False).head(300))
    top_links = tmp_links[tmp_links['weight'] > np.quantile(tmp_links['weight'], 0.95)]
    
    sub_links = pd.concat([sub_links, top_links])
    
   
sub_links = gdocker._remove_duplicate(sub_links)
sub_links['edges'] = sub_links.source + '_' + sub_links.target

sub_links['edges_type'] = [''.join(list(sub_links[sub_links.edges == x].cell_clusterid)) for x in sub_links.edges]
filter_sub_links = sub_links.iloc[[x for x in range(sub_links.shape[0]) if (len(sub_links.edges_type[x]) < 6)]]

#update links for pickle data
pickle_data.NetAttrs['links'] = filter_sub_links
gdocker.buildnet(pickle_data)
#pickle_data.NetAttrs['graph'] = nx.from_pandas_edgelist(filter_sub_links)
gt.get_centrality(pickle_data)
gt.community_detect(pickle_data)

colors = sns.color_palette().as_hex() + sns.color_palette('Paired', 100).as_hex()
unique_edge_list = filter_sub_links.drop_duplicates(subset= ['source', 'target'], keep = 'first')
edges_colors = dict(zip(set(unique_edge_list.edges_type), colors[:len(set(unique_edge_list.edges_type))]))

node_group = copy.deepcopy(pickle_data.NetAttrs['communities'])
node_color=[colors[n] for n in list(node_group.group)]

sn.static_netShow(pickle_data, filename = '/home/mwu/test.pdf', figure_size = [30, 30], scale = 6,
                  edge_color = [edges_colors['3'] if x == '3' else 'grey' for x in unique_edge_list.edges_type],
                  width = 4, font_size = 8, node_color = 'grey')

links = pd.read_csv('/home/mwu/links_all.csv', sep = ' ', header = 0, index_col = 0)
cell_info = pickle_data.CellAttrs
Expr = pickle_data.ExpMatrix
links.to_csv('/home/mwu/links_all.csv', sep = ' ')

for i in range(6):
        
        print(i+1)
        Tmp_Expr = Expr[list(cell_info.loc[cell_info.cluster_id.isin([i+1])].index)]
        top_Genes = list(Tmp_Expr.mean(1).sort_values(ascending = False).head(200).index)
        Tmp_Expr = Tmp_Expr.loc[top_Genes,]
        Tmp_gne = gnetdata.Gnetdata(ExpMatrix = Tmp_Expr)
        Tmp_gne_node2vec = gmodel.call_node2vec(Tmp_gne.deepcopy, 
                                                p =1, q = 1, dim_list = [20], walk_list = [50],
                                                num_walks_list = [1000], workers = 15)
        
        Tmp_Links = Tmp_gne_node2vec.NetAttrs['links']
        Tmp_Links['cell_clusterid'] = [str(i+1) for j in range(Tmp_Links.shape[0])]
        links = links.append(Tmp_Links)
    
