#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 10:04:57 2019

@author: mwu
"""
from __future__ import absolute_import
import sys
sys.path.append('/home/mwu/MING_V9T/PhD_Pro/SCNetEnrich')
#sys.path.append('/Users/angelawu/GitHub/SCNetEnrich')
from NetEnrich import measure_centrality as mc
from NetEnrich import random_walk as rw
from NetEnrich import top_gene_enrichment as genrich
from NetEnrich import de_bruijn_raw as debruijn
from NetEnrich import network_merge as nm
import networkx as nx

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import gseapy as gp


#link_table = pd.read_csv('/home/mwu/MING_V9T/PhD_Pro/SCNetEnrich/NetEnrich/links.txt', sep = '\t',
#                         header = -1)

link_table = pd.read_csv('/Users/angelawu/links_100_yeast3_large_dropout_low.txt', sep = '\t', header = -1)
link_table.columns = ['source', 'target', 'score']
link_table = link_table[link_table['score'] > 0.02]
link_table_2 = pd.read_csv('/Users/angelawu/links_100_yeast3_medium.txt', sep = '\t', header = -1)
link_table_2.columns = ['source', 'target', 'score']
link_table_2 = link_table_2[link_table_2['score'] > 0.02]

Net_test = mc.PyNet(link = link_table_2)
Net_test.get_centrality(source = Net_test.link.columns[0],
                        target = Net_test.link.columns[1])

Net_test.community_detect(filepath = '/Users/angelawu/links_100_yeast3_medium.pdf', font_size = 5)
path, links = rw.supervised_random_walk(Net_test, start='G12', steps=10, supervisedby='betweenness')

#pageRank_table = Net_test.pageRank
#random_gene_table = pd.DataFrame()
#random_gene_table['node'] = pd.Series((range(1,100))).apply(lambda x: 'G'+str(x))
#random_gene_table['pageRank'] = pd.Series(np.random.normal(0,1,size = (100)))
#
#pageRank_table['node']
#random_gene_table['node']
#tmp = genrich.gene_compare(list(pageRank_table['node']),list(random_gene_table['node']))


path = ['G12', 'G9', 'G29', 'G18', 'G16']
path1 = ['G12', 'G9', 'G8', 'G6', 'G50', 'G39']
tmp_dict = {}
import string
i = 0
for item in list(set(path + path1)):
    tmp_dict[item] = string.ascii_lowercase[i]
    i += 1

tmp = ""
tmp+=tmp.join(tmp_dict[x] for x in path)

tmp1 = ""
tmp1+=tmp1.join(tmp_dict[x] for x in path1)

g = debruijn.construct_graph([tmp1, tmp], 2)
print(debruijn.output_contigs(g))

node_list = set(list(link_table['source']) + list(link_table['target']) + list(link_table_2['source']) + list(link_table_2['target']))
#test with network merge functions
MergeNet_test = nm.MergedNetwork(link_1 = link_table, link_2 = link_table_2, node_list = node_list)

MergeNet_test.network_merge()
#MergeNet_test.knn_based_merge()
Net_test.Graph = MergeNet_test.Graph
Net_test.community_detect(filepath = '/Users/angelawu/links_100_yeast3_merged.pdf', font_size = 5)

















