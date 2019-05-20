#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 14:35:57 2019

@author: angelawu
"""



from __future__ import absolute_import
from BuildNet import gne_wgcna as wgcna
from BuildNet import gne_genie3 as genie3
from BuildNet import gne_bnn as bnn
#import scanpy.api as sc
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

#sim1_wgcna = wgcna.PyWGCNA(data = Sim_1, params = None)
#sim1_wgcna.adjacency(Nettype = 'unsigned')
#sim1_wgcna.Tom_dist_similarity()
#sim1_wgcna.generate_geneTree()
#sim1_wgcna.cutreeDynamic(minModuleSize = 5, maxTreeHeight=5)
#sim1_wgcna.module_eigen_gene()
#sim1_wgcna.plotDendroHeatmap(filepath = '/Users/angelawu/GitHub/SCNetEnrich/test.png')
#sim1_wgcna.plot_eigengene_network(filepath = '/Users/angelawu/GitHub/SCNetEnrich/test.png',
#                                  height = 800, width = 1000)

#tmp = wgcna_obj.plot_eigengene_network()

# =============================================================================
# test with genie3 in python
# =============================================================================
sim2_genie3 = genie3.PyGenie3(data = Sim_1, params = None)
sim2_genie3.run_genie3(nthreads = 4, filename = 'links_100_yeast3_medium.txt')

# =============================================================================
# test with bnn
# =============================================================================
params = {'output_path':'/Users/angelawu/GitHub/', 'filename':'sim2'}
sim2_bnn = bnn.Pybnn(data = Sim_2,
                     params = params)
sim2_bnn.run_bnn()


# =============================================================================
# test with bnlearn
# =============================================================================







