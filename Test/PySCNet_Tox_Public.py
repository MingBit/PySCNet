#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 14:02:55 2019

@author: mwu
"""
from __future__ import absolute_import
import scipy.io as si


import pandas as pd
import numpy as np
import sys
sys.path.append('/home/mwu/MING_V9T/PhD_Pro/pyscnet/')
import _pickle as pk
from pyscnet.Preprocessing import gnetdata
from pyscnet.Preprocessing import general_pipeline as pipeline
from pyscnet.BuildNet import gne_dockercaller as gdocker
from pyscnet.NetEnrich import graph_toolkit as gt
from pyscnet.Plotting import show_net as sn
import matplotlib.pyplot as plt

path = '/home/mwu/MING_V9T/PhD_Pro/Test/Tox_Public_scData/'

KO_Expr = pd.DataFrame(si.mmread(path + 'GSM3568594_scRNA_D7_BMC_TOXKO_GP33_matrix.mtx').todense())
KO_GeneInfo = pd.read_csv(path + 'GSM3568594_scRNA_D7_BMC_TOXKO_GP33_genes.tsv', delimiter = '\t',
                          header = None, index_col = 1)
KO_CellInfo = pd.read_csv(path + 'GSM3568594_scRNA_D7_BMC_TOXKO_GP33_barcodes.tsv', delimiter = '\t',
                       header = None, index_col = 0)
KO_Expr.index = KO_GeneInfo.index
KO_Expr.columns = [('KO_' + i) for i in KO_CellInfo.index]

WT_Expr = pd.DataFrame(si.mmread(path + 'GSM3568593_scRNA_D7_BMC_WT_GP33_matrix.mtx').todense())
WT_GeneInfo = pd.read_csv(path + 'GSM3568593_scRNA_D7_BMC_WT_GP33_genes.tsv', delimiter = '\t',
                       header = None, index_col = 1)
WT_CellInfo = pd.read_csv(path + 'GSM3568593_scRNA_D7_BMC_WT_GP33_barcodes.tsv', delimiter = '\t',
                       header = None, index_col = 0)
#WT_CellInfo['Condition'] = 'WT'
WT_Expr.index = [x.upper() for x in WT_GeneInfo.index]

WT_Expr.columns = [('WT_' + i) for i in WT_CellInfo.index]
WT_Expr = WT_Expr.loc[~WT_Expr.index.duplicated(keep='first')]

Run29_D8 = pd.read_excel('/home/mwu/MING_V9T/FA_DEG_Heatmap/output/DEG/Filter_DEG_Run29_D8_WT_KO.xlsx')
Run32_D20 = pd.read_excel('/home/mwu/MING_V9T/FA_DEG_Heatmap/output/DEG/Filter_DEG_Run32_Day20_WT_KO.xlsx')
Run32_D8_Tim3p = pd.read_excel('/home/mwu/MING_V9T/FA_DEG_Heatmap/output/DEG/Filter_DEG_Run32_Day8_TIM3p_WT_KO.xlsx')
Run32_D8_Tim3n = pd.read_excel('/home/mwu/MING_V9T/FA_DEG_Heatmap/output/DEG/Filter_DEG_Run32_Day8_TIM3n_WT_KO.xlsx')

Tox_Expr = KO_Expr.join(WT_Expr)
tmp = ['KO' for i in range(len(KO_CellInfo.index))]
tmp.extend(['WT' for i in range(len(WT_CellInfo.index))])
Tox_CellInfo = pd.DataFrame({'cellName': Tox_Expr.columns,
                             'Condition': tmp})

Mms_fator = pd.read_csv('/home/mwu/MING_V9T/PhD_Pro/pyscnet/Mus_TF_and_TFcofactor/Mus_musculus_TF.txt', sep = '\t')
Mms_fator['Symbol'] = [x.upper() for x in list(Mms_fator['Symbol'])]


Tox_Expr.index = [x.upper() for x in list(Tox_Expr.index)]

comm = list((set([x.upper() for x in list(Run29_D8.geneName)]) |
            set([x.upper() for x in list(Run32_D20.geneName)]) |
            set([x.upper() for x in list(Run32_D8_Tim3p.geneName)]) |
            set([x.upper() for x in list(Run32_D8_Tim3n.geneName)])) & set(Tox_Expr.index))

comm.append('TOX')
Expr = WT_Expr.loc[comm,]

Tox_GeneInfo = pd.DataFrame({'Gene': Expr.index,
                          'TF_Gene': ['TF' if x in list(Mms_fator['Symbol']) else 'Gene' for x in Expr.index]})


Tox_gne = gnetdata.Gnetdata(ExpMatrix = WT_Expr)
Tox_gne.CellAttrs['Design'] = WT_CellInfo
Tox_gne = pipeline.run_pipeline(Tox_gne, pipeline = 'scanpy')

Tox_gne.GeneAttrs = Tox_GeneInfo
top_Genes = list(Expr.mean(1).sort_values(ascending = False).head(300).index)
Tox_gne.ExpMatrix = Expr.loc[top_Genes,]

Tox_gne_GENIE3 = gdocker.rundocker(Tox_gne.deepcopy, method = 'GENIE3')
Tox_Links = Tox_gne_GENIE3.NetAttrs['links']
Tox_Links['cell_clusterid'] = ['All' for i in range(Tox_Links.shape[0])]

cell_info = Tox_gne_GENIE3.CellAttrs['Design']
top_Genes = list(Expr.mean(1).sort_values(ascending = False).head(100).index)

for i in range(5):
        print(i)
        Tmp_Expr = Expr[cell_info[cell_info.louvain == str(i)].index]
        top_Genes = list(Tmp_Expr.mean(1).sort_values(ascending = False).head(200).index)
        Tmp_Expr = Tmp_Expr.loc[top_Genes,]
        Tmp_gne = gnetdata.Gnetdata(ExpMatrix = Tmp_Expr)
        Tmp_gne_GENIE3 = gdocker.rundocker(Tmp_gne.deepcopy, method = 'GENIE3')
        Tmp_Links = Tmp_gne_GENIE3.NetAttrs['links']
        Tmp_Links['cell_clusterid'] = [str(i+1) for j in range(Tmp_Links.shape[0])]
        Tox_Links = Tox_Links.append(Tmp_Links)

Tox_gne_GENIE3.NetAttrs['links'] = Tox_Links
Tox_gne_GENIE3.CellAttrs['Design'].columns = ['n_genes', 'n_counts', 'cluster_id']
Tox_gne_GENIE3.CellAttrs['Design']['cell_type'] = ['WT' for i in range(4451)]
Tox_gne_GENIE3.save_as('/home/mwu/MING_V9T/PhD_Pro/Test/Tox_Public_scData/WT_Tox_Selected_MarkerGenes.pk')


tmp2 = gnetdata.load_Gnetdata_object('/home/mwu/MING_V9T/PhD_Pro/Test/Tox_Public_scData/Tox_Selected_MarkerGenes.pk')










