#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:32:32 2019
@author: mwu
"""

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from sklearn import preprocessing
import pandas as pd
import numpy as np


def geneHeatmap(gnetdata, gene, cell_clusterid, select_by, order_by=None, scale_data=True, cmap='RdBu', save_as=None,
                **kwargs):
    cell_info = gnetdata.CellAttrs['CellInfo'] if order_by is None else gnetdata.CellAttrs['CellInfo'].sort_values(
        order_by, ascending=True)
    cell = list(cell_info.loc[cell_info[select_by].isin([cell_clusterid])].index)
    sub_Expr = pd.DataFrame(preprocessing.scale(gnetdata.ExpMatrix),
                            index=gnetdata.ExpMatrix.index, columns=gnetdata.ExpMatrix.columns).loc[gene, cell] \
        if scale_data else gnetdata.ExpMatrix.loc[gene, cell]
    sns_plot = sns.clustermap(sub_Expr, cmap=cmap, xticklabels=False, **kwargs)

    if save_as is not None:
        sns_plot.savefig(save_as)

    plt.show()


def geneDynamic(gnetdata, gene, cell_clusterid, select_by, rolling=10, order_by=None, scale_data=True, save_as=None,
                colors=None, legend_size=None, **kwargs):
    cell_info = gnetdata.CellAttrs['CellInfo'] if order_by is None else gnetdata.CellAttrs['CellInfo'].sort_values(
        order_by, ascending=True)
    cell = list(cell_info.loc[cell_info[select_by].isin([cell_clusterid])].index)
    sub_Expr = pd.DataFrame(preprocessing.scale(gnetdata.ExpMatrix),
                            index=gnetdata.ExpMatrix.index, columns=gnetdata.ExpMatrix.columns).loc[gene, cell].T \
        if scale_data else gnetdata.ExpMatrix.loc[gene, cell].T

    colors = list(mcolors._colors_full_map.values()) if colors is None else colors
    fig, ax = plt.subplots()
    num = 0
    for gene in sub_Expr.columns:
        num = num + 1
        #         ax.plot(sub_Expr[gene], marker='o', color='0.6', linestyle=None)
        ax.plot(sub_Expr[gene].rolling(rolling, center=True).mean(), linewidth=10, color=colors[num],
                label=gene + '_' + str(rolling) + '-rolling mean', **kwargs)

        #         ax.set_xlabel('Pseudotime', fontsize=20)
        #         ax.set_ylabel('Scaled Expression', fontsize=20)
        ax.axis('off')
        ax.legend(prop={"size": 10 if legend_size is None else legend_size})

    if save_as is not None:
        plt.savefig(save_as)

    plt.show()


def geneCorrelation(gnetdata, gene, cell_clusterid, select_by, figsize=[10, 10], order_by=None, scale_data=True,
                    save_as=None, **kwargs):

    cell_info = gnetdata.CellAttrs['CellInfo'] if order_by is None else gnetdata.CellAttrs['CellInfo'].sort_values(
        order_by, ascending=True)
    cell = list(cell_info.loc[cell_info[select_by].isin([cell_clusterid])].index)
    sub_Expr = pd.DataFrame(preprocessing.scale(gnetdata.ExpMatrix),
                            index=gnetdata.ExpMatrix.index, columns=gnetdata.ExpMatrix.columns).loc[gene, cell].T \
        if scale_data else gnetdata.ExpMatrix.loc[gene, cell].T

    corr = sub_Expr.corr()
    mask = np.triu(np.ones_like(corr, dtype=np.bool))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(corr, mask=mask, cmap=cmap, vmax=.3, center=0,**kwargs)

    if save_as is not None:
        plt.savefig(save_as)

    plt.show()
