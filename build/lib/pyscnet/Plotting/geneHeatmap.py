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
    """
    Create Heatmap showing gene expression in individual cells along pseudotime
    ----------------------------------------------------------------------------
    :param gnetdata: Gnetdata object, default None
    :param gene: list, default None.
    :param cell_clusterid: str, default None. cell with cell_clusterid will be selected
    :param select_by: str, default None. key of filtering cells
    :param order_by: str, default None. key of ordering cells
    :param scale_data: bool, default True. whether or not scale the data
    :param cmap: str, default 'RdBu'. string denoting colors in clustermap
    :param save_as: str, default None. filepath+filename
    :param kwargs: additional parameters passed to seaborn.clustermap()
    :return: None
    """
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
    return None


def geneDynamic(gnetdata, gene, cell_clusterid, select_by, rolling=10, order_by=None, scale_data=True, save_as=None,
                colors=None, legend_size=None, **kwargs):
    """
    Create line plot showing gene dynamics
    ----------------------------------------------------------
    :param gnetdata: Gnetdata object
    :param gene: list, default None.
    :param cell_clusterid: str, default None. cell with cell_clusterid will be selected
    :param select_by: str, default None. key of filtering cells
    :param order_by: str, default None. key of ordering cells
    :param rolling: int, default 10. rolling window calculation
    :param scale_data: bool, default True. whether or not scale the data
    :param save_as: str, default None. filepath+filename
    :param colors: list, default None. list of string denoting color map
    :param legend_size: list, default None. specify legend size
    :param kwargs: additional parameters passed to matplotlib.pyplot.plot()
    :return: none

    """
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
    return None


def geneCorrelation(gnetdata, gene, cell_clusterid, select_by, order_by=None, scale_data=True,
                    save_as=None, figsize=None, **kwargs):

    """
    Create gene correlation heatmap
    --------------------------------------------------
    :param gnetdata: Gnetdata object
    :param gene: list, default None.
    :param cell_clusterid: str, default None. cell with cell_clusterid will be selected
    :param select_by: str, default None. key of filtering cells
    :param order_by: str, default None. key of ordering cells
    :param scale_data: bool, default True. whether or not scale the data
    :param save_as: str, default None. filepath+filename
    :param figsize: list, default None. a list of int defining figure size
    :param kwargs: additional parameters passed to seaborn.clustermap()
    :return: None
    """
    cell_info = gnetdata.CellAttrs['CellInfo'] if order_by is None else gnetdata.CellAttrs['CellInfo'].sort_values(
        order_by, ascending=True)
    cell = list(cell_info.loc[cell_info[select_by].isin([cell_clusterid])].index)
    sub_Expr = pd.DataFrame(preprocessing.scale(gnetdata.ExpMatrix),
                            index=gnetdata.ExpMatrix.index, columns=gnetdata.ExpMatrix.columns).loc[gene, cell].T \
        if scale_data else gnetdata.ExpMatrix.loc[gene, cell].T

    corr = sub_Expr.corr()
    mask = np.triu(np.ones_like(corr, dtype=np.bool))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    fig, ax = plt.subplots(figsize=[10, 10] if figsize is None else figsize)
    sns.heatmap(corr, mask=mask, cmap=cmap, vmax=.3, center=0, **kwargs)

    if save_as is not None:
        plt.savefig(save_as)

    plt.show()
    return None
