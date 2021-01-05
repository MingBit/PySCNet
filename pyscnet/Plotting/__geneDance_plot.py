#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:32:32 2019
@author: mwu
"""


# def geneDynamic(gnetdata, gene, cell_clusterid, select_by, rolling=10, order_by=None, scale_data=True, save_as=None,
#                 colors=None, legend_size=None, **kwargs):
#     """
#     Create line plot showing gene dynamics
#     ----------------------------------------------------------
#     :param gnetdata: Gnetdata object
#     :param gene: list, default None.
#     :param cell_clusterid: str, default None. cell with cell_clusterid will be selected
#     :param select_by: str, default None. key of filtering cells
#     :param order_by: str, default None. key of ordering cells
#     :param rolling: int, default 10. rolling window calculation
#     :param scale_data: bool, default True. whether or not scale the data
#     :param save_as: str, default None. filepath+filename
#     :param colors: list, default None. list of string denoting color map
#     :param legend_size: list, default None. specify legend size
#     :param kwargs: additional parameters passed to matplotlib.pyplot.plot()
#     :return: none
#
#     """
#     cell_info = gnetdata.CellAttrs['CellInfo'] if order_by is None else gnetdata.CellAttrs['CellInfo'].sort_values(
#         order_by, ascending=True)
#     cell = list(cell_info.loc[cell_info[select_by].isin([cell_clusterid])].index)
#     sub_Expr = pd.DataFrame(preprocessing.scale(gnetdata.ExpMatrix),
#                             index=gnetdata.ExpMatrix.index, columns=gnetdata.ExpMatrix.columns).loc[gene, cell].T \
#         if scale_data else gnetdata.ExpMatrix.loc[gene, cell].T
#
#     colors = list(mcolors._colors_full_map.values()) if colors is None else colors
#     fig, ax = plt.subplots()
#     num = 0
#
#     if save_as is not None:
#         plt.savefig(save_as)
#
#     # plt.show()
#     return ax


