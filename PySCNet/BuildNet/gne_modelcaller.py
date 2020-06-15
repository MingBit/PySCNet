#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 12:09:41 2019

@author: mwu
"""

# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:51:43 2019

@author: mwu
"""

from PySCNet.BuildNet.Models import model_node2vec as nv


def call_node2vec(gnetdata, reference_links, p, q, dim_list, walk_list, num_walks_list,
                  workers, n_pc, param_grid, use_ref=True, cell_clusterid=None, **kwargs):
    if cell_clusterid is not None:

        cell_info = gnetdata.CellAttrs
        Expr = gnetdata.ExpMatrix[list(cell_info.loc[cell_info.cluster_id.isin([cell_clusterid])].index)]
    else:
        Expr = gnetdata.ExpMatrix

    for dimensions in dim_list:
        for walk_length in walk_list:
            for num_walks in num_walks_list:
                node_matrix = nv.run_node2vec(Expr, p=p, q=q, size=dimensions, walk_len=walk_length,
                                              num_walks=num_walks, workers=workers, n_comp=n_pc, **kwargs)
                if use_ref:
                    links = nv._binary_classifier(embedded_node=node_matrix, reference_links=reference_links,
                                                  select_n=0,
                                                  use_ref=True, param_grid=param_grid)
                else:
                    links = nv._binary_classifier(embedded_node=node_matrix, reference_links=reference_links,
                                                  select_n=int(reference_links.shape[0] * 0.25),
                                                  param_grid=param_grid)
    gnetdata._add_netattr('links', links)
    return (gnetdata)
