#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mwu
"""

from __future__ import absolute_import
from .__nx2gt import nx2gt
import graph_tool.all as gt


def net_hierarchy_plot(gnetdata, vertex_size=None, filename=None, **kwarg):
    assert 'graph' in gnetdata.NetAttrs.keys(), 'graph is empty!'
    assert 'communities' in gnetdata.NetAttrs.keys(), 'node communities is empty!'

    graph = nx2gt(gnetdata.NetAttrs['graph'])
    node_group = gnetdata.NetAttrs['communities']

    deg = graph.degree_property_map('total') if vertex_size is None else vertex_size
    ngroup = graph.new_vertex_property('int')

    labels = dict(zip(list(range(graph.num_vertices())), list(graph.vertex_properties['id'])))
    for g in labels.keys():
        ngroup[g] = node_group.loc[node_group.node == labels[g], 'group']

    state = gt.minimize_nested_blockmodel_dl(graph, deg_corr=True)
    gt.draw_hierarchy(state, vertex_fill_color=ngroup, vertex_size=deg, vertex_anchor=0,
                      vertex_text=graph.vertex_properties['id'], vertex_font_size=vertex_size,
                      output=filename)
    return None



