#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mwu
"""

from __future__ import absolute_import
from ._nx2gt import nx2gt
import graph_tool.all as gt

import seaborn as sns
import numpy as np
from pyvis.network import Network


def dynamic_netShow(gnetdata, filename, node_size=50, node_community='all', html_size=["500px", "1800px"],
                    bgcolor="#222222", font_color="white", **kwarg):
    """
    create a GRN html
    -------------------------
    :param gnetdata: Gnetdata object
    :param filename: str, save as filename
    :param node_size: int, default 50.
    :param node_community: str or list, default all.
    :param html_size: list, default ["1200px", "1600px"]
    :param font_color: string, default white
    :param bgcolor: string, default #222222
    :return: None
    """
    assert 'graph' in gnetdata.NetAttrs.keys(), 'graph is empty!'
    assert 'communities' in gnetdata.NetAttrs.keys(), 'node communities is empty!'

    node_group = gnetdata.NetAttrs['communities']
    graph = gnetdata.NetAttrs['graph']

    net = Network(html_size[0], html_size[1], bgcolor=bgcolor, font_color=font_color, **kwarg)
    edge_data = zip(list(graph.edges), list(graph[x[0]][x[1]]['weight'] for x in graph.edges))
    colors = sns.color_palette().as_hex() + sns.color_palette('Paired', 100).as_hex()

    for e in edge_data:
        src = e[0][0]
        dst = e[0][1]
        w = e[1]

        if node_community == 'all':
            src_color = int(node_group[node_group['node'] == src]['group'])
            dst_color = int(node_group[node_group['node'] == dst]['group'])
        else:
            sub_node_group = node_group[node_group.group.isin(node_community)]
            src_color = int(sub_node_group[sub_node_group['node'] == src]['group']) if src in list(
                sub_node_group.node) else 'grey'
            dst_color = int(sub_node_group[sub_node_group['node'] == dst]['group']) if dst in list(
                sub_node_group.node) else 'grey'

        net.add_node(src, src, title=src, color='grey' if src_color == 'grey' else colors[src_color])
        net.add_node(dst, dst, title=dst, color='grey' if dst_color == 'grey' else colors[dst_color])
        net.add_edge(src, dst, value=w)

    neighbor_map = net.get_adj_list()

    # add neighbor data to node hover data
    for node in net.nodes:
        node["title"] += " Neighbors:<br>" + "<br>".join(neighbor_map[node["id"]])
        node["value"] = len(neighbor_map[node["id"]])
        node['size'] = node_size

    net.show_buttons(filter_='physics')
    net.show(filename)

    return None


def net_hierarchy_plot(gnetdata, filename=None, **kwarg):
    """
     create a hierarchy gene net plot
    -------------------------
    :param gnetdata: Gnetdata object
    :param filename: str, default None.
    :param kwarg: additional parameters passed to graph_tool.all.draw_hierarchy()
    :return: None
    """

    assert 'graph' in gnetdata.NetAttrs.keys(), 'graph is empty!'
    assert 'communities' in gnetdata.NetAttrs.keys(), 'node communities is empty!'

    graph = nx2gt(gnetdata.NetAttrs['graph'])
    node_group = gnetdata.NetAttrs['communities']

    # deg = graph.degree_property_map('total') if vertex_size is None else vertex_size
    ngroup = graph.new_vertex_property('int')

    labels = dict(zip(list(range(graph.num_vertices())), list(graph.vertex_properties['id'])))
    for g in labels.keys():
        ngroup[g] = node_group.loc[node_group.node == labels[g], 'group']

    state = gt.minimize_nested_blockmodel_dl(graph, deg_corr=True)
    gt.draw_hierarchy(state, vertex_fill_color=ngroup, vertex_anchor=0,
                      vertex_text=graph.vertex_properties['id'],
                      output=filename, **kwarg)
    return None


def net_matrix_plot(gnetdata, filename=None, highlight_path=None, **kwarg):
    """
    create a matrix gene net plot
    --------------------------------
    :param gnetdata: Gnetdata object
    :param filename: str, default None.
    :param highlight_path: list, default None. a list of gene Nodes
    :param kwarg: additional parameters passed to graph_tool.all.graph_draw()
    :return: None
    """

    assert 'graph' in gnetdata.NetAttrs.keys(), 'graph is empty!'
    assert 'communities' in gnetdata.NetAttrs.keys(), 'node communities is empty!'

    graph = nx2gt(gnetdata.NetAttrs['graph'])
    node_group = gnetdata.NetAttrs['communities']

    ngroup = graph.new_vertex_property('int')

    labels = dict(zip(list(range(graph.num_vertices())), list(graph.vertex_properties['id'])))
    for g in labels.keys():
        ngroup[g] = node_group.loc[node_group.node == labels[g], 'group']

    position = graph.new_vertex_property("vector<double>")
    edge_color = graph.new_edge_property("string")

    dim = int(np.sqrt(len(labels))) + 1
    node_pos = list()

    node_group_dict = dict(zip(list(range(graph.num_vertices())), list(node_group.group)))

    for i in range(len(set(node_group.group))):
        for key, item in node_group_dict.items():
            if item == i:
                node_pos.append(key)

    for i in range(len(node_pos)):
        position[graph.vertex(node_pos[i])] = (int(i / dim), i % dim)

    for e in graph.edges():
        if (highlight_path is not None) and (labels[list(e)[0]] in highlight_path) and (
                labels[list(e)[1]] in highlight_path) and abs(
                highlight_path.index(labels[list(e)[0]]) - highlight_path.index(labels[list(e)[1]])) == 1:
            #         print((labels[list(e)[0]], labels[list(e)[1]]))
            edge_color[e] = 'darkorange'
        else:
            edge_color[e] = 'lightgrey'

    gt.graph_draw(graph, pos=position, vertex_fill_color=ngroup, edge_color=edge_color,
                  vertex_text=graph.vertex_properties['id'], output=filename, **kwarg)

    return None
