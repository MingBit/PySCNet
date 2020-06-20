#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:32:32 2019
@author: mwu
"""

import seaborn as sns
from pyvis.network import Network
import networkx as nx
import matplotlib.pyplot as plt
import copy
import matplotlib as mpl


def dynamic_netShow(gnetdata, filename, **kwargs):
    """ it returns a html file representing interactive network
        """
    if (gnetdata.NetAttrs['parameters']['threshold'] != 'None'):
        filterby = float(gnetdata.NetAttrs['parameters']['threshold'])
        link_table = gnetdata.NetAttrs['links'].loc[gnetdata.NetAttrs['links'].weight > filterby]
    elif (gnetdata.NetAttrs['parameters']['top'] != 'None'):
        filterby = int(gnetdata.NetAttrs['parameters']['top'])
        link_table = gnetdata.NetAttrs['links'].sort_values('weight', ascending=False).head(filterby)
    else:
        link_table = gnetdata.NetAttrs['links']
    node_group = gnetdata.NetAttrs['communities']

    net = Network("1200px", '1600px', bgcolor="#222222", font_color="white")
    # net = Network('1200px', '1600px')
    #        net.barnes_hut(**kwargs)
    edge_data = zip(link_table['source'], link_table['target'], link_table['weight'])

    colors = sns.color_palette().as_hex() + sns.color_palette('Paired', 100).as_hex()

    for e in edge_data:
        src = e[0]
        dst = e[1]
        w = e[2]

        if node_group is not None:
            src_color = int(node_group[node_group['node'] == src]['group'])
            dst_color = int(node_group[node_group['node'] == dst]['group'])
        else:
            src_color = 2
            dst_color = 2

        net.add_node(src, src, title=src, color=colors[src_color])
        net.add_node(dst, dst, title=dst, color=colors[dst_color])
        net.add_edge(src, dst, value=w)

    neighbor_map = net.get_adj_list()

    # add neighbor data to node hover data
    for node in net.nodes:
        node["title"] += " Neighbors:<br>" + "<br>".join(neighbor_map[node["id"]])
        node["value"] = len(neighbor_map[node["id"]])
        node['size'] = 50

    net.show_buttons(filter_='physics')
    net.show(filename)


def static_netShow(gnetdata, filename, scale=4, figure_size=[20, 10],
                   random_path=None, path_highlight=False, neighhours_highlight=False,
                   start=None, **kwargs):
    if gnetdata.NetAttrs['parameters']['threshold'] != 'None':
        filterby = float(gnetdata.NetAttrs['parameters']['threshold'])
        link_table = gnetdata.NetAttrs['links'].loc[gnetdata.NetAttrs['links'].weight > filterby]
    elif gnetdata.NetAttrs['parameters']['top'] != 'None':
        filterby = int(gnetdata.NetAttrs['parameters']['top'])
        link_table = gnetdata.NetAttrs['links'].sort_values('weight', ascending=False).head(filterby)
    else:
        link_table = gnetdata.NetAttrs['links']

    node_group = copy.deepcopy(gnetdata.NetAttrs['communities'])
    colors = sns.color_palette().as_hex() + sns.color_palette('Paired', 100).as_hex()
    node_color = [colors[n] for n in list(node_group.group)]

    net = nx.from_pandas_edgelist(link_table)
    pos = nx.kamada_kawai_layout(net, scale=scale)

    if path_highlight:
        if random_path is None:
            raise Exception('random_path must not be Null if path_highlight is True!')
        else:
            nx.draw_networkx(net, pos, node_color='grey', **kwargs)
            nx.draw_networkx_nodes(net, pos, nodelist=random_path, node_color=colors[0])
            for i in range(len(random_path) - 1):
                nx.draw_networkx_edges(net, pos, {(random_path[i], random_path[i + 1]): str(i + 1)},
                                       edge_color=colors[1], width=3)
                nx.draw_networkx_edge_labels(net, pos, {(random_path[i], random_path[i + 1]): str(i + 1)},
                                             font_color=colors[2], font_size=4)

    elif neighhours_highlight:
        if start is None:
            raise Exception('start node must not be None!')
        else:
            neighbors = [start]
            neighbors.extend(list(net[start].keys()))
            nx.draw_networkx(net, pos, node_color='grey', **kwargs)
            nx.draw_networkx_nodes(net, pos, nodelist=neighbors, node_color=colors[0])
            for i in range(len(neighbors)):
                nx.draw_networkx_edges(net, pos, {(start, neighbors[i]): str(i)},
                                       edge_color=colors[3], width=4)
                nx.draw_networkx_edge_labels(net, pos, {(start, neighbors[i]): str(i)},
                                             font_color=colors[2], font_size=4)
    else:
        nx.draw_networkx(net, pos, **kwargs)

    mpl.rcParams['figure.figsize'] = figure_size
    plt.savefig(filename)
    plt.show()
