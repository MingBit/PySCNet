#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 20:50:39 2019

@author: mwu
"""
from __future__ import absolute_import
import networkx as nx
import pandas as pd
from community import community_louvain
import snf
import numpy as np
import warnings
import networkx.algorithms.traversal as nextra
from ._de_bruijn import construct_graph, output_contigs
from ._random_walk import supervised_random_walk
from ._ensemble_classifier import _generate_x_y, ensemble_classifier


def __init__():
    warnings.simplefilter("ignore")


def __linkage_to_adjlink(linkage_table, node_list):
    """
    convert linkage table to weighted adjacency matrix
    """
    adjlink_matrix = pd.DataFrame(0, columns=node_list, index=node_list, dtype=np.float)
    #        source, target, score = list(linkage_table.columns)
    for i in range(0, len(linkage_table)):
        if (linkage_table['source'][i] in node_list) & (linkage_table['target'][i] in node_list):
            adjlink_matrix.loc[linkage_table['source'][i]][linkage_table['target'][i]] = linkage_table['weight'][i]
            adjlink_matrix.loc[linkage_table['target'][i]][linkage_table['source'][i]] = linkage_table['weight'][i]
        else:
            break
    return np.array(adjlink_matrix)


def __snf_based_merge(link_1, link_2):
    """
    snf network merge based on Wang, Bo, et al. Nature methods 11.3 (2014): 333.
    """
    warnings.simplefilter("ignore")
    node_list = list(set(link_1['source']) & set(link_1['target']) & set(link_2['source']) & set(link_2['target']))

    adjlinks = list()
    adjlinks.append(__linkage_to_adjlink(link_1, node_list))
    adjlinks.append(__linkage_to_adjlink(link_2, node_list))
    affinity_matrix = snf.make_affinity(adjlinks)
    fused_network = snf.snf(affinity_matrix)
    Graph = nx.from_pandas_adjacency(pd.DataFrame(fused_network, index=node_list, columns=node_list))
    return pd.DataFrame(Graph.edges, columns=['source', 'target'])


def buildnet(gnetdata, key_links, top=None):
    """
    Given linkage table, build gene correlation graph
    ----------------------------------------------------
    :param gnetdata: Gnetdata object.
    :param key_links: str, key of links referring which linkage table for buidling graph
    :param top: int, default None. top ranked links
    :return: Gnetdata object with graph added into NetAttrs
    """
    top = gnetdata.NetAttrs[key_links].shape[0] if top is None else top
    links_filter = gnetdata.NetAttrs[key_links].sort_values('weight', ascending=False).head(top)

    G = nx.from_pandas_edgelist(links_filter,
                                source="source",
                                target="target",
                                edge_attr=True)
    gnetdata._add_netattr('graph', G)
    gnetdata._add_netattr_para('top', str(top))

    print('graph added into NetAttrs')
    return gnetdata


def get_centrality(gnetdata):
    """
    Measure node centrality in the network.
    ---------------------------------------
    :param gnetdata: Gnetdata object.
    :return: gnetData object with 'centralities' added into NetAttrs
    """
    G = gnetdata.NetAttrs['graph']
    centralities = pd.DataFrame(list(G.nodes), columns=['node'])
    centralities['betweenness'] = pd.DataFrame.from_dict(list(nx.betweenness_centrality(G).items()))[1]
    centralities['closeness'] = pd.DataFrame.from_dict(list(nx.closeness_centrality(G).items()))[1]
    centralities['degree'] = pd.DataFrame.from_dict(list(nx.degree_centrality(G).items()))[1]
    centralities['pageRank'] = pd.DataFrame.from_dict(list(nx.pagerank(G).items()))[1]

    gnetdata.NetAttrs['centralities'] = centralities
    print('node centralities added into NetAttrs.')
    return gnetdata


def detect_community(gnetdata, **kwargs):
    """
    Detect gene modules via louvain community dection algorithm.
    -------------------------------------------------------------
    :param gnetdata: Gnetdata object.
    :param kwargs: additional parameters passed to community_louvain.best_partition()
    :return: Gnetdata object with 'communities' added into NetAttrs
    """
    G = gnetdata.NetAttrs['graph']
    partition = community_louvain.best_partition(G, **kwargs)
    communities = pd.DataFrame.from_dict(list(partition.items()))
    communities.columns = ['node', 'group']

    gnetdata.NetAttrs['communities'] = communities
    print('gene communities added into NetAttrs')
    return gnetdata


def find_consensus_graph(gnetdata, link_key='all', method='intersection', toprank=100, threshold=None, **kwargs):
    """
    Given multiple linkage tables, it predicts consensus links.
    ------------------------------------------------------------
    :param gnetdata: Gnetdata object.
    :param link_key: str, default all. key referring to linkage table.
    :param method: str, default intersection. methods for detecting consensus links. Note: intersection is recommended when there are less than 3 linkage tables.
    :param toprank: int, default 100. top ranked edges for intersection method.
    :param threshold: int, default None. set threshold for ensemble method.
    :return: Gnetdata object with consensus links added into NetAttrs.
    """
    assert method in ['intersection', 'snf', 'ensemble'], 'only following methods acceptable : intersection, snf, ensemble'
    keys = list(filter(lambda x: 'links' in x, gnetdata.NetAttrs.keys())) if link_key == 'all' else link_key

    if method == 'intersection':
        merged_links = gnetdata.NetAttrs[keys[0]]
        for i in range(1, len(keys)):
            merged_links = graph_merge(merged_links, gnetdata.NetAttrs[keys[i]], method=method, toprank=toprank)

    elif method == 'ensemble':
        if threshold is None:
            raise Exception('threshold cannot be none!')
        links_dict = dict(filter(lambda i: i[0] in keys, gnetdata.NetAttrs.items()))
        X, Y = _generate_x_y(links_dict, threshold)
        merged_links = ensemble_classifier(X, Y, **kwargs)

    print('there are {} consensus edges found!'.format(merged_links.shape[0]))
    gnetdata._add_netattr('consensus', merged_links)

    return gnetdata


def graph_merge(link_1, link_2, toprank=None, method='union'):
    """
    Given two graphs, it returns merged graph.
    ----------------------------------------------
    :param link_1: dataframe. linkage table of graph_1
    :param link_2: dataframe. linkage table of graph_2
    :param toprank: int, default None. top edges from each methods for graph_merge
    :param method: str, default union. methods:[union, intersection, snf]. snf refers to similarity network fusion.
    :return: dataframe, merged linkage
    """
    assert method in ['union', 'intersection', 'snf'], 'valid method parameter: union, intersection, knn!'

    if method in ['union', 'intersection']:
        toprank = min(link_1.shape[0], link_2.shape[0]) if toprank is None else toprank
        link_1 = link_1.sort_values('weight', ascending=False).head(toprank)
        link_2 = link_2.sort_values('weight', ascending=False).head(toprank)
        mergedlinks = pd.merge(link_1.iloc[:, :-1], link_2.iloc[:, :-1], how='outer' if method == 'union' else 'inner')
        mergedlinks['weight'] = np.repeat(1, mergedlinks.shape[0])

    elif method == 'snf':
        mergedlinks = __snf_based_merge(link_1, link_2)

    return mergedlinks


def graph_traveral(graph, start, threshold, method='bfs'):
    """
    Given a graph, it provides graph traversal techniques including breadth-first search (bsf) and depth-first search (dfs) to explore hidden gene/tf associations.
    -----------------------------------------------------------------------------
    :param graph: network graph object.
    :param start: str. starting point of graph.
    :type start: str
    :param threshold: int. the depth-limit
    :param method: str. bfs or dfs
    :return: explored graph.
    """
    assert method in ['bfs', 'dfs'], 'valid method parameters: bfs, dfs!'
    if method == 'bfs':
        res_path = nextra.bfs_tree(graph, start, threshold)

    elif method == 'dfs':
        res_path = nextra.dfs_tree(graph, start, threshold)

    return res_path


def random_walk(gnetdata, start, supervisedby, steps):
    """
    Supervised random walk guided by node centrality attribute.
    ------------------------------------------------------------
    :param gnetdata: Gnetdata object
    :param start: str, starting point of graph.
    :param supervisedby: str, 'betweenness', 'closeness', 'degree' or 'pageRank'
    :param steps: int, number of steps.
    :return: a list of travelled nodes.
    """
    path = supervised_random_walk(gnetdata=gnetdata, start=start, supervisedby=supervisedby, steps=steps)
    return path


def path_merge(path_1, path_2, k_mer=3, path='Eulerian'):
    """
    TODO: perform de bruijn graph mapping for ginve two path lists
    """
    g = construct_graph([path_1, path_2], k_mer)
    merged_path = output_contigs(g)

    return merged_path
