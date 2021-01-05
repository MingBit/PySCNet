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
from ._random_walk import greedy_walk, supervised_random_walk
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
                                edge_attr='weight')
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


def find_consensus_graph(gnetdata, link_key='all', method='intersection',
                         top_rank=100, set_train=[0.05, 0.5], **kwargs):
    """
    Given multiple linkage tables, it predicts consensus links.
    ------------------------------------------------------------
    :param gnetdata: Gnetdata object.
    :param link_key: list, default all. key referring to linkage table.
    :param method: str, default intersection. methods for detecting consensus links. Note: intersection is recommended when there are less than 3 linkage tables.
    :param top_rank: int, default 100. top ranked edges for intersection method.
    :param set_train: list, default [0.05, 0.5]. Edges are ranked by weight obtained from individual methods. To train the classification model, we set top 5% edges as 'consensus edges (1)' and bottom 50% edges as 'non-consensus edges (0)' individual methods
    :return: Gnetdata object with consensus links added into NetAttrs.
    """
    assert method in ['intersection', 'snf',
                      'ensemble'], 'only following methods acceptable : intersection, snf, ensemble'
    keys = list(filter(lambda x: 'links' in x, gnetdata.NetAttrs.keys())) if link_key == 'all' else link_key

    if method == 'intersection':
        merged_links = gnetdata.NetAttrs[keys[0]].sort_values('weight', ascending=False, ignore_index=True).head(top_rank)
        for i in range(1, len(keys)):
            sorted_links = gnetdata.NetAttrs[keys[i]].sort_values('weight',
                                                                  ascending=False, ignore_index=True).head(top_rank)
            merged_links = graph_merge(merged_links, sorted_links, method=method)

    elif method == 'ensemble':
        links_dict = dict(filter(lambda i: i[0] in keys, gnetdata.NetAttrs.items()))
        X, Y, df_final = _generate_x_y(links_dict, top_rank=top_rank, set_train=set_train)
        merged_links = ensemble_classifier(X, Y, df_final, **kwargs)

    print('there are {} consensus edges found!'.format(merged_links.shape[0]))
    gnetdata._add_netattr('consensus_links', merged_links)

    return gnetdata


def graph_merge(link_1, link_2, method='union'):
    """
    Given two graphs, it returns merged graph.
    ----------------------------------------------
    :param link_1: dataframe. linkage table of graph_1
    :param link_2: dataframe. linkage table of graph_2
    :param method: str, default union. methods:[union, intersection, snf]. snf refers to similarity network fusion.
    :return: dataframe, merged linkage
    """
    assert method in ['union', 'intersection', 'snf'], 'valid method parameter: union, intersection, knn!'

    if method in ['union', 'intersection']:
        link_1['edge'] = ["_".join(sorted([link_1.source[i], link_1.target[i]])) for i in range(link_1.shape[0])]
        link_2['edge'] = ["_".join(sorted([link_2.source[i], link_2.target[i]])) for i in range(link_2.shape[0])]

        mergedlinks = pd.merge(link_1, link_2, on=['edge'], how='outer' if method == 'union' else 'inner').fillna(0)
        mergedlinks['weight'] = mergedlinks.mean(axis=1)

        mergedlinks['source'] = [x.split('_')[0] for x in mergedlinks.edge]
        mergedlinks['target'] = [x.split('_')[1] for x in mergedlinks.edge]

    elif method == 'snf':
        mergedlinks = __snf_based_merge(link_1, link_2)

    return mergedlinks[['source', 'target', 'weight']]


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


def self_guide_walk(gnetdata, start, method='greedy_walk', supervisedby='pageRank', steps=10, repeat=100):
    """
    For given a starting node, it provides greedy_walk and supervised random walk guided by node centrality.
    ------------------------------------------------------------
    :param gnetdata: Gnetdata object
    :param start: str, starting point of graph.
    :param method: str, 'greedy_walk' or 'supervised_random_walk'. default: greedy walk -- it prefers neighbors with maximum product of node centrality and weight. supervided_random_walk -- it prefers neighbors with randomized probability weighted by the product of node centrality and weight.
    :param supervisedby: str, 'betweenness', 'closeness', 'degree' or 'pageRank'. default: pageRank
    :param steps: int, number of steps. default: 10
    :param repeat: int, repeat time for supervised random walk. default: 100
    :return: a list of travelled nodes.
    """
    assert method in ['greedy_walk',
                      'supervised_random_walk'], 'method must be either greedy_walk or supervised_random_walk'
    if method == 'greedy_walk':
        path = greedy_walk(gnetdata=gnetdata, start=start, supervisedby=supervisedby, steps=steps)
    else:
        path = supervised_random_walk(gnetdata=gnetdata, start=start, supervisedby=supervisedby, steps=steps,
                                      repeat=repeat)
    return path


