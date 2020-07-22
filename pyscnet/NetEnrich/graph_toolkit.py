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
import copy
import networkx.algorithms.traversal as nextra
from ._de_bruijn import construct_graph, output_contigs
from ._random_walk import supervised_random_walk
from ._ensemble_classifier import _generate_x_y, ensemble_classifier


def __init__():
    warnings.simplefilter("ignore")


def _linkage_to_adjlink(linkage_table, node_list):
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


def _snf_based_merge(link_1, link_2):
    """
    :param link_1:
    :param link_2:
    :return:
    """
    warnings.simplefilter("ignore")
    node_list = list(set(link_1['source']) & set(link_1['target']) & set(link_2['source']) & set(link_2['target']))

    adjlinks = list()
    adjlinks.append(_linkage_to_adjlink(link_1, node_list))
    adjlinks.append(_linkage_to_adjlink(link_2, node_list))
    affinity_matrix = snf.make_affinity(adjlinks)
    fused_network = snf.snf(affinity_matrix)
    Graph = nx.from_pandas_adjacency(pd.DataFrame(fused_network, index=node_list, columns=node_list))
    return pd.DataFrame(Graph.edges, columns=['source', 'target'])


def get_centrality(gnetdata):
    """
    Measure node centrality in the network.
    :param gnetdata: gnetData object
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
    :param gnetdata: gnetData object
    :param kwargs: parameters for community_louvain.best_partition
    :return: gnetData object with 'communities' added into NetAttrs
    """
    G = gnetdata.NetAttrs['graph']
    partition = community_louvain.best_partition(G, **kwargs)
    communities = pd.DataFrame.from_dict(list(partition.items()))
    communities.columns = ['node', 'group']

    gnetdata.NetAttrs['communities'] = communities
    print('gene communities added into NetAttrs')
    return gnetdata


def find_consensus_graph(gnetdata, link_key='all', method='intersection', threshold=None, **kwargs):
    """
    Given multiple linkage tables, it predicts consensus links.
    :param gnetdata: gnetData object
    :param link_key: key of linkage tables
    :param method: methods for detecting consensus links. three methods provided: intersection, snf, ensemble
    :param threshold: set threshold for ensemble classifier training
    :return: return gnetData object with consensus links updated.
    """
    keys = list(filter(lambda x: 'links' in x, gnetdata.NetAttrs.keys())) if link_key == 'all' else link_key

    if method in ['intersection', 'snf']:
        merged_links = gnetdata.NetAttrs[keys[0]]
        for i in range(1, len(keys)):
            merged_links = graph_merge(merged_links, gnetdata.NetAttrs[keys[i]], method=method)

    elif method == 'ensemble':
        if threshold is None:
            raise Exception('threshold cannot be none!')
        links_dict = dict(filter(lambda i:i[0] in keys, gnetdata.NetAttrs.items()))
        merged_links = ensemble_classifier(_generate_x_y(links_dict, threshold=threshold), **kwargs)

    else:
        raise Exception('only following methods acceptable : intersection, snf, ensemble')

    gnetdata._add_netattr('consensus', merged_links)
    return gnetdata


def graph_merge(link_1, link_2, method='union'):
    """
    Given two graphs, it returns merged graph
    :param link_1: linkage table of graph_1
    :param link_2: linkage table of graph_2
    :param method: method can be 'union', 'intersection' or 'snf'. 'snf' refers to similarity network fusion algorithm.
    :return: merged graph
    """
    if method == 'union':
        union_links = pd.merge(link_1, link_2, how='outer')
        mergedlinks = union_links.groupby(['source', 'target'], as_index=False).mean().reindex()
        # Graph = nx.from_pandas_edgelist(mergedlinks,
        #                                 source='source',
        #                                 target='target',
        #                                 edge_attr=True)

    elif method == 'intersection':
        link_1_cp = copy.deepcopy(link_1)
        link_1_cp.columns = ['target', 'source', 'weight']
        inter_links_1 = pd.merge(link_1, link_2, how='inner')
        inter_links_2 = pd.merge(link_1_cp, link_2, how='inner')
        inter_links = pd.merge(inter_links_1, inter_links_2, how='outer')
        mergedlinks = inter_links.groupby(['source', 'target'], as_index=False).mean().reindex()
        # Graph = nx.from_pandas_edgelist(mergedlinks,
        #                                 source='source',
        #                                 target='target',
        #                                 edge_attr=True)
    elif method == 'snf':

        mergedlinks = _snf_based_merge(link_1, link_2)

    else:

        raise Exception('valid method parameter: union, intersection, knn!')

    return mergedlinks


def graph_traveral(graph, start, threshold, method='bfs'):
    """for given network and start point, it generates a path in a specific manner
        """

    if method == 'bfs':

        res_path = nextra.bfs_tree(G=graph, source=start, depth_limit=threshold)

    elif method == 'dfs':

        res_path = nextra.dfs_tree(G=graph, source=start, depth_limit=threshold)
    else:
        raise Exception('valid method parameter: bfs, dfs!')

    return res_path


def random_walk(gnetdata, start, supervisedby, steps):
    """
    perform supervided random walk for given steps and weights
    """
    path = supervised_random_walk(gnetdata=gnetdata, start=start, supervisedby=supervisedby, steps=steps)
    return path


def path_merge(path_1, path_2, k_mer=3, path='Eulerian'):
    """
    perform de bruijn graph mapping for ginve two path lists
    """
    g = construct_graph([path_1, path_2], k_mer)
    merged_path = output_contigs(g)

    return merged_path
