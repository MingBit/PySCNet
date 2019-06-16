#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 20:50:39 2019

@author: mwu
"""

import networkx as nx
import pandas as pd
import community
import snf
import numpy as np

import networkx.algorithms.traversal as nextra
from NetEnrich import de_bruijn_raw as debruijn
from NetEnrich import random_walk as random_walk


def get_centrality(gedata):

        """ returns betweeness, closeness, degree and pageRank
        """
        G = gedata.NetAttrs['graph']
        centralities = pd.DataFrame(list(G.node), columns=['node'])
        centralities['betweenness'] = pd.DataFrame.from_dict(list(nx.betweenness_centrality(G).items()))[1]
        centralities['closeness'] = pd.DataFrame.from_dict(list(nx.closeness_centrality(G).items()))[1]
        centralities['degree'] = pd.DataFrame.from_dict(list(nx.degree_centrality(G).items()))[1]
        centralities['pageRank'] = pd.DataFrame.from_dict(list(nx.degree_centrality(G).items()))[1]

        gedata.NetAttrs['centralities'] = centralities
        return(gedata)

def community_detect(gedata):

        """return predicted communities
        """
        G = gedata.NetAttrs['graph']
        subpara = {}
#        colors = sns.color_palette() + sns.color_palette('Paired', 100)
        partition = community.best_partition(G, **subpara)
        communities = pd.DataFrame.from_dict(list(partition.items()))
        communities.columns = ['node', 'group']

        gedata.NetAttrs['communities'] = communities


        return(gedata)

def __linkage_to_adjlink(linkage_table, node_list):
        """convert linkage table to weighted adjacency matrix
        """
        #TODO: debug
        adjlink_matrix = pd.DataFrame(0, columns=node_list, index = node_list, dtype = np.float)
#        source, target, score = list(linkage_table.columns)
        for i in range(0, len(linkage_table)):

                adjlink_matrix.loc[linkage_table['source'][i]][linkage_table['target'][i]] = linkage_table['score'][i]
                adjlink_matrix.loc[linkage_table['target'][i]][linkage_table['source'][i]] = linkage_table['score'][i]

        return np.array(adjlink_matrix)



def __knn_based_merge(link_1, link_2):

        node_list = list(set(link_1['source'], link_1['target']) & set(link_2['source'], link_2['target']))
        adjlink_1 = __linkage_to_adjlink(link_1, node_list)
        adjlink_2 = __linkage_to_adjlink(link_2, node_list)
        adjlinks = list()
        adjlinks.append(adjlink_1)
        adjlinks.append(adjlink_2)
        affinity_matrix = snf.make_affinity(adjlinks)
        fused_network = snf.snf(affinity_matrix)
        Graph = nx.from_numpy_matrix(fused_network)

        return(Graph)



def graph_merge(link_1, link_2, method = 'union'):

        """it returns the merged network
        """
        if method == 'union':
                union_links = pd.merge(link_1, link_2, how = 'outer')
                mergedlinks = union_links.groupby(['source', 'target'], as_index = False).mean().reindex()
                Graph = nx.from_pandas_edgelist(mergedlinks,
                                    source = 'source',
                                    target= 'target',
                                    edge_attr = True)
        elif method == 'intersection':
                inter_links = pd.merge(link_1, link_2, how = 'inner')
                mergedlinks = inter_links.groupby(['source', 'target'], as_index = False).mean().reindex()
                Graph = nx.from_pandas_edgelist(mergedlinks,
                                    source = 'source',
                                    target= 'target',
                                    edge_attr = True)
        elif method == 'knn':
                #TODO: debug
                Graph = __knn_based_merge(link_1, link_2)

        else:
            print('recheck the parameter method')

        return(Graph)


def graph_traveral(graph, start, threshold, method = 'bfs'):
        """for given network and start point, it generates a path in a specific manner
        """

        if method == 'bfs':
                return(nextra.bfs_tree(graph, source=start, depth_limit=threshold))

        elif method == 'dfs':
                return nextra.dfs_tree(graph, source = start, depth_limit = threshold)

        #TODO: debug
        elif method == 'random':

                graph = random_walk.supervised_random_walk(graph, start = start, steps = threshold)
                return(graph)



def path_merge(path_1, path_2, k_mer = 3, path = 'Eulerian'):
        """ perform de bruijn graph mapping for ginve two path lists
        """
        g = debruijn.construct_graph([path_1, path_2], k_mer)
        merged_path = debruijn.output_contigs(g)

        return merged_path









