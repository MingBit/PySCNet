#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:50:35 2019

@author: mwu
"""

import pandas as pd
import snf
import numpy as np
import networkx as nx
import scipy
import pdb


def linkage_to_adjlink(linkage_table, node_list):
        """convert linkage table to weighted adjacency matrix
        """
        adjlink_matrix = pd.DataFrame(0, columns=node_list, index = node_list, dtype = np.float)
#        source, target, score = list(linkage_table.columns)
        for i in range(0, len(linkage_table)):
                adjlink_matrix.loc[linkage_table['source'][i]][linkage_table['target'][i]] = linkage_table['score'][i]
                adjlink_matrix.loc[linkage_table['target'][i]][linkage_table['source'][i]] = linkage_table['score'][i]

        return np.array(adjlink_matrix)


class MergedNetwork:
        """ build merged network for given two linkage tables or paths
        """
        def __init__(self, link_1, link_2, node_list, params = None):
                self.link_1 = link_1
                self.link_2 = link_2
                self.node_list = node_list
#                self.mergedlinks = None
                self.params = params
                self.Graph = None

                return None

        def network_merge(self, source = 'source', target = 'target', method = 'Union'):
            """for given two networks, it returns a merged network
            """
        #    merged_links = pd.DataFrame(columns = ['source', 'target', 'weightScore'])
            if method == 'Union':
                union_links = pd.merge(self.link_1, self.link_2, how = 'outer')
                mergedlinks = union_links.groupby([source, target], as_index = False).mean().reindex()
                self.Graph = nx.from_pandas_edgelist(mergedlinks,
                                            source = source,
                                            target= target,
                                            edge_attr = True)
            elif method == 'Intersection':
                inter_links = pd.merge(self.link_1, self.link_2, how = 'inner')
                mergedlinks = inter_links.groupby(['source', 'target'], as_index = False).mean().reindex()
                self.Graph = nx.from_pandas_edgelist(mergedlinks,
                                            source = source,
                                            target= target,
                                            edge_attr = True)
#            elif method == 'de bruijn':
#                merged_links = pd.merge()
            else:
                    print('recheck the parameter method')

        def knn_based_merge(self, **kwarg):

                """for given two linkage tables, merge networks via KNN (snf)
                """
                adjlink_1 = linkage_to_adjlink(self.link_1, self.node_list)
                adjlink_2 = linkage_to_adjlink(self.link_2, self.node_list)
                adjlinks = list()
                adjlinks.append(adjlink_1)
                adjlinks.append(adjlink_2)
                affinity_matrix = snf.make_affinity(adjlinks, **kwarg)
                fused_network = snf.snf(affinity_matrix)
                self.Graph = nx.from_numpy_matrix(fused_network)
#                fused_linkage = np.matrix(fused_network.data)
#
#                x, y = scipy.where(fused_linkage > filter_by)
#                fused_linkage[:,:] = 0
#                fused_linkage[x,y] = 1
#                self.mergedlinks = fused_linkage
#                merged_network = nx.from_numpy_matrix(fused_linkage)
#                self.Graph = merged_network







