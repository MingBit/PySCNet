#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 22:49:50 2019

@author: angelawu
"""

import networkx as nx
import numpy as np
import pandas as pd
import community
import matplotlib.pyplot as plt
import seaborn as sns


class PyNet():
        """ for given a network, it returns node centralities
        parameters: source, target,
        """
        def __init__(self, link, params = None):
                self.link = link
                self.Graph = None
                self.params = {} if params == None else params
                self.betweenness = None
                self.closeness = None
                self.degree = None
                self.pageRank = None
                self.communities = None

                return None

        def get_centrality(self, source, target):
                """ returns betweeness, closeness, degree and pageRank
                """
                self.params['source'] = source
                self.params['target'] = target

                G = nx.from_pandas_edgelist(self.link,
                                            source = source,
                                            target= target,
                                            edge_attr = True)

                self.Graph = G
                self.betweenness = pd.DataFrame.from_dict(list(nx.betweenness_centrality(G).items()))
                self.betweenness.columns = ['node', 'betweenness']
                self.closeness = pd.DataFrame.from_dict(list(nx.closeness_centrality(G).items()))
                self.closeness.columns = ['node', 'closeness']
                self.degree = pd.DataFrame.from_dict(list(nx.degree_centrality(G).items()))
                self.degree.columns = ['node', 'degree']
                self.pageRank = pd.DataFrame.from_dict(list(nx.degree_centrality(G).items()))
                self.pageRank.columns = ['node', 'pageRank']

        def community_detect(self, filepath = None, node_size = 50, font_size = 20):
                """return predicted communities
                """
                self.params['filepath'] = filepath
                subpara = {}
                colors = sns.color_palette() + sns.color_palette('Paired', 100)
                partition = community.best_partition(self.Graph, **subpara)
                self.communities = pd.DataFrame.from_dict(list(partition.items()))
                self.communities.columns = ['node', 'group']
#                size = len(set(partition.values()))
                pos = nx.spring_layout(self.Graph)
                count = 0.
                print(set(partition.values()))
                for com in set(partition.values()) :
                    count = count + 1.
                    list_nodes = [nodes for nodes in partition.keys()
                                                if partition[nodes] == com]
                    nx.draw_networkx_nodes(self.Graph, pos,
                                                   list_nodes, node_size = node_size,
                                                   node_color = colors[com])


                nx.draw_networkx_edges(self.Graph, pos, alpha=0.8)
                labels = dict(zip(list(partition.keys()), list(partition.keys())))

                nx.draw_networkx_labels(self.Graph, pos, labels = labels, font_size = font_size)
                plt.savefig(filepath)
                plt.show()
                plt.close()








