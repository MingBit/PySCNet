#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 15:11:22 2019

@author: mwu
"""

import networkx.algorithms.traversal as nextra
from NetEnrich import de_bruijn_raw as debruijn


def breadth_first(graph, start, threshold):
        """ for given network and start point, it generates a path in breadth first search manner
        """
        return nextra.bfs_tree(graph, source = start, depth_limit = threshold)

def depth_first(graph, start, threshold):
        """for give network and start point, it returns a path in depth first search manner
        """
        return nextra.dfs_tree(graph, source = start, depth_limit = threshold)

def de_bruijn_merge(path_1, path_2, k_mer=3, path = 'Eulerian'):

        """ perform de bruijn graph mapping for ginve two path lists
        """
        g = debruijn.construct_graph([path_1, path_2], k_mer)
        merged_path = debruijn.output_contigs(g)
        return merged_path
