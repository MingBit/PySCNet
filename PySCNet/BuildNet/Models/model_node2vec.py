#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:24:10 2019

@author: mwu
"""

import numpy as np
import networkx as nx
from node2vec import Node2Vec
from node2vec.edges import HadamardEmbedder
import pandas as pd
import sklearn.metrics as sm

def _build_coexp_graph(ExpMatrix, method = 'pearson'):
    if method == 'cosine':
        corrMatrix = pd.DataFrame(sm.pairwise.cosine_similarity(ExpMatrix), columns = ExpMatrix.index, index = ExpMatrix.index)
    else:
        corrMatrix = ExpMatrix.T.corr(method = method)
    
    corrMatrix[corrMatrix == 1] = 0
    #convert adjMatrix to edgelist
    List = list()
    for source in corrMatrix.index.values:
        for target in corrMatrix.index.values:
            if(source != target):
                List.append((target,source,corrMatrix[source][target]))
    
    edgeList = pd.DataFrame(List, columns = ['source', 'target', 'weight'])
    return(edgeList)
     


def _build_node2vec_model(edgeList, method = 'pearson', **kwarg):
    #build weighted graph
    G = nx.from_pandas_edgelist(edgeList, 'source', 'target', edge_attr = 'weight')

    #pre-compute the probabilities and generate walks
    node2vec = Node2Vec(G, **kwarg)
    
    #embed the nodes: https://radimrehurek.com/gensim/models/word2vec.html
    #windon = Maximum distance between the current and predicted word within a sentence
    #min_count = Ignores all words with total frequency lower than this.
    model = node2vec.fit(window = 10, min_count = 1, batch_words = 4)
    node_vector = dict()
    for node in list(G.nodes):
        node_vector[node] = model.wv.get_vector(node)
    
    node2vec_corr = _build_coexp_graph(pd.DataFrame.from_dict(node_vector).T)
    return(node2vec_corr)
    
 
    
    

