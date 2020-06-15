#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:52:31 2019

@author: angelawu
"""
import networkx as nx
import numpy as np
import gensim
import pandas as pd
import sklearn.metrics as sm
from tqdm import tqdm, trange
from functools import reduce, partial
import multiprocessing as mp
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.decomposition import PCA



def _build_coexp_graph(ExpMatrix, method='pearson', return_graph=False):
    if method == 'cosine':
        corrMatrix = pd.DataFrame(sm.pairwise.cosine_similarity(ExpMatrix), columns=ExpMatrix.index,
                                  index=ExpMatrix.index)
    else:
        corrMatrix = ExpMatrix.T.corr(method=method)

    corrMatrix[corrMatrix == 1] = 0
    # convert adjMatrix to edgelist
    List = list()
    for source in corrMatrix.index.values:
        for target in corrMatrix.index.values:
            if (source != target):
                List.append((target, source, corrMatrix[source][target]))

    edgeList = pd.DataFrame(List, columns=['source', 'target', 'weight'])
    # keep abosulte value of weight???
    #    edgeList['weight'] = abs(edgeList['weight'])

    if return_graph:
        return (nx.from_pandas_edgelist(edgeList, 'source', 'target', edge_attr='weight'))
    else:
        return (edgeList)


def _biased_randomWalk(graph, walk_len, num_walks, p=1, q=1):
    """biased random walk to generate multiple vectors from high-dimentional graph
    p, q: hyperparameters guiding the random walk
    walk_len: the length of vector
    num_walks: the number of vectors
    """
    pagerank_dict = dict(zip(graph.node, pd.DataFrame.from_dict(list(nx.pagerank_numpy(graph).items()))[1]))
    vec_nodes = list()

    for start_node in list(graph.nodes):
        vec_each_node = list()

        for i in range(num_walks):
            vec_tmp = list()
            #            start_node = np.random.choice(list(graph.nodes))
            vec_len = 0
            vec_tmp.append(start_node)

            while (vec_len < walk_len - 1):

                neighours = list(nx.neighbors(graph, vec_tmp[-1]))
                prob_list = list()

                for node in neighours:
                    if len(vec_tmp) < 2:
                        alpha = 1 / q
                    else:
                        alpha = 1 / p if node == vec_tmp[-2] else 1 if node in list(
                            nx.neighbors(graph, vec_tmp[-2])) else 1 / q

                    prob = float(alpha) * abs(graph.get_edge_data(node, vec_tmp[-1])['weight']) * abs(
                        pagerank_dict['node'])
                    prob_list.append(prob)

                #            next_node = neighours[prob_list.index(max(prob_list))]

                vec_tmp.append(np.random.choice(neighours, size=1,
                                                p=[x / sum(prob_list) for x in prob_list])[0])

                vec_len = vec_len + 1

            vec_each_node.append(vec_tmp)

        vec_nodes = vec_nodes + vec_each_node

    return (vec_nodes)


def _biased_randomWalk_2(args):
    Expr, graph, start_node, walk_len, num_walks, p, q = args

    """biased random walk to generate multiple vectors from high-dimentional graph
    p, q: hyperparameters guiding the random walk
    walk_len: the length of vector
    num_walks: the number of vectors
    """

    gene_var = Expr.mean(1)
    outer = tqdm(total=num_walks, desc='nodes', position=1)
    vec_nodes = list()
    for i in range(num_walks):
        vec_tmp = list()
        #            start_node = np.random.choice(list(graph.nodes))
        vec_len = 0
        vec_tmp.append(start_node)
        #        inner = tqdm(total = walk_len, desc = 'walk length', position = 1)
        while (vec_len < walk_len - 1):

            neighours = list(nx.neighbors(graph, vec_tmp[-1]))
            prob_list = list()

            for node in neighours:
                if len(vec_tmp) < 2:
                    alpha = 1 / q
                else:
                    alpha = 1 / p if node == vec_tmp[-2] else 1 if node in list(
                        nx.neighbors(graph, vec_tmp[-2])) else 1 / q

                prob = float(alpha) * abs(graph.get_edge_data(node, vec_tmp[-1])['weight']) * gene_var[node]
                prob_list.append(prob)

            #            next_node = neighours[prob_list.index(max(prob_list))]
            #            inner.update(1)
            vec_tmp.append(np.random.choice(neighours, size=1, replace=False,
                                            p=[x / sum(prob_list) for x in prob_list])[0])
            vec_len = vec_len + 1
        outer.update(1)
        vec_nodes.append(vec_tmp)

    return (vec_nodes)


def _build_NN_Model(vector_list, dimensions, **parameter_list) -> gensim.models.Word2Vec:
    """build neural network for graph embedding
    vector_list: a list vectors generated from Random walk
    layers: number of layers in nerual network model
    nodes: number of nodes in each layer
    """
    return (gensim.models.Word2Vec(vector_list, size=dimensions, **parameter_list))


def _binary_classifier(embedded_node, reference_links):
    """
    build binary classifier for embedded node vector
    embedded_node: a dataframe of nodes with embedded vector
    reference_links: reference adjacency matrix of nodes
    -------
    return: the probabilities of binary classification
    """

    embedded_edges = pd.DataFrame(list(embedded_node.loc[reference_links.iloc[i]['source']]) +
                                  list(embedded_node.loc[reference_links.iloc[i]['target']])
                                  for i in range(reference_links.shape[0]))

    train_features, test_features, train_targets, test_targets = train_test_split(
        embedded_edges, reference_links['weight'].abs(),
        train_size=0.8,
        test_size=0.2,
        random_state=1,
        stratify=reference_links['weight'].abs())

    param_grid = {
        'n_estimators': [50, 100, 500],
        'max_depth': [10, 20, 50],
        'max_features': [3, 5],
        'min_samples_leaf': [3, 4, 5],
        'criterion': ['gini', 'entropy']
    }
    classifier = RandomForestClassifier(random_state=1)

    CV_rfc = GridSearchCV(estimator=classifier, param_grid=param_grid, cv=10, n_jobs=10)
    CV_rfc.fit(train_features, train_targets)

    prediction_training_targets = CV_rfc.predict(train_features)
    self_accuracy = accuracy_score(train_targets, prediction_training_targets)
    print("Accuracy for training data (self accuracy):", self_accuracy)

    # predict the 'target' for 'test data'
    prediction_test_targets = CV_rfc.predict(test_features)
    test_accuracy = accuracy_score(test_targets, prediction_test_targets)
    print("Accuracy for test data:", test_accuracy)

    prediction_prob = CV_rfc.predict_proba(embedded_edges)
    prediction_type = CV_rfc.predict(embedded_edges)

    prediction_all_df = pd.DataFrame({'source': reference_links.source,
                                      'target': reference_links.target,
                                      'weight': [(1 - prediction_prob[i][0]) for i in range(len(prediction_type))]})
    return (prediction_all_df)


def run_node2vec(Expr, reference_links, method,  p, q, walk_len, num_walks, size,
                 workers, n_pc, **parameters_list):
    """call above functions to build node2vec """
    print('-----running PCA')
    pca = PCA(n_components=n_pc)
    pca.fit(Expr.T)
    top_pcs = pd.DataFrame(pca.components_).T

    print('-----running node2vec---------')
    graph = _build_coexp_graph(Expr, method=method, return_graph=True)
    pool = mp.Pool(workers)
    num_of_nodes = len(graph.nodes)

    graph_list = [graph for _ in range(num_of_nodes)]
    walk_len_list = [walk_len for _ in range(num_of_nodes)]
    num_walks_list = [num_walks for _ in range(num_of_nodes)]
    p_list = [p for _ in range(num_of_nodes)]
    q_list = [q for _ in range(num_of_nodes)]
    Expr = [Expr for _ in range(num_of_nodes)]
    node_list = list(graph.nodes)

    vector_list = reduce(lambda x, y: x + y, pool.map(_biased_randomWalk_2,
                                                      zip(Expr, graph_list, node_list,
                                                          walk_len_list, num_walks_list,
                                                          p_list, q_list)))

    node2vec_model = _build_NN_Model(vector_list, dimensions=size, **parameters_list)
    node_vector = dict()
    for node in list(graph.nodes):
        node_vector[node] = node2vec_model.wv.get_vector(node)

    node_vector = pd.DataFrame.from_dict(node_vector).T
    top_pcs.index = node_vector.index
    merge_node_vector = pd.concat([node_vector, top_pcs], axis=1)

    node2vec_prob = _binary_classifier(merge_node_vector, reference_links)
    pool.close()

    return (node2vec_prob)
