#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 16:26:19 2020

@author: mwu
"""

import networkx as nx
import numpy as np
import gensim
import pandas as pd
import sklearn.metrics as sm
from tqdm import tqdm
from functools import reduce
import multiprocessing as mp
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Input, Dense, Conv2D, Flatten, MaxPool2D
import tensorflow.keras
from sklearn.model_selection import train_test_split
from sklearn.manifold import MDS


def _build_coexp_graph(ExpMatrix, method = 'pearson', return_graph = False):
    
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
    #keep abosulte value of weight???
    edgeList['weight'] = abs(edgeList['weight'])
    
    if return_graph:
        return (nx.from_pandas_edgelist(edgeList, 'source', 'target', edge_attr = 'weight'))
    else:
        return(edgeList)
    


def _biased_randomWalk(graph, walk_len, num_walks, p = 1, q = 1):
    
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
            
            while(vec_len < walk_len-1):
    
                neighours = list(nx.neighbors(graph, vec_tmp[-1]))
                prob_list = list()
                
                for node in neighours:
                    if len(vec_tmp) < 2:
                        alpha = 1/q
                    else:
                        alpha = 1/p if node == vec_tmp[-2] else 1 if node in list(nx.neighbors(graph, vec_tmp[-2])) else 1/q
                    
    
                    prob = float(alpha) * abs(graph.get_edge_data(node, vec_tmp[-1])['weight']) * abs(pagerank_dict['node'])
                    prob_list.append(prob)
                    
    #            next_node = neighours[prob_list.index(max(prob_list))]
                
                vec_tmp.append(np.random.choice(neighours, size = 1, 
                                                p = [x / sum(prob_list) for x in prob_list])[0])
                
                vec_len = vec_len + 1
                
            vec_each_node.append(vec_tmp)
                
        vec_nodes = vec_nodes + vec_each_node
        
    return(vec_nodes)

def _biased_randomWalk_2(args):
    graph, start_node, walk_len, num_walks, p, q = args
     
    """biased random walk to generate multiple vectors from high-dimentional graph
    p, q: hyperparameters guiding the random walk
    walk_len: the length of vector
    num_walks: the number of vectors
    """
#    pagerank_dict = dict(zip(graph.node, pd.DataFrame.from_dict(list(nx.pagerank_scipy(graph).items()))[1]))
    outer = tqdm(total = num_walks, desc = 'nodes', position = 1)
    vec_nodes = list()
    for i in range(num_walks):
        vec_tmp = list()
#            start_node = np.random.choice(list(graph.nodes))
        vec_len = 0
        vec_tmp.append(start_node)
#        inner = tqdm(total = walk_len, desc = 'walk length', position = 1)
        while(vec_len < walk_len-1):

            neighours = list(nx.neighbors(graph, vec_tmp[-1]))
            prob_list = list()
            
            for node in neighours:
                if len(vec_tmp) < 2:
                    alpha = 1/q
                else:
                    alpha = 1/p if node == vec_tmp[-2] else 1 if node in list(nx.neighbors(graph, vec_tmp[-2])) else 1/q
                

                prob = float(alpha) * abs(graph.get_edge_data(node, vec_tmp[-1])['weight'])
                prob_list.append(prob)
                
#            next_node = neighours[prob_list.index(max(prob_list))]
#            inner.update(1)
            vec_tmp.append(np.random.choice(neighours, size = 1, replace = False,
                                            p = [x / sum(prob_list) for x in prob_list])[0])        
            vec_len = vec_len + 1            
        outer.update(1)    
        vec_nodes.append(vec_tmp)
    
    return(vec_nodes)

def _build_NN_Model(vector_list, dimensions, **parameter_list) -> gensim.models.Word2Vec:
    """build neural network for graph embedding
    vector_list: a list vectors generated from Random walk
    layers: number of layers in nerual network model
    nodes: number of nodes in each layer
    """
    return(gensim.models.Word2Vec(vector_list, size = dimensions, **parameter_list))
    
def _x_y_generator(Expr_df, gene, cell = None, test_size = .6):

    if(cell is None):
        return train_test_split(Expr_df.drop(columns = gene), Expr_df[gene], test_size = test_size)

    elif(gene is None):
        return train_test_split(Expr_df.drop(index = cell).transpose(), Expr_df.loc[cell], test_size = test_size)

    else:
        print('check if gene or cell is NONE!')
        
        

def _build_DNN_Model(Expr, vector_list, dimensions, ref_list, **parameter_list):
    """build deep neural network to convert vectors from random walk
    """
    batch_size = 128
    epochs = 100
    encoding_dim = 32
    
    num_genes = Expr.shape[0]
    embedding = MDS(1)
    gene_var = dict(zip(Expr.index, embedding.fit_transform(Expr))) 
#    vector_df = np.array(list([gene_var[i] for i in sub_list] for sub_list in vector_list)).reshape(300,1000,15,1)
    vector_df = np.array(list([gene_var[i] for i in sub_list] for sub_list in vector_list)).reshape(300,1000,15)
    ref_list = matrix
    x_train, x_test, y_train, y_test = train_test_split(vector_df, ref_list, test_size = .2)
    

    model = Sequential()

#########################build CNN##################################
#    model.add(Conv2D(32, kernel_size=(3, 3),
#                     activation='relu',
#                     input_shape=(1000,15,1)))
#    
#    model.add(Conv2D(64, (3, 3), activation='relu'))
#    model.add(MaxPooling2D(pool_size=(2, 2)))
#    model.add(Dropout(0.25))
#    model.add(Flatten())
#    model.add(Dense(128, activation='relu'))
#    model.add(Dropout(0.5))
#    model.add(Dense(300, activation='softmax'))

#########################build fully connected NN##################################
    model.add(Dense(128, input_dim = (15), activation = 'relu'))
    model.add(Dense(64, activation = 'relu'))
    model.add(Dense(32, activation = 'relu'))
    model.add(Dense(16, activation = 'sigmoid'))

    model.compile(loss=tensorflow.keras.losses.categorical_crossentropy,
                  optimizer=tensorflow.keras.optimizers.Adadelta(),
                  metrics=['accuracy'])
    
    model.fit(x_train, y_train,
              batch_size=batch_size,
              epochs=epochs,
              verbose=1,
              validation_data=(x_test, y_test))
    score = model.evaluate(x_test, y_test, verbose=0)
    print('Test loss:', score[0])
    print('Test accuracy:', score[1])
    
    
def run_node2vec(Expr, method = 'pearson', p = 1, q = 1, 
                walk_len = 20, num_walks = 100, size = 32, 
                workers = 4, **parameters_list):
    """call above functions to build node2vec """
    
    graph = _build_coexp_graph(Expr, method = method, return_graph = True)
    pool = mp.Pool(workers)
    num_of_nodes = len(graph.nodes)
    
    graph_list = [graph for _ in range(num_of_nodes)]
    walk_len_list = [walk_len for _ in range(num_of_nodes)]
    num_walks_list = [num_walks for _ in range(num_of_nodes)]
    p_list = [p for _ in range(num_of_nodes)]
    q_list = [q for _ in range(num_of_nodes)]
    node_list = list(graph.nodes)
    
        
    vector_list = reduce(lambda x,y: x+y, pool.map(_biased_randomWalk_2, 
                                                   zip(graph_list, node_list,
                                                       walk_len_list, num_walks_list, 
                                                       p_list, q_list)))
    
#    pdb.set_trace()
#    node2vec_model = _build_NN_Model(vector_list, dimensions = size, **parameters_list)
#    node_vector = dict()
#    for node in list(graph.nodes):
#        node_vector[node] = node2vec_model.wv.get_vector(node)
#    
#    node2vec_corr = _build_coexp_graph(pd.DataFrame.from_dict(node_vector).T)
#    pool.close()
    return(vector_list)
    

Sim_gne = gnetdata.Gnetdata(Sim)    
vector_list = run_node2vec(Sim, walk_len = 15, num_walks = 1000, size = 5, workers = 12)

import _pickle as pk
with open('/home/mwu/G300_vector_list.pk', 'wb') as f:
    pk.dump(vector_list, f)
    
    
matrix = pd.DataFrame(np.zeros((Sim.shape[0], Sim.shape[0])), index = list(Sim.index), columns = list(Sim.index))

for i in range(Sim_Ref.shape[0]):
    matrix[Sim_Ref.loc[i].node1][Sim_Ref.loc[i].node2] = 1
    matrix[Sim_Ref.loc[i].node2][Sim_Ref.loc[i].node1] = 1

matrix = (np.array(matrix)).reshape(300,300)

 
    
  
