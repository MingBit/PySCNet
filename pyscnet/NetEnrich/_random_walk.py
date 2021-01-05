#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 22:50:56 2019

@author: angelawu
"""
import numpy as np
import pandas as pd


def greedy_walk(gnetdata, start, supervisedby, steps):
    """
    for given network and node centralities, it performs greedy walk guided by node centrality
    """

    if supervisedby in list(gnetdata.NetAttrs['centralities'].columns)[1:]:
        supervise_matrix = gnetdata.NetAttrs['centralities'][['node', supervisedby]]
    else:
        raise Exception("valid superviseby: betweenness, closeness, degree, pageRank")

    path = [start]
    for i in range(steps):
        neighbors = gnetdata.NetAttrs['graph'][start]
        remove_traverse = [x for x in list(neighbors.keys()) if x not in path]
        new_neighbors = {key: neighbors[key] for key in remove_traverse}

        if len(new_neighbors) == 0:
            break
        else:
            new_neighbors_score = {}
            for key, score in new_neighbors.items():
                new_score = float(supervise_matrix.loc[supervise_matrix['node'] == key, supervisedby]) * float(
                    list(score.values())[0])
                new_neighbors_score[key] = new_score
            go = sorted(new_neighbors_score.items(), key=lambda x: x[1], reverse=True)[0][0]
            path.append(go)
            start = go
    return path


def supervised_random_walk(gnetdata, start, supervisedby, steps, repeat=1000):
    """
    for given network and node centralities, it performs supervised random walk
    """

    if supervisedby in list(gnetdata.NetAttrs['centralities'].columns)[1:]:
        supervise_matrix = gnetdata.NetAttrs['centralities'][['node', supervisedby]]
    else:
        raise Exception("valid superviseby: betweenness, closeness, degree, pageRank")
    path_list = list()

    n = 0
    while n < repeat:
        n += 1
        path = [start]
        start_node = start
        for i in range(steps):
            neighbors = gnetdata.NetAttrs['graph'][start_node]
            remove_traverse = [x for x in list(neighbors.keys()) if x not in path]

            rand_pro = dict(zip(remove_traverse, list(np.random.dirichlet(np.ones(len(remove_traverse))))))
            new_neighbors = {key: (neighbors[key], rand_pro[key]) for key in remove_traverse}

            if len(new_neighbors) == 0:
                break
            else:
                new_neighbors_score = {}
                for key, score in new_neighbors.items():
                    new_score = float(supervise_matrix.loc[supervise_matrix['node'] == key, supervisedby]) * float(
                        list(score[0].values())[0]) * float(score[1])
                    new_neighbors_score[key] = new_score

                go = sorted(new_neighbors_score.items(), key=lambda x: x[1], reverse=True)[0][0]
                path.append(go)
                start_node = go

        path_list.append(path)

    path_list = pd.DataFrame.from_records(path_list)
    final_path = [start]

    for i in path_list.columns[1:]:
        tmp = path_list[i].value_counts()
        j = 0
        next_node = tmp.idxmax()
        while j < len(tmp):
            j += 1
            if next_node in final_path or next_node not in list(gnetdata.NetAttrs['graph'][final_path[-1]]):
                if j == len(tmp):
                    break
                else:
                    next_node = tmp.iloc[j:].idxmax()

            else:
                final_path.append(next_node)
                break

    return final_path

