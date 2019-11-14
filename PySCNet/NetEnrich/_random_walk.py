#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 22:50:56 2019

@author: angelawu
"""

import numpy as np
import pandas as pd

def supervised_random_walk(gnetdata, start, supervisedby = 'pageRank', steps = 10):
        """for given network and node centralities, it performs supervised random walk
        """
        if supervisedby in list(gnetdata.NetAttrs['centralities'].columns)[1:]:
                supervise_matrix = gnetdata.NetAttrs['centralities'][['node', supervisedby]]
        else:
                raise Exception("valid superviseby: betweenness, closeness, degree, pageRank")

        path = [start]
#        links = pd.DataFrame(columns = ('source', 'target', 'weightScore'))
        for i in range(steps):
                neighbors = gnetdata.NetAttrs['graph'][start]
                remove_traverse = [x for x in list(neighbors.keys()) if x not in path]
                new_neighbors = {key: neighbors[key] for key in remove_traverse}

                if len(new_neighbors) == 0:
                        break
                else:
                        new_neighbors_score = {}
                        for key, score in new_neighbors.items():
                                new_score = float(supervise_matrix.loc[supervise_matrix['node'] == key, supervisedby]) * float(list(score.values())[0])
                                new_neighbors_score[key] = new_score
                        go = sorted(new_neighbors_score.items(), key=lambda x: x[1], reverse=True)[0][0]
#                        links.loc[links.shape[0] + 1] = [start, go, new_neighbors_score[go]]
                        path.append(go)
                        start = go
        return(path)



# =============================================================================
#   To do list: random walk based on starting and ending point
# =============================================================================
