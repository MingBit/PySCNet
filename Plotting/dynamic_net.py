#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:32:32 2019

@author: mwu
"""

import networkx as nx
import pandas as pd
#import matplotlib.pyplot as plt
#from networkx.readwrite import json_graph
from pyvis.network import Network

def dynamic_netShow(netx, filepath, PyNetObj = None, filterby = 'edges'):

        """ it returns a html file representing interactive network
        """
        if PyNetObj is not None:
                link_table = PyNetObj.link
        else:
                tmp = pd.DataFrame.from_dict(nx.get_edge_attributes(netx, 'score'), oriented = 'index', columns=['score'])
                link_table = pd.concat(tmp, pd.DataFrame(list(tmp.index)))

        net = Network()
        #Gdy.from_nx(Net_test.Graph)
        edge_data = zip(link_table['source'], link_table['target'], link_table['score'])

        for e in edge_data:
            src = e[0]
            dst = e[1]
            w = e[2]

            net.add_node(src, src, title=src)
            net.add_node(dst, dst, title=dst)
            net.add_edge(src, dst, value=w)

        neighbor_map = net.get_adj_list()

        # add neighbor data to node hover data
        for node in net.nodes:
            node["title"] += " Neighbors:<br>" + "<br>".join(neighbor_map[node["id"]])
            node["value"] = len(neighbor_map[node["id"]])

        net.show_buttons(filter_=filterby)
        net.show(filepath)



