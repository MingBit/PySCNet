#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:32:32 2019
@author: mwu
"""

import seaborn as sns
from pyvis.network import Network
import pandas as pd

def dynamic_netShow(gnetdata, filepath):

        """ it returns a html file representing interactive network
        """
        if(gnetdata.NetAttrs['parameters']['threshold'] != 'None'):
                filterby = float(gnetdata.NetAttrs['parameters']['threshold'])
                link_table = gnetdata.NetAttrs['links'].loc[gnetdata.NetAttrs['links'].weight > filterby]
        elif(gnetdata.NetAttrs['parameters']['top'] != 'None'):
                filterby = int(gnetdata.NetAttrs['parameters']['top'])
                link_table = gnetdata.NetAttrs['links'].sort_values('weight', ascending = False).head(filterby)
        else:
                link_table = gnetdata.NetAttrs['links']
        node_group = gnetdata.NetAttrs['communities']

        net = Network("800px", '1600px',bgcolor="#222222", font_color="white")
        edge_data = zip(link_table['source'], link_table['target'], link_table['weight'])

        colors = sns.color_palette().as_hex() + sns.color_palette('Paired', 100).as_hex()

        for e in edge_data:
            src = e[0]
            dst = e[1]
            w = e[2]

            if node_group is not None:
                    src_color = int(node_group[node_group['node'] == src]['group'])
                    dst_color = int(node_group[node_group['node'] == dst]['group'])
            else:
                    src_color = 2
                    dst_color = 2

            net.add_node(src, src, title=src, color = colors[src_color])
            net.add_node(dst, dst, title=dst, color = colors[dst_color])
            net.add_edge(src, dst, value=w)
        neighbor_map = net.get_adj_list()

        # add neighbor data to node hover data
        for node in net.nodes:
            node["title"] += " Neighbors:<br>" + "<br>".join(neighbor_map[node["id"]])
            node["value"] = len(neighbor_map[node["id"]])
            node['size'] = 20

        #net.show_buttons(filter_= 'edges')
        net.show(filepath)

def Genes_modules_Cell_Clusters(gnetdata):
        """it returns summary of genes in each modules"""

        node_group = gnetdata.NetAttrs['communities']
        gene_module = node_group.groupby('group').count()

        return(gene_module)

def static_netShow(graph, filepath, node_group = None):


        net = Network("800px", '1600px',bgcolor="#222222", font_color="white")
        edge_data = graph.edges()
        colors = sns.color_palette().as_hex() + sns.color_palette('Paired', 100).as_hex()

        for e in edge_data:
            src = e[0]
            dst = e[1]


            if node_group is not None:
                    src_color = int(node_group[node_group['node'] == src]['group'])
                    dst_color = int(node_group[node_group['node'] == dst]['group'])
            else:
                    src_color = 2
                    dst_color = 2

            net.add_node(src, src, title=src, color = colors[src_color])
            net.add_node(dst, dst, title=dst, color = colors[dst_color])
            net.add_edge(src, dst, value=1)
        neighbor_map = net.get_adj_list()

        # add neighbor data to node hover data
        for node in net.nodes:
            node["title"] += " Neighbors:<br>" + "<br>".join(neighbor_map[node["id"]])
            node["value"] = len(neighbor_map[node["id"]])
            node['size'] = 20

        #net.show_buttons(filter_= 'edges')
        net.show(filepath)













