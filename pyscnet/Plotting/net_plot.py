#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mwu
"""

from __future__ import absolute_import
from ._nx2gt import nx2gt
import graph_tool.all as gt
import seaborn as sns
import pandas as pd
import numpy as np
from pyvis.network import Network
import matplotlib.pyplot as plt

from sklearn import preprocessing
from .__dash_network import __update_filter_link, __update_sub_network, __update_object

from jupyter_plotly_dash import JupyterDash
import dash_cytoscape as cyto
import dash_html_components as html
from dash.dependencies import Input, Output


def geneHeatmap(gnetdata, gene, cell_clusterid, select_by, order_by=None, scale=True, cmap='RdBu', save_as=None,
                **kwargs):
    """
    Create Heatmap showing gene expression in individual cells along pseudotime
    ----------------------------------------------------------------------------
    :param gnetdata: Gnetdata object, default None
    :param gene: list, default None.
    :param cell_clusterid: str, default None. cell with cell_clusterid will be selected
    :param select_by: str, default None. key of filtering cells
    :param order_by: str, default None. key of ordering cells
    :param scale: bool, default True. whether or not scale the data
    :param cmap: str, default 'RdBu'. string denoting colors in clustermap
    :param save_as: str, default None. filepath+filename
    :param kwargs: additional parameters passed to seaborn.clustermap()
    :return: None
    """
    cell_info = gnetdata.CellAttrs['CellInfo'] if order_by is None else gnetdata.CellAttrs['CellInfo'].sort_values(
        order_by, ascending=True)
    cell = list(cell_info.loc[cell_info[select_by].isin([cell_clusterid])].index)
    sub_Expr = pd.DataFrame(preprocessing.scale(gnetdata.ExpMatrix),
                            index=gnetdata.ExpMatrix.index, columns=gnetdata.ExpMatrix.columns).loc[gene, cell] \
        if scale else gnetdata.ExpMatrix.loc[gene, cell]
    sns_plot = sns.clustermap(sub_Expr, cmap=cmap, xticklabels=False, **kwargs)

    if save_as is not None:
        sns_plot.savefig(save_as)

    # plt.show()
    return sns_plot


def geneCorrelation(gnetdata, gene, cell_clusterid, select_by, order_by=None, scale_data=True,
                    save_as=None, figsize=None, **kwargs):
    """
    Create gene correlation heatmap
    --------------------------------------------------
    :param gnetdata: Gnetdata object
    :param gene: list, default None.
    :param cell_clusterid: str, default None. cell with cell_clusterid will be selected
    :param select_by: str, default None. key of filtering cells
    :param order_by: str, default None. key of ordering cells
    :param scale_data: bool, default True. whether or not scale the data
    :param save_as: str, default None. filepath+filename
    :param figsize: list, default None. a list of int defining figure size
    :param kwargs: additional parameters passed to seaborn.clustermap()
    :return: None
    """
    cell_info = gnetdata.CellAttrs['CellInfo'] if order_by is None else gnetdata.CellAttrs['CellInfo'].sort_values(
        order_by, ascending=True)
    cell = list(cell_info.loc[cell_info[select_by].isin([cell_clusterid])].index)
    sub_Expr = pd.DataFrame(preprocessing.scale(gnetdata.ExpMatrix),
                            index=gnetdata.ExpMatrix.index, columns=gnetdata.ExpMatrix.columns).loc[gene, cell].T \
        if scale_data else gnetdata.ExpMatrix.loc[gene, cell].T

    corr = sub_Expr.corr()
    #     mask = np.triu(np.ones_like(corr, dtype=np.bool))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    fig, ax = plt.subplots(figsize=[10, 10] if figsize is None else figsize)
    sns_heatmap = sns.heatmap(corr, cmap=cmap, center=0, **kwargs)

    if save_as is not None:
        plt.savefig(save_as)

    # plt.show()
    return sns_heatmap


def dynamic_netShow(gnetdata, filename, node_size=50, node_community='all', html_size=["500px", "1800px"],
                    bgcolor="#222222", font_color="white", **kwarg):
    """
    create a GRN html
    ------------------------------------------
    :param gnetdata: Gnetdata object
    :param filename: str, save as filename
    :param node_size: int, default 50.
    :param node_community: str or list, default all.
    :param html_size: list, default ["1200px", "1600px"]
    :param font_color: string, default white
    :param bgcolor: string, default #222222
    :return: None
    """
    assert 'graph' in gnetdata.NetAttrs.keys(), 'graph is empty!'
    assert 'communities' in gnetdata.NetAttrs.keys(), 'node communities is empty!'

    node_group = gnetdata.NetAttrs['communities']
    graph = gnetdata.NetAttrs['graph']

    net = Network(html_size[0], html_size[1], bgcolor=bgcolor, font_color=font_color, **kwarg)
    edge_data = zip(list(graph.edges), list(graph[x[0]][x[1]]['weight'] for x in graph.edges))
    colors = sns.color_palette().as_hex() + sns.color_palette('Paired', 100).as_hex()

    for e in edge_data:
        src = e[0][0]
        dst = e[0][1]
        w = e[1]

        if node_community == 'all':
            src_color = int(node_group[node_group['node'] == src]['group'])
            dst_color = int(node_group[node_group['node'] == dst]['group'])
        else:
            sub_node_group = node_group[node_group.group.isin(node_community)]
            src_color = int(sub_node_group[sub_node_group['node'] == src]['group']) if src in list(
                sub_node_group.node) else 'grey'
            dst_color = int(sub_node_group[sub_node_group['node'] == dst]['group']) if dst in list(
                sub_node_group.node) else 'grey'

        net.add_node(src, src, title=src, color='grey' if src_color == 'grey' else colors[src_color])
        net.add_node(dst, dst, title=dst, color='grey' if dst_color == 'grey' else colors[dst_color])
        net.add_edge(src, dst, value=w)

    neighbor_map = net.get_adj_list()

    # add neighbor data to node hover data
    for node in net.nodes:
        node["title"] += " Neighbors:<br>" + "<br>".join(neighbor_map[node["id"]])
        node["value"] = len(neighbor_map[node["id"]])
        node['size'] = node_size

    net.show_buttons(filter_='physics')
    net.show(filename)

    return None


def net_hierarchy_plot(gnetdata, filename=None, **kwarg):
    """
     create a hierarchy gene net plot
    ---------------------------------------------
    :param gnetdata: Gnetdata object
    :param filename: str, default None.
    :param kwarg: additional parameters passed to graph_tool.all.draw_hierarchy()
    :return: None

    """

    assert 'graph' in gnetdata.NetAttrs.keys(), 'graph is empty!'
    assert 'communities' in gnetdata.NetAttrs.keys(), 'node communities is empty!'

    graph = nx2gt(gnetdata.NetAttrs['graph'])
    node_group = gnetdata.NetAttrs['communities']

    # deg = graph.degree_property_map('total')
    ngroup = graph.new_vertex_property('int')

    labels = dict(zip(list(range(graph.num_vertices())), list(graph.vertex_properties['id'])))
    for g in labels.keys():
        ngroup[g] = node_group.loc[node_group.node == labels[g], 'group']

    state = gt.minimize_nested_blockmodel_dl(graph, deg_corr=True)
    gt.draw_hierarchy(state, vertex_fill_color=ngroup, vertex_anchor=0,
                      vertex_text=graph.vertex_properties['id'],
                      output=filename, **kwarg)
    return None


def net_matrix_plot(gnetdata, filename=None, highlight_path=None, **kwarg):
    """
    create a matrix gene net plot
    --------------------------------------
    :param gnetdata: Gnetdata object.
    :param filename: str, default None.
    :param highlight_path: list, default None. a list of gene Nodes.
    :param kwarg: additional parameters passed to graph_tool.all.graph_draw()
    :return: None

    """

    assert 'graph' in gnetdata.NetAttrs.keys(), 'graph is empty!'
    assert 'communities' in gnetdata.NetAttrs.keys(), 'node communities is empty!'

    graph = nx2gt(gnetdata.NetAttrs['graph'])
    node_group = gnetdata.NetAttrs['communities']

    ngroup = graph.new_vertex_property('int')

    # deg = graph.degree_property_map('total')
    labels = dict(zip(list(range(graph.num_vertices())), list(graph.vertex_properties['id'])))
    for g in labels.keys():
        ngroup[g] = node_group.loc[node_group.node == labels[g], 'group']

    position = graph.new_vertex_property("vector<double>")
    edge_color = graph.new_edge_property("string")

    dim = int(np.sqrt(len(labels))) + 1
    node_pos = list()

    node_group_dict = dict(zip(list(range(graph.num_vertices())), list(node_group.group)))

    for i in range(len(set(node_group.group))):
        for key, item in node_group_dict.items():
            if item == i:
                node_pos.append(key)

    for i in range(len(node_pos)):
        position[graph.vertex(node_pos[i])] = (int(i / dim), i % dim)

    for e in graph.edges():
        if (highlight_path is not None) and (labels[list(e)[0]] in highlight_path) and (
                labels[list(e)[1]] in highlight_path) and abs(
            highlight_path.index(labels[list(e)[0]]) - highlight_path.index(labels[list(e)[1]])) == 1:
            edge_color[e] = 'darkorange'
        else:
            edge_color[e] = 'lightgrey'

    gt.graph_draw(graph, pos=position, vertex_fill_color=ngroup, edge_color=edge_color,
                  vertex_text=graph.vertex_properties['id'], output=filename, **kwarg)

    return None


def create_app(gnetdata, grn_method, top_links, resolution=0.5, layout='cose'):
    """
    Create dash-plotly in Jupyter notebook
    ------------------------------------------------------------
    :param gnetdata: Gnetdata object.
    :param grn_method: str, defualt None. It refers to the links table stored in NetAttrs.
    :param top_links: int, default None.
    :param resolution: float, default 0.5. resolution for gene module detection.
    :param layout: str, default cose. Network layout.
    :return: dashboard html

    """
    app = JupyterDash('pyscnet-plotly-dash')
    elements = __update_filter_link(gnetdata, grn_method, top_links, resolution)
    neighbours, sub_element_1, sub_element_2 = __update_sub_network(click_node=None)
    def_text = 'please click on the gene node!'
    FONT_STYLE = {
        "color": '#343a40',
        'font-size': '30'
    }
    new_stylesheet = [
        {
            'selector': 'node',
            'style': {
                'label': 'data(id)',
                'background-color': 'data(color)',
                'color': '#343a40'}
        }]
    app.layout = html.Div([
        html.H3("pyscnet-plotly-dash"),
        cyto.Cytoscape(
            id='gene_network',
            layout={'name': layout},
            style={'width': '100%', 'height': '800px', 'background-color': '#eddcd2'},
            stylesheet=new_stylesheet,
            elements=elements
        ),

        html.H3(id='node_neighbors', children=def_text, style=FONT_STYLE),
        cyto.Cytoscape(
            id='selected_node_neighbors',
            layout={'name': layout},
            style={'width': '100%', 'height': '800px', 'background-color': '#eddcd2'},
            stylesheet=new_stylesheet,
            elements=sub_element_1
        ),

        html.H3(id='node_module', children=def_text, style=FONT_STYLE),
        cyto.Cytoscape(
            id='selected_node_module',
            layout={'name': 'grid'},
            style={'width': '100%', 'height': '800px', 'background-color': '#eddcd2'},
            stylesheet=new_stylesheet,
            elements=sub_element_2)

    ])

    @app.callback([Output('selected_node_neighbors', 'elements'),
                   Output('selected_node_module', 'elements'),
                   Output('selected_node_neighbors', 'stylesheet'),
                   Output('selected_node_module', 'stylesheet'),
                   Output('node_neighbors', 'children'),
                   Output('node_module', 'children')],
                  [Input('gene_network', 'tapNodeData')])
    def update_sub_net(data):
        if data:
            neighbours, new_sub_elements_1, new_sub_elements_2 = __update_sub_network(data['id'])
            new_stylesheet_1 = [{
                'selector': 'node',
                'style': {
                    'label': 'data(id)',
                    'color': '#343a40',
                    'background-color': 'data(color)'
                }
            }]

            neighbour_text = 'Genes connected to ' + data['id']
            module_text = 'Genes assigned to the same module as ' + data['id']

        return [new_sub_elements_1, new_sub_elements_2, new_stylesheet_1,
                new_stylesheet_1, neighbour_text, module_text]

    return app
