import random
import numpy as np
import pandas as pd
import pyscnet.NetEnrich as ne


def __update_object(gnetdata, grn_method, top_links, resolution=0.5):
    new_object = ne.buildnet(gnetdata, key_links=grn_method, top=int(top_links))
    new_object = ne.get_centrality(new_object)
    new_object = ne.detect_community(new_object, resolution=resolution)

    return new_object


def __update_filter_link(gnetdata, grn_method, top_links, resolution=0.5):
    global new_object
    new_object = __update_object(gnetdata, grn_method, top_links, resolution)
    filtered_link = new_object.NetAttrs[grn_method].sort_values('weight', ascending=False).head(top_links)
    gene_module = new_object.NetAttrs['communities']

    color = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(len(np.unique(gene_module.group)))]

    gene_module['color'] = [color[i] for i in gene_module.group]
    gene_module = pd.concat([new_object.NetAttrs['centralities'], gene_module[['color', 'group']].reset_index(drop=True)],
                            axis=1)
    new_object.NetAttrs['communities'] = gene_module
    nodes = [{'data': {'id': name, 'label': name, 'betweenness': betweenness, 'closeness': closeness,
                       'degree': degree, 'pageRank': pageRank, 'color': color, 'group': group}} for
             name, betweenness, closeness, degree, pageRank, color, group in
             list(gene_module.itertuples(index=False, name=None))]
    edges = [{'data': {'source': source, 'target': target, 'weight': weight}} for source, target, weight in
             list(filtered_link.itertuples(index=False, name=None))]

    new_elements = nodes + edges

    return new_elements


def __update_sub_network(click_node=None):
    click_node = list(new_object.NetAttrs['graph'].nodes)[0] if click_node is None else click_node
    neighbours = list(new_object.NetAttrs['graph'].neighbors(click_node))
    gene_module = new_object.NetAttrs['communities']

    sub_link = pd.DataFrame({'source': np.repeat(click_node, len(neighbours)),
                             'target': neighbours})
    sub_gene_module = gene_module[gene_module.node.isin([click_node] + neighbours)][['node', 'color']]

    sub_nodes = [{'data': {'id': name, 'label': name, 'color': color}} for name, color in
                 sub_gene_module.itertuples(index=False, name=None)]
    sub_edges = [{'data': {'source': source, 'target': target}} for source, target in
                 list(sub_link.itertuples(index=False, name=None))]

    sub_element = sub_nodes + sub_edges

    color = gene_module.loc[gene_module.node == click_node, 'color'].to_list()
    sub_module = gene_module[gene_module.color == color[0]][['node', 'color']]

    inter_node = set(neighbours) & set(sub_module)
    sub_nodes_2 = [{'data': {'id': name, 'label': name, 'color': color}} for name, color in
                   sub_module.itertuples(index=False, name=None)]
    sub_edges_2 = [{'data': {'source': source, 'target': target}}
                   for source, target in
                   list(sub_link.loc[sub_link.target.isin(inter_node)].itertuples(index=False, name=None))]

    sub_element_module = sub_nodes_2 + sub_edges_2

    return [[click_node] + neighbours, sub_element, sub_element_module]

