from __future__ import absolute_import
import scnode2vec_raw as scnv
import _pickle as pk
import pandas as pd
import inspect


def run_scnode2vec(Expr, reference_links, p, q, size, walk_len, num_walks,
                   workers, n_pc, param_grid, filename='links.txt', use_ref=True, **kwargs):

    run_randmwalk_args = [k for k, v in inspect.signature(scnv.run_randmwalk).parameters.items()]
    run_randmwalk_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in run_randmwalk_args}

    binary_classifier_args = [k for k, v in inspect.signature(scnv.binary_classifier).parameters.items()]
    binary_classifier_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in binary_classifier_args}

    node_matrix = scnv.run_randmwalk(Expr, p=p, q=q, size=size, walk_len=walk_len,
                                     num_walks=num_walks, workers=workers, n_comp=n_pc, **run_randmwalk_dict)
    if use_ref:
        links = scnv.binary_classifier(embedded_node=node_matrix, reference_links=reference_links,
                                       select_n=0, use_ref=True, param_grid=param_grid, **binary_classifier_dict)
    else:
        links = scnv.binary_classifier(embedded_node=node_matrix, reference_links=reference_links,
                                       select_n=int(reference_links.shape[0] * 0.25), param_grid=param_grid, **binary_classifier_dict)

    links.to_csv(filename, sep='\t', index=False, header=False)


Expr = pd.read_csv('Expr.txt', sep='\t', header=0, index_col=0)

with open('paras.pk', 'rb') as parameter:
    paras = pk.load(parameter)['parameters']

run_scnode2vec(Expr, reference_links=paras['reference_links'], p=paras['p'], q=paras['q'],
               size=paras['size'], walk_len=paras['walk_len'], num_walks=paras['num_walks'],
               workers=paras['workers'], n_pc=paras['n_pc'], param_grid=paras['param_grid'])
