#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:51:34 2019

@author: mwu
"""
from __future__ import absolute_import
import docker
import pandas as pd
import os
import networkx as nx
import tarfile
import warnings
import _pickle as pk

global path
path = os.path.join(os.path.dirname(__file__)) + '/Docker_App/'


def _copy_to(container_id, src, dst):
    client = docker.from_env()
    container = client.containers.get(container_id)
    strm, stat = container.get_archive(src)

    with open(os.getenv('HOME') + '/temp.tar', 'w') as f:
        for line in strm:
            f.write(str(line, 'utf-8'))
        f.seek(0)

        thisTar = tarfile.TarFile(f.name)
        thisTar.extract('links.txt', path=dst)
    os.remove(os.getenv('HOME') + '/temp.tar')


def _remove_duplicate(links):
    links_list = sorted(links[['source', 'target']].values.tolist())
    for i in range(len(links_list)):
        links_list[i] = tuple(sorted(links_list[i]))
    nodes = pd.DataFrame(list(set(links_list)), columns=('source', 'target'))
    links = pd.merge(links, nodes, how='right')

    return links


def _rundocker(gnetdata, method, cell_clusterid=None, Mms_TF=None, **kwargs):
    client = docker.from_env()
    # TODO: remove cell_clusterid???
    if cell_clusterid is not None:
        #            pd.DataFrame.to_csv(cell_clusterid, path + method + '/cell_clusterid.txt',  sep= '\t', header=False, index=False)
        cell_info = gnetdata.CellAttrs
        tmp_expr = gnetdata.ExpMatrix[list(cell_info.loc[cell_info.cluster_id.isin([cell_clusterid])].index)]
    else:
        tmp_expr = gnetdata.ExpMatrix

    if Mms_TF is not None:
        pd.DataFrame.to_csv(Mms_TF, path + method + '/Mms_TF.txt', sep='\t', header=False, index=False)

    with open(path + method + '/paras.pk', 'wb') as outfile:
        pk.dump(kwargs, outfile)

    pd.DataFrame.to_csv(tmp_expr, path + method + '/Expr.txt', sep='\t')
    client.images.build(path=path + method, dockerfile='Dockerfile', tag=method.lower())
    container = client.containers.run(method.lower(), detach=True)
    _copy_to(container_id=container.short_id, src='/' + method + '/links.txt', dst=os.getenv('HOME'))

    #        client.remove_container(container.short_id)
    container.stop()
    client.containers.prune()
    client.images.prune()
    os.system('rm ' + path + method + '/Expr.txt | rm ' + path + method + '/paras.pk')
    raw_links = pd.read_csv(os.getenv('HOME') + '/links.txt', sep='\t', header=None)
    raw_links.columns = ['source', 'target', 'weight']
    #    raw_links = _remove_duplicate(raw_links).fillna(0)
    gnetdata._add_netattr('links', raw_links)
    gnetdata._add_netattr_para('method', method)

    return gnetdata


def rundocker(gnetdata, method, cell_clusterid=None, Mms_TF=None, **kwargs):
    if method == 'GENIE3':
        gnetdata = _rundocker(gnetdata, method='GENIE3',
                              cell_clusterid=cell_clusterid, Mms_TF=Mms_TF)

    elif method == 'GRNBOOST2':
        gnetdata = _rundocker(gnetdata, method='GRNBOOST2',
                              cell_clusterid=cell_clusterid, Mms_TF=Mms_TF)

    elif method == 'PIDC':
        # remove genes with 0 counts
        gnetdata = _rundocker(gnetdata, method='PIDC',
                              cell_clusterid=cell_clusterid, Mms_TF=Mms_TF)

    elif method == 'SCNODE2VEC':
        gnetdata = _rundocker(gnetdata, method='SCNODE2VEC',
                              cell_clusterid=cell_clusterid, Mms_TF=Mms_TF, **kwargs)

    elif method == "CORR":
        gnetdata = _rundocker(gnetdata, method='CORR',
                              cell_clusterid=cell_clusterid, Mms_TF=Mms_TF)

    else:
        raise Exception("valid method: GENIE3, CORR, PIDC, GRNBOOST2, SCNODE2VEC")

    return gnetdata


def buildnet(gnetdata, threshold=None, top=None):
    if (top is None) & (threshold is not None):
        links_filter = gnetdata.NetAttrs['links'].loc[gnetdata.NetAttrs['links']['weight'] > threshold]
    elif (top is not None) & (threshold is None):
        links_filter = gnetdata.NetAttrs['links'].sort_values('weight', ascending=False).head(top)
    elif (top is None) & (threshold is None):
        links_filter = gnetdata.NetAttrs['links']
    else:
        raise Exception("Cannot filter by threshold and top!")
    G = nx.from_pandas_edgelist(links_filter,
                                source="source",
                                target="target",
                                edge_attr=True)
    gnetdata._add_netattr('graph', G)
    gnetdata._add_netattr_para('threshold', str(threshold))
    gnetdata._add_netattr_para('top', str(top))

    return gnetdata


if __name__ == '__main__':
    warnings.simplefilter("ignore")
