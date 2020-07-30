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
import tarfile
import warnings
import _pickle as pk

global path
path = os.path.join(os.path.dirname(__file__)) + '/Docker_App/'


def __copy_to(container_id, src, dst):
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


def __remove_duplicate(links):
    links_list = sorted(links[['source', 'target']].values.tolist())
    for i in range(len(links_list)):
        links_list[i] = tuple(sorted(links_list[i]))
    nodes = pd.DataFrame(list(set(links_list)), columns=('source', 'target'))
    links = pd.merge(links, nodes, how='right')

    return links


def __rundocker(gnetdata, method, cell=None, feature=None, cell_clusterid=None, select_by=None, Mms_TF=None, **kwargs):
    client = docker.from_env()

    if feature is None:
        feature = gnetdata.ExpMatrix.index

    if cell_clusterid is None:
        cell = gnetdata.ExpMatrix.columns if cell is None else cell
    else:
        cell_info = gnetdata.CellAttrs['CellInfo']
        cell = list(cell_info.loc[cell_info[select_by].isin([cell_clusterid])].index)

    tmp_expr = gnetdata.ExpMatrix.loc[feature, cell]

    if Mms_TF is not None:
        pd.DataFrame.to_csv(Mms_TF, path + method + '/TF_Names.txt', sep='\t', header=False, index=False)

    with open(path + method + '/paras.pk', 'wb') as outfile:
        pk.dump(kwargs, outfile)

    pd.DataFrame.to_csv(tmp_expr, path + method + '/Expr.txt', sep='\t')
    client.images.build(path=path + method, dockerfile='Dockerfile', tag=method.lower())
    container = client.containers.run(method.lower(), detach=True)
    __copy_to(container_id=container.short_id, src='/' + method + '/links.txt', dst=os.getenv('HOME'))

    #        client.remove_container(container.short_id)
    container.stop()
    client.containers.prune()
    client.images.prune()
    os.system('rm ' + path + method + '/Expr.txt | rm ' + path + method + '/paras.pk')
    raw_links = pd.read_csv(os.getenv('HOME') + '/links.txt', sep='\t', header=None)
    raw_links.columns = ['source', 'target', 'weight']
    gnetdata._add_netattr(method + '_links', raw_links)
    print(method + '_links added into NetAttrs')
    return gnetdata


def rundocker(gnetdata, method, cell=None, feature=None, cell_clusterid=None, select_by=None, Mms_TF=None, **kwargs):
    """
    Call GRN methods via docker.
    -------------------------------------
    :param gnetdata: Gnetdata object, default None.
    :param method: str, default None. methods: [GENIE3, GRNBOOST2, PIDC, CORR]
    :param cell: list, default None. a list of cell names
    :param feature: list, default None. a list of gene names
    :param cell_clusterid: str, default None. cell with cell_clusterid will be selected
    :param select_by: str, default None. key of filtering cells
    :param Mms_TF: list, default None. a list of transcription factor names
    :param kwargs: additional parameters passed to scnode2vec()
    :return: Gnetdata object with links saved in NetAttrs
    """
    if method == 'GENIE3':
        gnetdata = __rundocker(gnetdata, method='GENIE3', cell=cell, feature=feature,
                               cell_clusterid=cell_clusterid, select_by=select_by, Mms_TF=Mms_TF)

    elif method == 'GRNBOOST2':
        gnetdata = __rundocker(gnetdata, method='GRNBOOST2', cell=cell, feature=feature,
                               cell_clusterid=cell_clusterid, select_by=select_by, Mms_TF=Mms_TF)

    elif method == 'PIDC':
        # remove genes with 0 counts
        gnetdata = __rundocker(gnetdata, method='PIDC', cell=cell, feature=feature,
                               cell_clusterid=cell_clusterid, select_by=select_by, Mms_TF=Mms_TF)

    elif method == 'SCNODE2VEC':
        gnetdata = __rundocker(gnetdata, method='SCNODE2VEC', cell=cell, feature=feature,
                               cell_clusterid=cell_clusterid, select_by=select_by, Mms_TF=Mms_TF, **kwargs)

    elif method == "CORR":
        gnetdata = __rundocker(gnetdata, method='CORR', cell=cell, feature=feature,
                               cell_clusterid=cell_clusterid, select_by=select_by, Mms_TF=Mms_TF)

    else:
        raise Exception("valid method: GENIE3, CORR, PIDC, GRNBOOST2, SCNODE2VEC")

    return gnetdata


if __name__ == '__main__':
    warnings.simplefilter("ignore")
