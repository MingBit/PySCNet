#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:51:34 2019

@author: mwu
"""
from __future__ import absolute_import
import os
import sys

global path
path = os.path.join(os.path.dirname(__file__)) + '/Docker_App/'
sys.path.append("..")

import docker
import shutil
import pandas as pd
import tarfile
import warnings
import _pickle as pk
from ..utils import *

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


def __rundocker(gnetdata, method, cell, feature, **kwargs):

    client = docker.from_env()
    cell_clusterid = kwargs.get('cell_clusterid', None)
    select_by = kwargs.get('select_by', None)
    Mms_TF = kwargs.get('Mms_TF', None)
    directed = kwargs.get('directed', True)
    
    shutil.copyfile(os.path.join('/', *(os.path.dirname(__file__).split(os.path.sep))[:-1], 'utils.py'), 
                    path + method + '/utils.py')
    
    tmp_expr = gnetdata_subset(gnetdata, cell, feature, **kwargs)

    if Mms_TF is not None:
        pd.DataFrame.to_csv(pd.DataFrame(Mms_TF), path + method + '/TF_Names.txt', sep='\t', header=False, index=False)

    with open(path + method + '/paras.pk', 'wb') as outfile:
        pk.dump(kwargs, outfile)
        
    pd.DataFrame.to_csv(tmp_expr, path + method + '/Expr.txt', sep='\t')
    client.images.build(path=path + method, dockerfile='Dockerfile', tag=method.lower())
    container = client.containers.run(method.lower(), detach=True)
    
    __copy_to(container_id=container.short_id, src='/' + method + '/links.txt', dst=os.getenv('HOME'))

    #client.remove_container(container.short_id)
    container.stop()
    client.containers.prune()
    client.images.prune()
    raw_links = pd.read_csv(os.getenv('HOME') + '/links.txt', sep='\t', header=None)
    raw_links.columns = ['source', 'target', 'weight']
    os.system('rm ' + path + method + '/*.txt | rm ' + path + method + '/paras.pk')

    if directed is False: raw_links = order_source_target(raw_links).drop_duplicates(subset=['source', 'target'])

    gnetdata._add_netattr(method + '_links', raw_links)
    print(method + '_links added into NetAttrs')
    return gnetdata


def rundocker(gnetdata, method, cell=None, feature=None, **kwargs):
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
    :param directed: bool, default True.
    :param kwargs: additional parameters passed to scnode2vec()
    :return: Gnetdata object with links saved in NetAttrs
    """

    valid_methods = ['GENIE3', 'CORR', 'PIDC', 'GRNBOOST2', 'SCNODE2VEC']
    assert method in valid_methods, f"valid method: {', '.join(valid_methods)}"

    gnetdata = __rundocker(gnetdata, method=method, cell = cell, feature = feature, **kwargs)
    return gnetdata


if __name__ == '__main__':
    warnings.simplefilter("ignore")
