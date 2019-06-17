#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:51:34 2019

@author: mwu
"""
from __future__ import absolute_import
import sys
if('/home/mwu/MING_V9T/PhD_Pro/PySCNet/' not in sys.path): sys.path.append('/home/mwu/MING_V9T/PhD_Pro/PySCNet/')

import docker
import pandas as pd
import os
import networkx as nx
import tarfile
import warnings

global path
path = sys.path[-1] + 'BuildNet/Docker_App/'


def __copy_to(container_id, src, dst):
    client = docker.from_env()
    container = client.containers.get(container_id)
    strm, stat = container.get_archive(src)

    with open(os.getenv('HOME') + '/temp.tar', 'w') as f:
            for line in strm:
                    f.write(str(line, 'utf-8'))
            f.seek(0)

            thisTar = tarfile.TarFile(f.name)
            thisTar.extract('links.txt', path = dst)
    os.remove(os.getenv('HOME') + '/temp.tar')


def __remove_duplicate(links):
        links_list = sorted(links[['source', 'target']].values.tolist())
        for i in range(len(links_list)):
                links_list[i] = tuple(sorted(links_list[i]))
        nodes = pd.DataFrame(list(set(links_list)), columns=('source', 'target'))
        links = pd.merge(links, nodes, how='right')
        return(links)



def __rundocker(gnetdata,method, path = path):

        client = docker.from_env()
        pd.DataFrame.to_csv(gnetdata.GeneMatrix, path + method + '/Expr.txt',  sep= '\t')
        client.images.build(path = path + method, dockerfile = 'Dockerfile', tag = method.lower())
        container = client.containers.run(method.lower(), detach = True)
        __copy_to(container_id=container.short_id, src = '/' + method + '/links.txt', dst=os.getenv('HOME'))
#        client.remove_container(container.short_id)
        container.stop()
        client.containers.prune()
        client.images.prune()
        os.system('rm ' + path + method + '/Expr.txt')
        raw_links = pd.read_csv(os.getenv('HOME') + '/links.txt', sep = '\t', header = 0)
        raw_links.columns = ['source', 'target', 'weight']
        raw_links = __remove_duplicate(raw_links)
        gnetdata._add_netattr('links', raw_links)
        gnetdata._add_netattr_para('method', method)
        return(gnetdata)


def rundocker(gnetdata, method):


        if method == 'GENIE3':
                gnetdata = __rundocker(gnetdata, 'GENIE3')

        elif method == 'PIDC':
                gnetdata = __rundocker(gnetdata, 'PIDC')

        #TODO: check the output from SCODE
        elif method == "SCODE":
                gnetdata = __rundocker(gnetdata, 'SCODE')

        elif method == "CORR":
                gnetdata = __rundocker(gnetdata, 'CORR')
        #TODO: Input data with clusterid
        elif method == "SINCERA":
                gnetdata = __rundocker(gnetdata, 'SINCERA')
        #TODO: permission issue
        elif method == "SJARACNE":
                gnetdata = __rundocker(gnetdata, 'SJARACNE')

        else:
                raise Exception("valid method: GENIE3, PIDC, SCODE, CORR, SINCERA, SJARACNE")

        return(gnetdata)


def buildnet(gnetdata, threshold):

        links_filter = gnetdata.NetAttrs['links'].loc[gnetdata.NetAttrs['links']['weight'] > threshold]
        G = nx.from_pandas_edgelist(links_filter,
                                    source = "source",
                                    target= "target",
                                    edge_attr = True)
        gnetdata._add_netattr('graph', G)
        gnetdata._add_netattr_para('threshold', str(threshold))
        return(gnetdata)

if __name__ == '__main__':

        warnings.simplefilter("ignore")
