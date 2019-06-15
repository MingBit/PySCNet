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


def __rundocker(gedata,method):

        client = docker.from_env()
        pd.DataFrame.to_csv(gedata.GeneMatrix, path + method + '/Expr.txt')
        client.images.build(path = path + method, dockerfile = 'Dockerfile', tag = method.lower())
        container = client.containers.run(method.lower(), detach = True)
        __copy_to(container_id=container.short_id, src = '/' + method + '/links.txt', dst=os.getenv('HOME'))
#        client.remove_container(container.short_id)
        container.stop()
        client.containers.prune()
        client.images.prune()
        raw_links = pd.read_csv(os.getenv('HOME') + '/links.txt', sep = '\t', header = 0)
        raw_links.columns = ['source', 'target', 'weight']
        gedata._add_netattr('links', raw_links)

        return(gedata)


def rundocker(gedata, method):


        if method == 'GENIE3':
                gedata = __rundocker(gedata, 'GENIE3')

        elif method == 'PIDC':
                gedata = __rundocker(gedata, 'PIDC')

        elif method == "SCNS":
                gedata = __rundocker(gedata, 'SCNS')

        elif method == "CORR":
                gedata = __rundocker(gedata, 'CORR')

        elif method == "SINCERA":
                gedata = __rundocker(gedata, 'SINCERA')

        elif method == "SJARACNE":
                gedata = __rundocker(gedata, 'SJARACNE')

        else:
                print("invalid method!")

        return(gedata)


def buildnet(gedata, threshold):

        links_filter = gedata.NetAttrs['links'].loc[gedata.NetAttrs['links']['weight'] > threshold]
        G = nx.from_pandas_edgelist(links_filter,
                                    source = "source",
                                    target= "target",
                                    edge_attr = True)
        gedata._add_netattr('graph', G)

        return(gedata)

