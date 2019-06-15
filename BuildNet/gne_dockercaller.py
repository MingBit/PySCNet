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

def copy_to(container_id, src, dst):
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


def buildnet(gedata, threshold):

        links_filter = gedata.NetAttrs['links'].loc[gedata.NetAttrs['links']['weight'] > threshold]
        G = nx.from_pandas_edgelist(links_filter,
                                            source = "source",
                                            target= "target",
                                            edge_attr = True)
        gedata._add_netattr('graph', G)

        return(gedata)

def rundocker(gedata, method):


        if method == 'GENIE3':
                client = docker.from_env()
                pd.DataFrame.to_csv(gedata.GeneMatrix, path + 'GENIE3/Expr.txt')
                client.images.build(path = path + "/GENIE3", dockerfile = 'Dockerfile', tag = 'genie3')
                container = client.containers.run('genie3', detach = True)
                copy_to(container_id=container.short_id, src = '/GENIE3/links.txt', dst=os.getenv('HOME'))

        elif method == 'PIDC':
                clinet = docker.from_env()
                pd.DataFrame.to_csv(gedata.GeneMatrix, path + 'PIDC/Expr.csv')
                client.images.build(path = path + 'PIDC/', dockerfile = 'Dockerfile', tag = 'pidc')
                clinet.containers.run('pidc', detach = True)
                copy_to(container_id=container.short_id, src = '/PIDC/links.txt', dst=os.getenv('HOME'))

        elif method == "SCNS":
                client = docker.from_env()
                pd.DataFrame.to_csv(gedata.GeneMatrix, path + 'SCNS/Expr.csv')
                client.images.build(path + 'SCNS/', 'scns')
                client.containers.run('scns', detach = True)
                copy_to(client, 'scns:/links.txt', '/links.txt')

        elif method == "CORR":
                client = docker.from_env()
                pd.DataFrame.to_csv(gedata.GeneMatrix, path + 'CORR/Expr.csv')
                client.images.build(path + 'CORR/', 'corr')
                client.containers.run('corr', detach = True)
                copy_to(client, 'corr:/links.txt', '/links.txt')

        elif method == "SINCERA":
                client = docker.from_env()
                pd.DataFrame.to_csv(gedata.GeneMatrix, path + 'SINCERA/Expr.csv')
                client.images.build(path + 'SINCERA/', 'sincera')
                client.containers.run('sincera', detach = True)
                copy_to(client, 'sincera:/links.txt', '/links.txt')

        elif method == "SJARACNE":
                client = docker.from_env()
                pd.DataFrame.to_csv(gedata.GeneMatrix, path + 'SJARACNE/Expr.csv')
                client.images.build(path + 'SJARACNE/', 'sjaracne')
                client.containers.run('sjaracne', detach = True)
                copy_to(client, 'sjaracne:/links.txt', '/links.txt')

        raw_links = pd.read_csv(os.getenv('HOME') + '/links.txt', sep = '\t', header = 0)
        raw_links.columns = ['source', 'target', 'weight']
        gedata._add_netattr('links', raw_links)


        return(gedata)



