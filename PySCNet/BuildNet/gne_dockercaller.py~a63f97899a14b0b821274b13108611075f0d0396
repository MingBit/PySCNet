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



global path
#path = sys.path[-1] + 'BuildNet/Docker_App/'
path = os.path.join(os.path.dirname(__file__)) + '/Docker_App/'
print(path)

def _copy_to(container_id, src, dst):
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


def _remove_duplicate(links):
        links_list = sorted(links[['source', 'target']].values.tolist())
        for i in range(len(links_list)):
                links_list[i] = tuple(sorted(links_list[i]))
        nodes = pd.DataFrame(list(set(links_list)), columns=('source', 'target'))
        links = pd.merge(links, nodes, how='right')
        return(links)



def _rundocker(gnetdata, method, path = path, time_point = None, cell_clusterid = None):

        client = docker.from_env()
        pd.DataFrame.to_csv(gnetdata.ExpMatrix, path + method + '/Expr.txt',  sep= '\t')
        
        if time_point is not None:
            pd.DataFrame.to_csv(time_point, path + method + '/time_point.txt',  sep= '\t', header=False, index=False)
        
        if cell_clusterid is not None:
            pd.DataFrame.to_csv(cell_clusterid, path + method + '/cell_clusterid.txt',  sep= '\t', header=False, index=False)
            
        
        client.images.build(path = path + method, dockerfile = 'Dockerfile', tag = method.lower())
        container = client.containers.run(method.lower(), detach = True)
        _copy_to(container_id=container.short_id, src = '/' + method + '/links.txt', dst=os.getenv('HOME'))

#        client.remove_container(container.short_id)
        container.stop()
        client.containers.prune()
        client.images.prune()
        os.system('rm ' + path + method + '/Expr.txt')
        raw_links = pd.read_csv(os.getenv('HOME') + '/links.txt', sep = '\t', header = 0)
        raw_links.columns = ['source', 'target', 'weight']
        raw_links = _remove_duplicate(raw_links)
        gnetdata._add_netattr('links', raw_links)
        gnetdata._add_netattr_para('method', method)

        return(gnetdata)


def rundocker(gnetdata, method, time_point = None, cell_clusterid = None):


        if method == 'GENIE3':
                gnetdata = _rundocker(gnetdata, 'GENIE3')

        elif method == 'PIDC':
                gnetdata = _rundocker(gnetdata, 'PIDC')

        elif method == 'JUMP3':
                gnetdata = _rundocker(gnetdata, 'JUMP3')

        elif method == "SCODE":
                gnetdata = _rundocker(gnetdata, 'SCODE', time_point = time_point)

        elif method == "CORR":
                gnetdata = _rundocker(gnetdata, 'CORR')
        #TODO: Input data with clusterid
        elif method == "SINCERA":
                gnetdata = _rundocker(gnetdata, 'SINCERA')
        #TODO: permission issue
        elif method == "SJARACNE":
                gnetdata = _rundocker(gnetdata, 'SJARACNE')

        else:
                raise Exception("valid method: GENIE3, PIDC, SCODE, CORR, SINCERA, SJARACNE, JUMP3")

        return(gnetdata)


def buildnet(gnetdata, threshold = None, top = None):


        if((top is None) & (threshold is not None)):
                links_filter = gnetdata.NetAttrs['links'].loc[gnetdata.NetAttrs['links']['weight'] > threshold]
        elif((top is not None) & (threshold is None)):
                links_filter = gnetdata.NetAttrs['links'].sort_values('weight', ascending = False).head(top)
        elif((top is None) & (threshold is None)):
                links_filter = gnetdata.NetAttrs['links']
        else:
                raise Exception("Cannot filter by threshold and top!")
        G = nx.from_pandas_edgelist(links_filter,
                                    source = "source",
                                    target= "target",
                                    edge_attr = True)
        gnetdata._add_netattr('graph', G)
        gnetdata._add_netattr_para('threshold', str(threshold))
        gnetdata._add_netattr_para('top', str(top))

        return(gnetdata)

if __name__ == '__main__':

        warnings.simplefilter("ignore")
