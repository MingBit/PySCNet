#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:51:34 2019

@author: mwu
"""

import docker
import pandas as pd


path = '/home/mwu/MING_V9T/PhD_Pro/PySCNet/BuildNet/Docker_App/'

def buildnet(gedata, method = 'GENIE3', paras = None):

        pd.DataFrame.to_csv(gedata.GeneMatrix, path + 'GENIE3/Expr.csv')

        if method == 'GENIE3':
                client = docker.from_env()
                client.images.build(path + 'GENIE3/', 'genie3')
                client.containers.run('genie3', detach = True)

        elif method == 'PIDC':
                clinet = docker.from_env()
                client.images.build(path + 'PIDC/', 'pidc')
                clinet.containers.run('pidc', detach = True)

        elif method == "SCNS":
                client = docker.from_env()
                client.images.build(path + 'SCNS/', 'scns')
                client.containers.run('scns', detach = True)

        elif method == "SCENIC":
                client = docker.from_env()
                client.images.build(path + 'SCENIC/', 'scenic')
                client.containers.run('scenic', detach = True)
