#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 12:23:37 2019

@author: mwu
"""

def Genes_modules_Cell_Clusters(gnetdata):
        """it returns summary of genes in each modules"""

        node_group = gnetdata.NetAttrs['communities']
        gene_module = node_group.groupby('group').count()

        return(gene_module)