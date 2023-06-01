#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:51:34 2019

@author: mwu
"""

import pandas as pd


def gnetdata_subset(gnetdata, cell, feature, **kwargs):
    
    cell_clusterid = kwargs.get('cell_clusterid', None)
    select_by = kwargs.get('select_by', None)
    
    feature = gnetdata.Exp['feature'] if feature is None else feature

    if cell_clusterid is None:
        cell = gnetdata.Exp['cell'] if cell is None else cell
    else:
        cell_info = gnetdata.CellAttrs['CellInfo']
        # cell = list(cell_info.loc[cell_info[select_by].isin([cell_clusterid])].index)
        cell = cell_info.loc[cell_info[select_by].eq(cell_clusterid)].index.tolist()


    subExpr = pd.DataFrame(gnetdata.Exp['matrix'],
                            index = gnetdata.Exp['cell'],
                            columns = gnetdata.Exp['feature']).loc[cell, feature].T
    
    return subExpr


def order_source_target(df):
    
    df_new = df[['source', 'target']].apply(sorted, axis=1, result_type='expand').rename(columns={0: 'source', 1: 'target'})
    df_new['weight'] = df['weight']

    return df_new


