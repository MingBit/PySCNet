#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:25:56 2019

@author: mwu
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn import metrics
import pdb



def mapping_edges(df_1, df_2, df_1_col_1, df_1_col_2, df_2_col_1, df_2_col_2):

    
    df_1['tmp1'] = df_1[[df_1_col_1, df_1_col_2]].apply(lambda x: '_'.join(x), axis = 1)
    # df_1['tmp2'] = df_1[[df_1_col_2, df_1_col_1]].apply(lambda x: '_'.join(x), axis = 1)
    df_2['tmp1'] = df_2[[df_2_col_1, df_2_col_2]].apply(lambda x: '_'.join(x), axis = 1)

  
    return (len(set(df_1['tmp1']) & set(df_2['tmp1'])))


def evaluation(links, Ref_links, Num_Genes):

#    links_filtered=links.loc[(abs(links['value']) > threshold) & (links['var1'] != links['var2'])]

    Detected = links.shape[0]
    Ref_links = Ref_links[Ref_links.weight != 0]
    TP = mapping_edges(links, Ref_links, 'source', 'target', 'source', 'target')
    FN =  Ref_links.shape[0] - TP
    FP = Detected - TP
    TN = (Num_Genes * Num_Genes) - Ref_links.shape[0] - Detected + TP

    Precision = TP/(TP + FP)
    Recall = TP/(TP + FN)
    FDR = FP/(TN + FP)
    
    F1_Score = (2*Precision*Recall)/(Precision + Recall)
    print('Detected:', Detected)
    print('TP:', TP, '\n', 'FN:', FN, '\n', 'FP:', FP , '\n', 'TN:', TN )

    print('Precision:', Precision, '\n',  'Recall:',  Recall, '\n',
          'F1 Score:', F1_Score, '\n', 'FDR:', FP/(TN+FP) )
    return(Detected, TP, FN, FP, TN, Precision, Recall, FDR, F1_Score)


def test_run(links, Ref_links, input_dataset, filename = None):

    Detected, TP, FN, FP, TN, Precision, Recall, FDR, F1_Score = evaluation(links, Ref_links, input_dataset.shape[0])
    
    Comp_Links = pd.merge(links, Ref_links, on = ['source', 'target'], how = 'right').fillna(0)

    
    auc = metrics.roc_auc_score(np.array(Comp_Links['weight_y'].abs()), np.array(Comp_Links['weight_x'].abs()))
    fpr, tpr, threshold_1 = metrics.roc_curve(Comp_Links['weight_y'].abs(), Comp_Links['weight_x'].abs())
    pre, recall, threshold_2 = metrics.precision_recall_curve(Comp_Links['weight_y'].abs(), Comp_Links['weight_x'].abs())
    avg_pre = metrics.average_precision_score(Comp_Links['weight_y'].abs(), Comp_Links['weight_x'].abs())
        

    return ([fpr, tpr, pre, recall, auc, avg_pre, Precision, Recall, F1_Score])
    



