#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:25:56 2019

@author: mwu
"""
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import metrics

def mapping_edges(df_1, df_2, df_1_col_1, df_1_col_2, df_2_col_1, df_2_col_2):


    df_1['tmp1'] = df_1[[df_1_col_1, df_1_col_2]].apply(lambda x: '_'.join(x), axis = 1)
    df_1['tmp2'] = df_1[[df_1_col_2, df_1_col_1]].apply(lambda x: '_'.join(x), axis = 1)
    df_2['tmp1'] = df_2[[df_2_col_1, df_2_col_2]].apply(lambda x: '_'.join(x), axis = 1)

    return (len(set(df_1['tmp1']) & set(df_2['tmp1'])) + len(set(df_1['tmp2']) & set(df_2['tmp1'])))


def evaluation(links, Ref_links, Num_Genes):

#    links_filtered=links.loc[(abs(links['value']) > threshold) & (links['var1'] != links['var2'])]

    Detected = links.shape[0]
    TP = mapping_edges(links, Ref_links, 'source', 'target', 'node1', 'node2')
    FN =  Ref_links.shape[0] - TP
    FP = Detected - TP
    TN = (Num_Genes * (Num_Genes -1)/2) - Ref_links.shape[0] - Detected + TP

    Precision = TP/(TP + FP)
    Recall = TP/(TP + FN)
    FDR = FP/(TN + FP)
    print()
    F1_Score = (2*Precision*Recall)/(Precision + Recall)
    print('Detected:', Detected)
    print('TP:', TP, '\n', 'FN:', FN, '\n',
          'FP:',FP , '\n', 'TN:', TN )

    print('Precision:', Precision, '\n',  'Recall:',  Recall, '\n',
          'F1 Score:', F1_Score, '\n', 'FDR:', FP/(TN+FP) )
    return(Detected, TP, FN, FP, TN, Precision, Recall, FDR, F1_Score)


def test_run(links, Ref_links, input_dataset, filename = None):

        Detected, TP, FN, FP, TN, Precision, Recall, FDR, F1_Score = evaluation(links, Ref_links, input_dataset.shape[0])
        actual_ref = pd.DataFrame(0, columns = ['var1', 'var2', 'value'], index = links.index)
        actual_ref[['var1', 'var2']] = links[['source', 'target']]

        for i in range(0, len(Ref_links.index)):
#        print(Ref_links.iloc[i]['node1'], Ref_links.iloc[i]['node2'])
                actual_ref.loc[(actual_ref['var1'] == Ref_links.iloc[i]['node1']) & (actual_ref['var2'] == Ref_links.iloc[i]['node2']), 'value'] = 1
                actual_ref.loc[(actual_ref['var2'] == Ref_links.iloc[i]['node1']) & (actual_ref['var1'] == Ref_links.iloc[i]['node2']), 'value'] = 1

        actual_ref.loc[actual_ref['value'] < 1, 'value'] = 0

        auc = metrics.roc_auc_score(actual_ref['value'], links['weight'])

        fpr, tpr, thresholds = metrics.roc_curve(actual_ref['value'], links['weight'])

        return (fpr, tpr, auc)
    
#        plt.plot([0,1], [0,1], linestyle = '--')
#        plt.plot(fpr, tpr, marker = '.')
#        plt.title('ROC: AUC = {0:0.2f}'.format(auc))
#
#        plt.savefig(filename + '_ROC.png')
#        plt.close()

        #PR curve
#        average_precision = metrics.average_precision_score(actual_ref['value'], links['weight'])
#
#        precision, recall, threshold = metrics.precision_recall_curve(actual_ref['value'], links['weight'])
#        plt.step(recall, precision, color='b', alpha=0.2,
#         where='post')
#        plt.fill_between(recall, precision, step='post', alpha=0.2,
#                     color='b')
#
#        plt.xlabel('Recall')
#        plt.ylabel('Precision')
#        plt.ylim([0.0, 1.05])
#        plt.xlim([0.0, 1.0])
#        plt.title('Precision-Recall curve: AP={0:0.2f}'.format(average_precision))
#
#        plt.savefig(filename + '_PR.png')
#        plt.close()


