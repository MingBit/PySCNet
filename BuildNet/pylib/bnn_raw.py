#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 12:41:21 2019

@author: angela
"""

import pymc3 as pm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from warnings import filterwarnings
from pymc3 import floatX
filterwarnings('ignore')
sns.set_style('white')
from sklearn.model_selection import train_test_split
from sklearn import metrics
import seaborn as sns



def Generate_Zero(sim_df):

    for cell in sim_df.index:
        tmp = sim_df.loc[cell, (sim_df.loc[cell] < sim_df.mean(axis = 0))]
        tmp_np = np.random.binomial(1, 0.5, len(tmp))
        for i in range(0, len(tmp)):
            tmp[i] = tmp[i] if tmp_np[i] == 1 else 0
        sim_df.loc[cell, (sim_df.loc[cell] < sim_df.mean(axis = 0))] = tmp

    return (sim_df)


def binary_converter(m):
    if type(m) == list:
        ml = []
        for l in m:
            ml.append(((l > l.mean(0))*1).astype(int))
        return ml
    else:
        return ((m > m.mean(0))*1).astype(int)

def x_y_generator(Expr_df, gene, cell = None, test_size = .6):

    if(cell is None):
        return train_test_split(Expr_df.drop(columns = gene), Expr_df[gene], test_size = test_size)

    elif(gene is None):
        return train_test_split(Expr_df.drop(index = cell).transpose(), Expr_df.loc[cell], test_size = test_size)

    else:
        print('check if gene or cell is NONE!')


def construct_bnn(ann_input, ann_output, input_nodes, hidden_layers, total_size):


    #initialize random weights between layers
    init_1 = floatX(np.random.randn(input_nodes, hidden_layers[0]))
    init_2 = floatX(np.random.randn(hidden_layers[0], hidden_layers[1]))
    init_out = floatX(np.random.randn(hidden_layers[1]))

    with pm.Model() as ann_model:

        #prior distribution for weights
        weights_in_1 = pm.Normal('w_in_1', 0, sd=1, testval = init_1, shape = (input_nodes, hidden_layers[0]))
        weights_1_2 = pm.Normal('w_1_2', 0, sd=1, testval = init_2, shape = (hidden_layers[0], hidden_layers[1]))
        weights_2_out = pm.Normal('w_2_out', 0, sd = 1, testval = init_out, shape = (hidden_layers[1],))

        #prior distribution for bias
        bias_in_1 = pm.Normal('bias_in_1', 0, sd=1)
        bias_1_2 = pm.Normal('bias_1_2', 0, sd=1)
        bias_2_out = pm.Normal('bias_2_out', 0, sd=1)

#        pdb.set_trace()
        #non-linear transformation
        act_1 = pm.math.tanh(pm.math.dot(ann_input, weights_in_1) + bias_in_1)
        act_2 = pm.math.tanh(pm.math.dot(act_1, weights_1_2) + bias_1_2)
        act_out = pm.math.sigmoid(pm.math.dot(act_2, weights_2_out) + bias_2_out)

        sd = pm.HalfNormal('sd', sd = 1)

#likelihood
#        out = pm.distributions.discrete.ZeroInflatedNegativeBinomial('out', psi = .2,
#                                                               alpha = 5, mu = act_out, observed = ann_output)
#        out = pm.Bernoulli('out', logit_p = act_out, observed = ann_output)
        out = pm.Normal('out', mu=act_out, sd = sd, observed = ann_output, total_size = total_size)

    return ann_model


#fit and evaluate
from theano.compile.sharedvalue import shared
from pymc3.theanof import set_tt_rng, MRG_RandomStreams
set_tt_rng(MRG_RandomStreams(42))
import time

def pred_eva(x_train, x_test, y_train, y_test, bnn_func, input_nodes, hidden_layers, bnn_kwargs = None, sample_kwargs = None):

    if bnn_kwargs is None:
        bnn_kwargs = {}

    if sample_kwargs is None:
        sample_kwargs = {}

    ann_input = shared(np.asarray(x_train))
    ann_output = shared(np.asarray(y_train))
#    ann_input = pm.Minibatch(x_train, batch_size=50)
#    ann_output = pm.Minibatch(y_train, batch_size=50)
    model = bnn_func(ann_input, ann_output, input_nodes, hidden_layers, y_train.shape[0], **bnn_kwargs)
#    print(model['out'].tag.test_value)

    start_time = time.clock()

    with model:

        minibatches = {
            ann_input: pm.Minibatch(x_train, batch_size=10),
            ann_output: pm.Minibatch(y_train, batch_size=10)
        }

        inference_args = {
            'n' : 10000,
            'callbacks': [pm.callbacks.CheckParametersConvergence()],
            'obj_optimizer': pm.adagrad_window(learning_rate=2e-4),
             'more_replacements': minibatches
             }

        #fit model
        approx = pm.fit(method=pm.ADVI(), **inference_args)
        trace = approx.sample(draws = 3000, **sample_kwargs)


    with model:

        ppc_train = pm.sample_posterior_predictive(trace, samples=500, model = model)
        #    pred_train = pd.Series((ppc_train['out'].mean(axis = 0) > 0.5) * 1, index = y_train.index)

        pred_train = pd.Series(ppc_train['out'].mean(0), index = y_train.index)

#    make prediction
    ann_input.set_value(x_test)
    ann_output.set_value(y_test)

    with model:
        ppc_test = pm.sample_posterior_predictive(trace, samples=500, model = model)

        #pred_test = pd.Series((ppc_test['out'].mean(0) > 0.5) * 1, index = y_test.index)
        pred_test = pd.Series(ppc_test['out'].mean(0), index = y_test.index)

    print(time.clock() - start_time, "seconds")
    return pred_train, pred_test




def mapping_edges(df_1, df_2, df_1_col_1, df_1_col_2, df_2_col_1, df_2_col_2):


    df_1['tmp1'] = df_1[[df_1_col_1, df_1_col_2]].apply(lambda x: '_'.join(x), axis = 1)
    df_1['tmp2'] = df_1[[df_1_col_2, df_1_col_1]].apply(lambda x: '_'.join(x), axis = 1)
    df_2['tmp1'] = df_2[[df_2_col_1, df_2_col_2]].apply(lambda x: '_'.join(x), axis = 1)

    return (len(set(df_1['tmp1']) & set(df_2['tmp1'])) + len(set(df_1['tmp2']) & set(df_2['tmp1'])))


def evaluation(links, Ref_links, threshold, Num_Genes):

    links_filtered=links.loc[(abs(links['value']) > threshold) & (links['var1'] != links['var2'])]

    Detected = links_filtered.shape[0]
    TP = mapping_edges(links_filtered, Ref_links, 'var1', 'var2', 'node1', 'node2')
    FN =  Ref_links.shape[0] - TP
    FP = Detected - TP
    TN = (Num_Genes * (Num_Genes -1)/2) - Ref_links.shape[0] - Detected + TP

    Precision = TP/(TP + FP)
    Recall = TP/(TP + FN)
    FDR = FP/(TN + FP)

    F1_Score = (2*Precision*Recall)/(Precision + Recall)

    return(Detected, TP, FN, FP, TN, Precision, Recall, FDR, F1_Score)


# =============================================================================
# Iteratively update gene matrix

def get_BNN_table(input_dataset, output_path, filename, hidden_layers):

    for gene in input_dataset.columns:

         X_train, X_test, Y_train, Y_test = x_y_generator(input_dataset, gene)
         gene_train, gene_test = pred_eva(x_train=X_train, x_test=X_test, y_train=Y_train, y_test=Y_test,
                                          hidden_layers = hidden_layers, bnn_func=construct_bnn, input_nodes=X_train.shape[1])
         
         input_dataset[gene] = (gene_test.append(gene_train)).reindex(index = input_dataset.index)    
    
    input_dataset.to_csv(output_path + filename + '.csv')

    return(input_dataset)
# # =============================================================================
# =============================================================================
# create new Genematrix, gene updated based on raw values of other genes' expressions
# =============================================================================
#     if(BNN_Expr_data.index) == 0:
#            BNN_Expr_data[gene] = gene_test.append(gene_train)
#     else:
#            BNN_Expr_data[gene] = (gene_test.append(gene_train)).reindex(index = BNN_Expr_data.index)
#
#    BNN_Expr_data.to_csv(path + output_path + filename + '.csv')

# create Matrix_Feature_1 and Matrix _Feature_2
#def get_BNN_table(input_dataset, output_path, filename, hidden_layers):
#
#        feature_1_Matrix = input_dataset.copy()
#        feature_2_Matrix = input_dataset.copy()
#
#        for gene in feature_1_Matrix.columns:
#                X_train, X_test, Y_train, Y_test = x_y_generator(feature_1_Matrix, gene)
#                gene_train, gene_test = pred_eva(x_train=X_train, x_test=X_test, y_train=Y_train, y_test=Y_test,
#                                         hidden_layers = hidden_layers, bnn_func=construct_bnn, input_nodes=X_train.shape[1])
#
#                feature_1_Matrix[gene] = (gene_test.append(gene_train)).reindex(index = feature_1_Matrix.index)
#        feature_1_Matrix.to_csv(path + output_path + filename + '_Feature_1.csv')
#
#        for cell in feature_2_Matrix.index:
#                X_train, X_test, Y_train, Y_test = x_y_generator(feature_2_Matrix, gene = None, cell = cell)
#                gene_train, gene_test = pred_eva(x_train=X_train, x_test=X_test, y_train=Y_train,
#                                         y_test=Y_test, hidden_layers = hidden_layers, bnn_func=construct_bnn,
#                                         input_nodes=X_train.shape[1])
#
#                feature_2_Matrix.loc[cell] = (gene_test.append(gene_train)).reindex(index =
#                            feature_2_Matrix.transpose().index)
#
#        feature_2_Matrix.to_csv(path + output_path + filename + '_Feature_2.csv')
#
#        Features_Avg = pd.concat([feature_1_Matrix, feature_2_Matrix]).groupby(level = 0).mean()
#        Features_Avg.to_csv(path + output_path + filename + '_Feature_Avg.csv')
#
#        return(Features_Avg)


def test_run(BNN_Expr_data, Ref_links, threshold, input_dataset, output_path, filename):

#    from sklearn.utils import shuffle
    links = BNN_Expr_data.corr(method = 'spearman').unstack().sort_values().drop_duplicates().reset_index()
    links.columns = ['var1', 'var2','value']
    links = links.loc[(links['var1'] != links['var2'])]

#    links = links.reindex(shuffle(np.array(links.index)))

    Detected, TP, FN, FP, TN, Precision, Recall, FDR, F1_Score = evaluation(links, Ref_links, threshold, input_dataset.shape[1])
    with open(output_path + filename + 'Reslog.txt', 'w') as text_file:
        text_file.write('Rlog\n Threshold :%s\n Detected: %s\n TP: %s\n FP: %s\n FN: %s\n TN:%s\n Precision: %s\n Recall: %s\n F1 Score: %s\n FDR: %s\n'
                     % (threshold, Detected, TP, FP, FN, TN, Precision, Recall, F1_Score, FDR))

    print('Detected:', Detected)
    print('TP:', TP, '\n', 'FN:', FN, '\n',
          'FP:',FP , '\n', 'TN:', TN )

    print('Precision:', Precision, '\n',  'Recall:',  Recall, '\n',
          'F1 Score:', F1_Score, '\n', 'FDR:', FP/(TN+FP) )

    #Performance AUC and PR

    actual_ref = pd.DataFrame(0, columns = ['var1', 'var2', 'value'], index = links.index)
    actual_ref[['var1', 'var2']] = links[['var1', 'var2']]

    for i in range(0, len(Ref_links.index)):
#        print(Ref_links.iloc[i]['node1'], Ref_links.iloc[i]['node2'])
        actual_ref.loc[(actual_ref['var1'] == Ref_links.iloc[i]['node1']) & (actual_ref['var2'] == Ref_links.iloc[i]['node2']), 'value'] = 1
        actual_ref.loc[(actual_ref['var2'] == Ref_links.iloc[i]['node1']) & (actual_ref['var1'] == Ref_links.iloc[i]['node2']), 'value'] = 1

    actual_ref.loc[actual_ref['value'] < 1, 'value'] = 0
#    pdb.set_trace()
    auc = metrics.roc_auc_score(actual_ref['value'], links['value'])
#    print('AUC: %.3f' % auc)
    #ROC curve
    fpr, tpr, thresholds = metrics.roc_curve(actual_ref['value'], links['value'])
    plt.plot([0,1], [0,1], linestyle = '--')
    plt.plot(fpr, tpr, marker = '.')
    plt.title('ROC: AUC = {0:0.2f}'.format(auc))

    plt.savefig(output_path + filename + '_ROC.png')
    plt.close()

    #PR curve
    average_precision = metrics.average_precision_score(actual_ref['value'], links['value'])

    precision, recall, threshold = metrics.precision_recall_curve(actual_ref['value'], links['value'])
    plt.step(recall, precision, color='b', alpha=0.2,
         where='post')
    plt.fill_between(recall, precision, step='post', alpha=0.2,
                     color='b')

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('Precision-Recall curve: AP={0:0.2f}'.format(average_precision))

    plt.savefig(output_path + filename + '_PR.png')
    plt.close()




