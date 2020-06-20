#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 18:58:48 2019

@author: mwu
"""

import pymc3 as pm
import pandas as pd
import numpy as np
import seaborn as sns
from warnings import filterwarnings
from pymc3 import floatX

filterwarnings('ignore')
sns.set_style('white')
from sklearn.model_selection import train_test_split
import seaborn as sns
from pymc3.theanof import set_tt_rng, MRG_RandomStreams

set_tt_rng(MRG_RandomStreams(42))
from theano.compile.sharedvalue import shared


def _construct_bnn(ann_input, ann_output, input_nodes, hidden_layers, total_size):
    # initialize random weights between layers
    init_1 = floatX(np.random.randn(input_nodes, hidden_layers[0]))
    init_2 = floatX(np.random.randn(hidden_layers[0], hidden_layers[1]))
    init_out = floatX(np.random.randn(hidden_layers[1]))

    with pm.Model() as ann_model:
        # prior distribution for weights
        weights_in_1 = pm.Normal('w_in_1', 0, sd=1, testval=init_1, shape=(input_nodes, hidden_layers[0]))
        weights_1_2 = pm.Normal('w_1_2', 0, sd=1, testval=init_2, shape=(hidden_layers[0], hidden_layers[1]))
        weights_2_out = pm.Normal('w_2_out', 0, sd=1, testval=init_out, shape=(hidden_layers[1],))

        # prior distribution for bias
        bias_in_1 = pm.Normal('bias_in_1', 0, sd=1)
        bias_1_2 = pm.Normal('bias_1_2', 0, sd=1)
        bias_2_out = pm.Normal('bias_2_out', 0, sd=1)

        #        pdb.set_trace()
        # non-linear transformation
        act_1 = pm.math.tanh(pm.math.dot(ann_input, weights_in_1) + bias_in_1)
        act_2 = pm.math.tanh(pm.math.dot(act_1, weights_1_2) + bias_1_2)
        act_out = pm.math.sigmoid(pm.math.dot(act_2, weights_2_out) + bias_2_out)

        out = pm.Bernoulli('out', act_out, observed=ann_output, total_size=total_size)

    return ann_model


def bnn_classifier(X, Y, test_size=0.4, hideen_layers=[50, 20]):
    x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=test_size)

    ann_input = shared(np.asarray(x_train.iloc[:, 2:]))
    ann_output = shared(np.asarray(y_train))
    neural_network = _construct_bnn(ann_input, ann_output, input_nodes=x_train.shape[1] - 2,
                                    hidden_layers=hideen_layers, total_size=y_train.shape[0])
    with neural_network:
        approx = pm.fit(method=pm.ADVI())
        trace = approx.sample(draws=5000)
        ppc_train = pm.sample_posterior_predictive(trace, samples=500, model=neural_network)
        #                pred_train = pd.Series((ppc_train['out'].mean(0) > 0.5) * 1, index = y_train.index)
        pred_train = pd.Series(ppc_train['out'].mean(0), index=y_train.index)

    ann_input.set_value(x_test.iloc[:, 2:])
    ann_output.set_value(y_test)

    with neural_network:
        ppc_test = pm.sample_posterior_predictive(trace, samples=500, model=neural_network)
        pred_test = pd.Series(ppc_test['out'].mean(0), index=y_test.index)
    #                pred_test = pd.Series((ppc_test['out'].mean(0) > 0.5) * 1, index = y_test.index)

    res_df = pd.DataFrame(
        {'source': x_train.source.append(x_test.source), 'target': x_train.target.append(x_test.target),
         'weight': list(pred_train) + list(pred_test)})

    return (res_df)
