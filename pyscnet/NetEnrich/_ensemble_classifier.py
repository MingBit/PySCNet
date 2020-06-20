#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 18:37:20 2019

@author: angelawu
"""

# from sklearn.model_selection import cross_val_score
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split
# from sklearn import model_selection
from functools import reduce
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier
import autosklearn.classification


# def _generate_x_y(links_dict, threshold):
#
#        for key in links_dict.keys():
#                links_dict[key] = links_dict[key].fillna(0)
#
#        dfs = list(links_dict.values())
#        df_final = reduce(lambda left,right: pd.merge(left,right,on=['source', 'target']), dfs)
#        rep = (df_final.shape[1] - 2)/len(links_dict)
#        df_final.columns = list(df_final.columns[:2]) + list(np.repeat(list(links_dict.keys()), rep))
#        avg = df_final.mean(axis = 1)
#        df_final['Y'] = [(lambda x: 1 if x > threshold else 0)(x) for x in avg]
#
##        X = df_final.iloc[:,:-1].iloc[:,2:]
##        Y = df_final.Y
##
##        node_pairs = X['source'] + '_' + X['target']
#        return (df_final)


def ensemble_classifier(X, Y, threshold=0.5, test_size=0.4, seed=3, model='RF', max_features=5, num_trees=100, **kwarg):
    #    XY = _generate_x_y(links_dict, threshold)
    #
    #    X = XY.iloc[:,:-1]
    #    Y = XY.Y
    x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=test_size)

    #    kfold = model_selection.KFold(n_splits=10, random_state=seed)
    if X.shape[1] - 2 < max_features: max_features = X.shape[1] - 2

    if model == 'RF':

        model = RandomForestClassifier(n_estimators=num_trees, max_features=max_features, **kwarg)


    elif model == 'BDT':
        cart = DecisionTreeClassifier()
        model = BaggingClassifier(base_estimator=cart, n_estimators=num_trees, random_state=seed, **kwarg)


    elif model == 'ET':
        model = ExtraTreesClassifier(n_estimators=num_trees, max_features=max_features, **kwarg)


    elif model == 'AdaB':
        model = AdaBoostClassifier(n_estimators=num_trees, random_state=seed, **kwarg)


    elif model == 'SGB':
        model = GradientBoostingClassifier(n_estimators=num_trees, random_state=seed, **kwarg)

    elif model == 'autoSK':
        model = autosklearn.classification.AutoSklearnClassifier(**kwarg)

    else:
        raise Exception('Valid model: RF, BDT, ET, AdaB, SGB, autoSK')

    model.fit(x_train.iloc[:, 2:], y_train)
    #    pred_train = model.predict(x_train.iloc[:,2:])
    pred_train_prob = model.predict_proba(x_train.iloc[:, 2:])[:, 1]
    #    pred_test = model.predict(x_test.iloc[:,2:])
    pred_test_prob = model.predict_proba(x_test.iloc[:, 2:])[:, 1]
    res_df = pd.DataFrame(
        {'source': x_train.source.append(x_test.source), 'target': x_train.target.append(x_test.target),
         'weight': list(pred_train_prob) + list(pred_test_prob)})
    return (res_df)
