#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 18:37:20 2019

@author: angelawu
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split
from functools import reduce
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn import metrics


def __order_source_target(df):
    
    df_new = pd.DataFrame(tuple(sorted(a)) for a in df[['source', 'target']].values.tolist()). \
            rename(columns={0: 'source',
                            1: 'target'})
    
    df_new['weight'] = df.weight
    
    return df_new


def _generate_x_y(links_dict, top_rank, set_train):
    for key in links_dict.keys():
        links_dict[key] = links_dict[key].fillna(0)
        links_dict[key].weight = links_dict[key].weight.abs()
        links_dict[key] = __order_source_target(links_dict[key])

    dfs = list(links_dict.values())
    df_final = reduce(lambda left, right: pd.merge(left, right, on=['source', 'target']), dfs)
    rep = (df_final.shape[1] - 2) / len(links_dict)
    df_final.columns = list(df_final.columns[:2]) + list(np.repeat(list(links_dict.keys()), rep))
    
    for col in df_final.columns[2:]:
        df_final[col + '_rank'] = list(df_final[col].rank())
        df_final = df_final.drop(col, axis=1)

    df_final['avg'] = df_final.mean(axis=1)
    df_final_sorted = df_final.sort_values('avg', ascending=False, ignore_index=True)

    train_head = df_final_sorted.head(int(df_final.shape[0] * set_train[0]))   
    train_bottom = df_final_sorted.tail(int(df_final.shape[0] * set_train[1]))

    train_df = train_head.append(train_bottom)
    X = train_df.iloc[:, :-1]
    Y = np.concatenate((np.repeat(1, train_head.shape[0]), np.repeat(0, train_bottom.shape[0])))

    df_final_filter = df_final.sort_values('avg', ascending=False, ignore_index=True).head(top_rank)

    return X, Y, df_final_filter.iloc[:, :-1]


def ensemble_classifier(X, Y, df_final, test_size=0.4, seed=3,
                        model='RF', max_features=5, num_trees=100, **kwarg):

    x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=test_size)

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

    else:
        raise Exception('Valid model: RF, BDT, ET, AdaB, SGB')

    model.fit(x_train.iloc[:, 2:].to_numpy(), y_train)
    y_pred = model.predict(x_test.iloc[:, 2:].to_numpy())
    print("Test Accuracy:", metrics.accuracy_score(y_test, y_pred))

    df_final_pred = model.predict(df_final.iloc[:, 2:].to_numpy())

    res_df = pd.DataFrame(
        {'source': df_final.source,
         'target': df_final.target,
         'weight': df_final_pred})

    return res_df[res_df.weight == 1]