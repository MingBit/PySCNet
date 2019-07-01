#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 18:37:20 2019

@author: angelawu
"""

#from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
#from sklearn.model_selection import train_test_split
from sklearn import model_selection
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier

def ensemble_classifier(X, Y, seed = 3, model = 'RF', max_features = 5, num_trees = 100):

#    x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size = test_size)    
    kfold = model_selection.KFold(n_splits=10, random_state=seed)    
    
    if model == 'RF':
        
        model = RandomForestClassifier(n_estimators=num_trees, max_features=max_features)
       
        
    elif model == 'BDT':
        cart = DecisionTreeClassifier()
        model = BaggingClassifier(base_estimator=cart, n_estimators=num_trees, random_state=seed)
       

    elif model == 'ET':
        model = ExtraTreesClassifier(n_estimators=num_trees, max_features=max_features)
       
        
    elif model == 'AdaB':
        model = AdaBoostClassifier(n_estimators=num_trees, random_state=seed)
       
        
    elif model == 'SGB':
        model = GradientBoostingClassifier(n_estimators=num_trees, random_state=seed)

    results = model_selection.cross_val_score(model, X, Y, cv=kfold)
        
    return(results)


