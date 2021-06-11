# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 23:49:08 2020

@author: XTY
"""


import lifelines
from lifelines.statistics import logrank_test
from sklearn.metrics import roc_auc_score
import math


def compute_p_value(pp, status, survive_time):
    for i in range(len(pp)):
        if pp[i] < 0.5:
            pp[i] = 0
        else:
            pp[i] = 1
    group1_survive = []
    group1_status = []
    group0_survive = []
    group0_status = []
    
    for i in range(len(pp)):
        if pp[i] == 1:
            group1_status.append(status[i])
            group1_survive.append(survive_time[i])
        else:
            group0_status.append(status[i])
            group0_survive.append(survive_time[i])
    result = logrank_test(group1_survive, group0_survive, group1_status, group0_status)
    p_value = result.summary['p'][0]
    return p_value

def predict_result(dataset, pp):
    result = dict()
    for i,j in enumerate(dataset.index):
        result[j] = pp[i].tolist()
    return result
            

def generate_result(train, train_ds, val, val_ds, test, test_ds, dataset2, dataset2_ds, model):
    pp_train = model.predict(train_ds)
    pp_val = model.predict(val_ds)
    pp_test = model.predict(test_ds)
    pp_dataset2 = model.predict(dataset2_ds)
    report = dict()
    
    
    
    train_status = train['status']
    val_status = val['status']
    test_status = test['status']
    dataset2_status = dataset2['status']
    
    train_result = predict_result(train, pp_train)
    test_result = predict_result(test, pp_test)
    val_result = predict_result(val, pp_val)
    dataset2_result = predict_result(dataset2, pp_dataset2)
    
    dataset1_result = dict()
    dataset1_result['train_result'] = train_result
    dataset1_result['test_result'] = test_result    
    dataset1_result['val_result'] = val_result
    
    report['dataset1_result'] = dataset1_result
    report['dataset2_result'] = dataset2_result
    
    train_auc = roc_auc_score(train_status, pp_train)
    val_auc = roc_auc_score(val_status, pp_val)
    test_auc = roc_auc_score(test_status, pp_test)
    dataset2_auc = roc_auc_score(dataset2_status, pp_dataset2)
    
    train_survive_time = train['time']
    test_survive_time = test['time']
    val_survive_time = val['time']
    dataset2_survive_time = dataset2['time']
    
    train_p_value = compute_p_value(pp_train, train_status, train_survive_time)
    test_p_value = compute_p_value(pp_test, test_status, test_survive_time)
    val_p_value = compute_p_value(pp_val, val_status, val_survive_time)
    datatset2_p_value = compute_p_value(pp_dataset2, dataset2_status, dataset2_survive_time)
    
    report['train_p_value'] = train_p_value
    report['test_p_value'] = test_p_value
    report['val_p_value'] = val_p_value
    report['dataset2_p_value'] = datatset2_p_value
          
    

    
    report['train_auc'] = train_auc
    report['val_auc'] = val_auc
    report['test_auc'] = test_auc
    report['dataset2_auc'] = dataset2_auc
    

    
    return report
    
    

         
    
    
    
    
    
    