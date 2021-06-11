# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 14:03:52 2020

@author: XTY
"""

import dataloader
import dataset_maker
import network
import report_generator
from sklearn.model_selection import train_test_split
import json
import time
import os

start = time.time()
GSE49710_new_T, E_MATAB_new_T, RNA_group0, RNA_group1 = dataloader.dataloader()
print('Data is loaded')
print('=========================================================================')
end = time.time() - start
print(end)
p_test = []
p_train = []
p_data2 = []
p_val = []
auc_test = []
auc_train = []
auc_data2 = []
auc_val = []


for i in range(20):

    train, test = train_test_split(GSE49710_new_T, test_size = 0.3,random_state=i)
    kfold_index = dataloader.Kfold_processor(train, i)
    j = 0
    for train_index, val_index in kfold_index:
        print('the ' + str(i) +'iter')
        print('the '+ str(j) + 'folder')
        start_1 = time.time()
        train_ds, val_ds, test_ds, dataset2_ds, train_1, val = dataset_maker.get_datasets(train, test, E_MATAB_new_T, train_index, val_index)
        model = network.bulid_network(RNA_group0, RNA_group1)
        model = network.train_network(model, train_ds, val_ds)
        report = report_generator.generate_result(train_1, train_ds, val, val_ds, test, test_ds, E_MATAB_new_T, dataset2_ds, model)
        p_train.append(report['train_p_value'])
        p_test.append(report['test_p_value'])
        p_val.append(report['val_p_value'])
        p_data2.append(report['dataset2_p_value'])
        auc_train.append(report['train_auc'])
        auc_val.append(report['val_auc'])
        auc_test.append(report['test_auc'])
        auc_data2.append(report['dataset2_auc'])
        end_1 = time.time() - start_1
        print(end_1)
        
        filefolder = 'C:/Users/XTY/Desktop/FCZ/result2/'
        filename = str(i) + 'iter' + str(j) + 'fold.json'
        report = json.dumps(report)
        path = filefolder + filename
        try:
            if not os.path.exists(path):
                with open(path, "w", encoding='utf-8') as f:
                    f.write(report + "\n")
                    print("^_^ write success")
            else:
                with open(path, "a", encoding='utf-8') as f:
                    f.write(report + "\n")
                    print("^_^ write success")
        except Exception as e:
            print("write error==>", e)
        j = j+1

final_report = dict()
final_report['p_test'] = p_test
final_report['p_train'] = p_train
final_report['p_data2'] = p_data2
final_report['p_val'] = p_val
final_report['auc_test'] = auc_test
final_report['auc_train'] = auc_train
final_report['auc_data2'] = auc_data2
final_report['auc_val'] = auc_val
path = 'C:/Users/XTY/Desktop/FCZ/result2/final_report.json'
final_report = json.dumps(final_report)
try:
    if not os.path.exists(path):
        with open(path, "w", encoding='utf-8') as f:
            f.write(final_report + "\n")
            print("^_^ write success")
    else:
        with open(path, "a", encoding='utf-8') as f:
            f.write(final_report + "\n")
            print("^_^ write success")
except Exception as e:
    print("write error==>", e)     
        
            
            
            
            
            