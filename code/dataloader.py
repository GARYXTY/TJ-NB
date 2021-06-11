# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 22:28:47 2020

@author: XTY
"""

import pandas as pd
import numpy as np
from sklearn.feature_selection import chi2
from sklearn.cluster import KMeans
from sklearn.model_selection import KFold

#### Download and Pre-process Data
def dataloader():
    E_MATAB = pd.read_csv('../data/E-MTAB/dat_EMTAB.csv')
    E_MATAB_ph = pd.read_csv('../data/E-MTAB/phenotype_EMTAB.csv')
    GSE49710 = pd.read_csv('../data/GSE49710/GSE49710_matrix.csv')   
    GSE49710_ph = pd.read_csv('../data/GSE49710/GSE49710_phenotype.csv')
    GSE85047 = pd.read_csv('../data/GSE85047/GSE85047_matrix.csv')
    
    
    E_MATAB_name = list(E_MATAB['Unnamed: 0'])
    GSE49710_name = list(GSE49710['Unnamed: 0'])
    GSE85047_name = list(GSE85047['Unnamed: 0'])
    cross = set(E_MATAB_name)&set(GSE49710_name)&set(GSE85047_name)
    
    droplist = []
    cross = list(cross)
    for i,j in enumerate(GSE49710_name):
        if j not in cross:
            droplist.append(i)
    GSE49710_new = GSE49710.drop(labels=droplist)
    droplist = []
    cross = list(cross)
    for i,j in enumerate(E_MATAB_name):
        if j not in cross:
            droplist.append(i)
    E_MATAB_new = E_MATAB.drop(labels=droplist)
    
    GSE49710_new_T = pd.DataFrame(GSE49710_new.T, index=GSE49710_new.columns, columns=GSE49710_new.index) 
    E_MATAB_new_T = pd.DataFrame(E_MATAB_new.T, index=E_MATAB_new.columns, columns=E_MATAB_new.index)
    E_MATAB_new_T.columns = E_MATAB_new['Unnamed: 0']
    GSE49710_new_T.columns = GSE49710_new['Unnamed: 0']
    
    GSE49710_status = GSE49710_ph.pop('os_bin')
    GSE49710_new_T = GSE49710_new_T.drop(labels='Unnamed: 0')
    E_MATAB_new_T = E_MATAB_new_T.drop(labels='Unnamed: 0')
    #E_MATAB_new_T = (E_MATAB_new_T-GSE49710_new_T.min())/(E_MATAB_new_T.max()-GSE49710_new_T.min())
    GSE49710_new_T = (GSE49710_new_T-GSE49710_new_T.min())/(GSE49710_new_T.max()-GSE49710_new_T.min())
    E_MATAB_new_T = (E_MATAB_new_T-E_MATAB_new_T.min())/(E_MATAB_new_T.max()-E_MATAB_new_T.min())
    
    for i in GSE49710_new_T.columns:
        temp = GSE49710_new_T[i]
        temp1 = np.array(temp)
        temp1 = np.expand_dims(temp1, axis =1)
        x = chi2(temp1, GSE49710_status)
        if x[1]>0.05:
            GSE49710_new_T.pop(i)
    E_MATAB_status = E_MATAB_ph['Characteristics.overall.survival..1.dead..0.alive.w.o.event..']
    E_MATAB_time = E_MATAB_ph['Characteristics.overall.survival.']

    for i in E_MATAB_new_T.columns:
        if i not in GSE49710_new_T.columns:
            E_MATAB_new_T.pop(i)
    
    
    test = np.asarray(GSE49710_new_T)
    test = np.transpose(test)
    kmeans = KMeans(n_clusters=2, random_state=0).fit(test)
    
    RNA_group0 = []
    RNA_group1 = []
    for i in range(len(list(kmeans.labels_))):
        if list(kmeans.labels_)[i] == 0:
            RNA_group0.append(E_MATAB_new_T.columns[i])
        else:
            RNA_group1.append(E_MATAB_new_T.columns[i])
    
    GSE49710_new_T['time'] = 0.0
    E_MATAB_new_T['time'] = 0.0
    GSE49710_new_T['status'] = 0
    E_MATAB_new_T['status'] = 0.0
    GSE49710_time = GSE49710_ph.pop('os_days')
    GSE49710_time.astype(np.float64)
    E_MATAB_time.astype(np.float64)
    for i in range(498):
        GSE49710_new_T['status'][i] = int(GSE49710_status[i])
        GSE49710_new_T['time'][i] = int(GSE49710_time[i]/365.0)
    for i in range(len(E_MATAB_status)):
        E_MATAB_new_T['time'][i] = E_MATAB_time[i]
        E_MATAB_new_T['status'][i] = E_MATAB_status[i]
    
    GSE49710_new_T = GSE49710_new_T.astype(np.float64)
    GSE49710_new_T['status'] = GSE49710_new_T['status'].astype('int')
    E_MATAB_new_T = E_MATAB_new_T.astype(np.float64)
    E_MATAB_new_T['status'] = E_MATAB_new_T['status'].astype('int')
    
    return GSE49710_new_T, E_MATAB_new_T, RNA_group0, RNA_group1


    
# K-fold process
def Kfold_processor(dataset, seed):
    kf = KFold(n_splits = 10, shuffle= True, random_state = seed)
    kf.get_n_splits(dataset)
    return kf.split(dataset)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    