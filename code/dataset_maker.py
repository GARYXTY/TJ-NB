# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 23:23:50 2020

@author: XTY
"""

import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow import feature_column
from tensorflow.keras import layers
import matplotlib.pyplot as plt
from tensorflow import keras
import dataloader




def df_to_dataset(dataframe, shuffle=True, batch_size=32):
    dataframe = dataframe.copy()
    labels = dataframe.pop('status')
    ds = tf.data.Dataset.from_tensor_slices((dict(dataframe), labels))
    if shuffle:
        ds = ds.shuffle(buffer_size=len(dataframe))
    ds = ds.batch(batch_size)
    return ds



def get_datasets(dataset1_train, dataset1_test, dataset2, train_index, val_index):
    train = dataset1_train.iloc[train_index]
    val = dataset1_train.iloc[val_index]
    train_ds = df_to_dataset(train, shuffle=False, batch_size=64)
    val_ds = df_to_dataset(val, shuffle=False, batch_size=64)
    test_ds = df_to_dataset(dataset1_test, shuffle=False, batch_size=64)
    dataset2_ds = df_to_dataset(dataset2, shuffle=False, batch_size=64)
    return train_ds, val_ds, test_ds, dataset2_ds, train, val


