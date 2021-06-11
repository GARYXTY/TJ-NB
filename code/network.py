# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 23:40:30 2020

@author: XTY
"""

import pandas as pd
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from tensorflow import feature_column
from tensorflow.keras import layers
import matplotlib.pyplot as plt
from tensorflow import keras
from sklearn.feature_selection import f_regression
import tensorflow.keras.backend as K
import dataloader

def bulid_network(RNA_group0, RNA_group1):
    feature_columns = []
    feature_layer_inputs = {}
    temp = RNA_group0 + RNA_group1
    for header in temp:
        feature_columns.append(feature_column.numeric_column(header))
        feature_layer_inputs[header] = tf.keras.Input(shape=(1,), name=header)
    
    feature_layer = tf.keras.layers.DenseFeatures(feature_columns)
    feature_layer_output = feature_layer(feature_layer_inputs)

    input1 = feature_layer_output[:,0:122]
    input2 = feature_layer_output[:,122:]
    
    x1 = layers.Dense(100, activation='relu', name='dense_1_1')(input1)
    x1 = layers.Dropout(0.2)(x1)
    x1= layers.Dense(50, activation='relu', name='dense_1_2')(x1)
    x1 = layers.Dropout(0.2)(x1)
    x1= layers.Dense(10, activation='relu', name='dense_1_3')(x1)
    
    x2 = layers.Dense(100, activation='relu', name='dense_2_1')(input2)
    x2 = layers.Dropout(0.2)(x2)
    x2 = layers.Dense(50, activation='relu', name='dense_2_2')(input2)
    x2 = layers.Dropout(0.2)(x2)
    x2= layers.Dense(10, activation='relu', name='dense_2_3')(x2)
    
    x = keras.layers.concatenate([tf.sigmoid(x1)*x2, tf.sigmoid(x2)*x1],axis=1)
    
    x = layers.Dense(10, activation='relu', name='dense_3_1')(x)
    x = layers.Dense(10, activation='relu', name='dense_3_2')(x)
    
    
    outputs = layers.Dense(1, activation = 'sigmoid')(x1)
    model = keras.Model(inputs=[v for v in feature_layer_inputs.values()], outputs=outputs)
    return model



def focal_loss(gamma=2., alpha=.25):
	def focal_loss_fixed(y_true, y_pred):
		pt_1 = tf.where(tf.equal(y_true, 1), y_pred, tf.ones_like(y_pred))
		pt_0 = tf.where(tf.equal(y_true, 0), y_pred, tf.zeros_like(y_pred))
		return -K.mean(alpha * K.pow(1. - pt_1, gamma) * K.log(pt_1)) - K.mean((1 - alpha) * K.pow(pt_0, gamma) * K.log(1. - pt_0))
	return focal_loss_fixed


def train_network(model, train_ds, val_ds):
    earlystop=tf.keras.callbacks.EarlyStopping(monitor='val_accuracy',min_delta = -0.04,patience=5,verbose=True,restore_best_weights=True)
    callbacks=[earlystop]
    model.compile(optimizer='adam',
              loss = [focal_loss(alpha=.2, gamma=2)],
              metrics=['mse','mae','accuracy'],
              run_eagerly=True)
    model.fit(train_ds,
          validation_data=val_ds,
          epochs=100)
    return model
    


