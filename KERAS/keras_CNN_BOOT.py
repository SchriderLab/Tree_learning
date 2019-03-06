#!/usr/bin/env python
import tensorflow as tf
sess = tf.Session(config=tf.ConfigProto(log_device_placement=True))
import pandas
from itertools import product
import sys, argparse, os
import numpy as np
from math import log, ceil
from scipy.stats import multinomial, chi2
from math import factorial
from keras.utils import plot_model, to_categorical
from keras.models import Model, Sequential, load_model
from keras.layers import Input, Dense, Flatten, Dropout, BatchNormalization, ZeroPadding2D
from keras.layers.convolutional import Conv2D
from keras.layers.pooling import MaxPooling2D, AveragePooling2D
from keras.layers.merge import concatenate
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.optimizers import Adam

#N unrooted trees given N taxa
def n_unroot(Ntaxa):
    N=factorial(2*Ntaxa-5)/(factorial(Ntaxa-3)*2**(Ntaxa-3))
    return(int(N))


def main():
    parser = argparse.ArgumentParser(description='Keras BOOT run')
    parser.add_argument( '--test', help = "Test dataset in FASTA",dest='TEST')
    parser.add_argument( '--lab', help = "Labels of TEST dataset",dest='LAB')
    parser.add_argument( '-w', help = "Weights file",dest='WEIGHTS')
    parser.add_argument( '-k', help = "Keras model",dest='JASON')
    parser.add_argument( '-b', help = "N bootstrap replicates", type=int, dest='NBoot')
    parser.add_argument( '-N', help = "N taxa", type=int, dest='Ntaxa')    
    args = parser.parse_args()
    
    #Read inputs
    print("Reading input")
    test_data1=np.load(args.TEST)
    print("Done")
    
    #Model and Weight load
    loaded_model = load_model(args.JASON)
    loaded_model.summary()
    loaded_model.load_weights(args.WEIGHTS)
    loaded_model.compile(loss='categorical_crossentropy',optimizer='adam',metrics=['accuracy'])
    
    #Generate labels
    class_lab=np.loadtxt(args.LAB)
    Nlabels=n_unroot(args.Ntaxa)
    #Bootstrapping 
    cnn_boot=[]
    for i in range(0,len(test_data1)):
        aln=test_data1[i]
        booty=[]
        for b in range(0,args.NBoot):
            site=np.argwhere(aln[0,]>=0)[:,0]
            pad=np.argwhere(aln[0,]<0)[:,0]
            boot_site=np.random.choice(site,len(site),replace=True)
            boot_site=np.concatenate((boot_site,pad),axis=None)
            boot_aln=aln[:,boot_site]
            booty.append(boot_aln)
        booty=np.array(booty)
        boot_lab=to_categorical(np.repeat(class_lab[i], args.NBoot),num_classes=Nlabels)
        boot_eval = loaded_model.evaluate(booty,boot_lab,batch_size=300, verbose=1, steps=None)
        cnn_boot.append(boot_eval)
    np.savetxt("test.boots.txt",np.array(cnn_boot)[:,1],fmt='%f')    
    
    
if __name__ == "__main__":
    main()     
    
    
    
    
    
    
    
    
    
    
    
    
    
