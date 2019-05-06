#!/usr/bin/env python
import tensorflow as tf
sess = tf.Session(config=tf.ConfigProto(log_device_placement=True))
import pandas
import time
from itertools import product
import sys, argparse, os
import numpy as np
from math import log, ceil
from scipy.stats import multinomial, chi2
from math import factorial
from keras.utils import plot_model, to_categorical
from keras.models import Model, Sequential
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

#Read FASTA convert to numeric
def fasta_pars(aln_file,seq_number,Lmax):
    aln=open(aln_file)
    dic={'A':'0','T':'1','C':'2','G':'3','-':'4'}
    matrix_out=[]
    fasta_dic={}
    for line in aln:
        if line[0]==">":
            header=line[1:].rstrip('\n').strip()
            fasta_dic[header]=[]
        elif line[0].isalpha() or line[0]=='-':
            for base, num in dic.items():
                line=line[:].rstrip('\n').strip().replace(base,num)
            line=list(line)
            line=[int(n) for n in line]
            #Mkae all matrices of equal length for CNN +[-15]*(Lmax-len(line)) 
            fasta_dic[header]+=line+[-15]*(Lmax-len(line)) 
            if len(fasta_dic)==seq_number:
                taxa_block=[]
                for taxa in sorted(list(fasta_dic.keys())):
                    taxa_block.append(fasta_dic[taxa.strip()])
                fasta_dic={}
                matrix_out.append(taxa_block)
    return(np.array(matrix_out))

#Read training, validation and test datasets to equalize sizes 
def tv_parse(train,valid,test,seq_number):
    tr=open(train)
    va=open(valid)
    te=open(test)
    LT=max([len(r.strip()) for r in tr])
    print("Training largest alignment: "+str(LT))
    LV=max([len(r.strip()) for r in va])
    print("Validation largest alignment: "+str(LV))
    LTE=max([len(r.strip()) for r in te])
    print("Testing largest alignment: "+str(LTE))
    Lmax=max([LT]+[LV]+[LTE])
    tr=fasta_pars(train,seq_number,Lmax)
    va=fasta_pars(valid,seq_number,Lmax)
    te=fasta_pars(test,seq_number,Lmax)
    return tr, va, te   


#Classification: tree topology
def build_standartCNN(X_train, Y_train, X_valid, Y_valid, Ntaxa,conv_pool_n,filter_n,droput_rates,batch_sizes):
    Aln_length=X_train.shape[2]
    Nlabels=n_unroot(Ntaxa)
    
    #Hyperparameters
    #Hight (horizontal)
    conv_x=[4,1,1,1,1,1,1,1]
    #Width (vertical)
    conv_y=[1,2,2,2,2,2,2,2]
    pool=[1,4,4,4,2,2,2,1]
    filter_s=[1024,1024,128,128,128,128,128,128]
    
    print(conv_pool_n,filter_n,droput_rates,batch_sizes)
   
    #Arhitecture
    visible = Input(shape=(Ntaxa,Aln_length,1))
    x = visible
    for l in list(range(0,conv_pool_n)):
        x = ZeroPadding2D(padding=((0, 0), (0,conv_y[l]-1)))(x)        
        x = Conv2D(filters=filter_s[l], kernel_size=(conv_x[l], conv_y[l]), strides=1, activation='relu')(x)
        x = BatchNormalization()(x)
        x = Dropout(rate=0.2)(x)
        x = AveragePooling2D(pool_size=(1,pool[l]))(x)
    flat = Flatten()(x)
    hidden1 = Dense(1024,activation='relu')(flat)
    drop1=Dropout(rate=0.2)(hidden1)
    output = Dense(Nlabels, activation='softmax')(drop1)
    model_cnn = Model(inputs=visible, outputs=output)
    model_cnn.compile(loss='categorical_crossentropy',optimizer='adam',metrics=['accuracy'])
    
    #Print model
    print(model_cnn.summary())
    
    #Model stopping criteria
    callback1=EarlyStopping(monitor='val_loss', min_delta=0.001, patience=10, verbose=1, mode='auto')
    callback2=ModelCheckpoint('best_weights_clas', monitor='val_loss', verbose=0, save_best_only=True, save_weights_only=False, mode='auto', period=1)
    #Run

    model_cnn.fit(x=X_train,y=Y_train,validation_data=(X_valid,Y_valid),batch_size=batch_sizes,callbacks=[callback1,callback2],epochs=200,verbose=1,shuffle=True)
    return(model_cnn)
 

def main():
    parser = argparse.ArgumentParser(description='Keras run')
    parser.add_argument( '-t', help = "Training dataset in .npy",dest='TRAIN')
    parser.add_argument( '-v', help = "Validation dataset in .npy",dest='VALID')
    parser.add_argument( '--test', help = "Test dataset in .npy",dest='TEST')
    parser.add_argument( '-N', help = "N taxa", type=int, dest='Ntaxa')
    args = parser.parse_args()
    
    #Read inputs
    print("Reading input")
 
    train_data1=np.load(args.TRAIN)
    valid_data1=np.load(args.VALID)
    test_data1=np.load(args.TEST)
    print("Done")
    
    #Generate labels
    Nlabels=n_unroot(args.Ntaxa)
    train_label=to_categorical(np.repeat(range(0,Nlabels),len(train_data1)/Nlabels), num_classes=None)
    valid_label=to_categorical(np.repeat(range(0,Nlabels),len(valid_data1)/Nlabels), num_classes=None)
    test_label=to_categorical(np.repeat(range(0,Nlabels),len(test_data1)/Nlabels), num_classes=None)
 
    
    #Model Run
    #Classification TOPO
    starttime=time.time()
    model_cnn=build_standartCNN(X_train=train_data1, Y_train=train_label, X_valid=valid_data1, Y_valid=valid_label,Ntaxa=args.Ntaxa,conv_pool_n=8,filter_n=500,droput_rates=0.20,batch_sizes=200)
    #Load best model
    model_cnn.load_weights('best_weights_clas')
    model_cnn.compile(loss='categorical_crossentropy',optimizer='adam',metrics=['accuracy'])
    model_cnn.save("keras_model.h5")
    endtime=time.time()
    print(endtime-starttime,"sec for training")
    print('Evaluate with best class weights')
   

    starttime=time.time()
    evals = model_cnn.evaluate(test_data1,test_label,batch_size=100, verbose=1, steps=None)
    classes = model_cnn.predict(test_data1, batch_size=100, verbose=1, steps=None)
    endtime=time.time()
    print(endtime-starttime,"sec for testing")
    np.savetxt("test.evals_class.txt",evals,fmt='%f')
    np.savetxt("test.classprobs_class.txt",classes,fmt='%f')
    class_lab = classes.argmax(axis=-1)
    np.savetxt("test.classeslab_class.txt",class_lab,fmt='%f')
    
if __name__ == "__main__":
    main()     
    
    
    
    
    
    
    
    
    
    
    
    
    
