#!/usr/bin/env python3
import tensorflow as tf
sess = tf.Session(config=tf.ConfigProto(log_device_placement=True))
from itertools import product
import sys, argparse, os
import numpy as np
from math import log
from scipy.stats import multinomial, chi2
from math import factorial
from keras.utils import plot_model, to_categorical
from keras.models import Model, Sequential
from keras.layers import Input, Dense, Flatten, Dropout, BatchNormalization
from keras.layers.convolutional import Conv2D
from keras.layers.pooling import MaxPooling2D
from keras.layers.merge import concatenate
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.models import model_from_json
from keras.models import load_model
#N unrooted trees given N taxa
def n_unroot(Ntaxa):
    N=factorial(2*Ntaxa-5)/(factorial(Ntaxa-3)*2**(Ntaxa-3))
    return(int(N))

#Generate Numpy array of labels for Keras
def label_gen(classes,n_labels):
    vec=[[0]*classes]*int(n_labels)
    vec=np.array(vec)
    section=int(int(n_labels)/int(classes))
    for cl in range(0,classes):
        vec[int(cl*section):int(section+section*cl),cl]=1
    return(vec)
    
#Read FASTA convert to numeric
def fasta_pars(aln_file,seq_number,Lmax):
    aln=open(aln_file)
    #Lmax=max([len(r) for r in aln])
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
def tv_parse(train,valid,test,seq_number,Lmax):
    tr=fasta_pars(train,seq_number,Lmax)
    va=fasta_pars(valid,seq_number,Lmax)
    te=fasta_pars(test,seq_number,Lmax)
    print("Training N: "+str(len(tr)))
    print("Validation N: "+str(len(va)))
    print("Testing N: "+str(len(te)))
    return tr, va, te   


#Classification: tree topology
def build_standartCNN(X_train, Y_train, X_valid, Y_valid, Ntaxa,batch_sizes,model,w):
    
    model_cnn=load_model(model)
    model_cnn.load_weights(w)
    #Model stopping criteria
    callback1=EarlyStopping(monitor='val_loss', min_delta=0.001, patience=10, verbose=1, mode='auto')
    callback2=ModelCheckpoint('best_weights_clas', monitor='val_loss', verbose=0, save_best_only=True, save_weights_only=False, mode='auto', period=1)
    #Run

    model_cnn.fit(x=X_train,y=Y_train,validation_data=(X_valid,Y_valid),batch_size=batch_sizes,callbacks=[callback1,callback2],epochs=200,verbose=1,shuffle=True)
    return(model_cnn)



def main():
    parser = argparse.ArgumentParser(description='Keras run')
    parser.add_argument( '-w', help = "Weights file",dest='WEIGHTS')
    parser.add_argument( '-k', help = "Keras model",dest='JASON')
    parser.add_argument( '-t', help = "Training dataset in .fasta",dest='TRAIN')
    parser.add_argument( '-v', help = "Validation dataset in .fasta",dest='VALID')
    parser.add_argument( '--test', help = "Test dataset in .fasta",dest='TEST')
    parser.add_argument( '-N', help = "N taxa", type=int, dest='Ntaxa')
    args = parser.parse_args()
    
    #Read inputs

    loaded_model = load_model(args.JASON)
    loaded_model.summary()
    Lmax=loaded_model.layers[0].get_output_at(0).get_shape().as_list()[2]
    
    print("Reading input")
    train_data1, valid_data1, test_data1 = tv_parse(args.TRAIN,args.VALID,args.TEST,args.Ntaxa,Lmax)    
    #Reshape for Keras
    train_data1=train_data1.reshape(train_data1.shape[0],train_data1.shape[1],train_data1.shape[2],1)
    valid_data1=valid_data1.reshape(valid_data1.shape[0],valid_data1.shape[1],valid_data1.shape[2],1)
    test_data1=test_data1.reshape(test_data1.shape[0],test_data1.shape[1],test_data1.shape[2],1)
    print("Done")
    
    #Generate labels
    Nlabels=n_unroot(args.Ntaxa)
    train_label=to_categorical(np.repeat(range(0,Nlabels),len(train_data1)/Nlabels), num_classes=None)
    valid_label=to_categorical(np.repeat(range(0,Nlabels),len(valid_data1)/Nlabels), num_classes=None)
    test_label=to_categorical(np.repeat(range(0,Nlabels),len(test_data1)/Nlabels), num_classes=None)
    
    #load weights into new model
    model_cnn=build_standartCNN(X_train=train_data1, Y_train=train_label, X_valid=valid_data1, Y_valid=valid_label,Ntaxa=args.Ntaxa,batch_sizes=100,model=args.JASON,w=args.WEIGHTS)
    model_cnn.load_weights('best_weights_clas')
    model_cnn.compile(loss='categorical_crossentropy',optimizer='adam',metrics=['accuracy'])
    
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