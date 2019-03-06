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


def main():
    parser = argparse.ArgumentParser(description='Keras run')
    parser.add_argument( '-t', help = "Evaluation dataset in FASTA",dest='EVAL')
    parser.add_argument( '-w', help = "Weights file",dest='WEIGHTS')
    parser.add_argument( '-k', help = "Keras model",dest='JASON')
    parser.add_argument( '-N', help = "N taxa", type=int, dest='Ntaxa')
    args = parser.parse_args()
    
    loaded_model = load_model(args.JASON)
    loaded_model.summary()
    Lmax=loaded_model.layers[0].get_output_at(0).get_shape().as_list()[2]
    #Read eval
    print("Reading input")
    eval_data=fasta_pars(args.EVAL,args.Ntaxa,Lmax)
    eval_data=eval_data.reshape(eval_data.shape[0],eval_data.shape[1],eval_data.shape[2],1)
    print("Done")
    
    #Generate labels
    Nlabels=n_unroot(args.Ntaxa)
    eval_label=to_categorical(np.repeat(range(0,Nlabels),len(eval_data)/Nlabels), num_classes=None)
    
    # load weights into new model
    print("Loading weights from file")
    loaded_model.load_weights(args.WEIGHTS)
    loaded_model.compile(loss='categorical_crossentropy',optimizer='adam',metrics=['accuracy'])
    
    evals = loaded_model.evaluate(eval_data,eval_label,verbose=1, steps=None)
    classes = loaded_model.predict(eval_data, verbose=1, steps=None)
    np.savetxt("test.evals.txt",evals,fmt='%f')
    np.savetxt("test.classes.txt",classes,fmt='%f')
    class_lab = classes.argmax(axis=-1)
    np.savetxt("test.classeslab.txt",class_lab,fmt='%f')
    
if __name__ == "__main__":
    main() 
