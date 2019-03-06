#!/usr/bin/env python3
import pandas
from itertools import product
import sys, argparse, os
import numpy as np
from math import log, ceil
from scipy.stats import multinomial, chi2
from math import factorial

###Convert TRAIN VALID TEST to numeric and save as nupy array

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
    print("Training N: "+str(len(tr)))
    print("Validation N: "+str(len(va)))
    print("Testing N: "+str(len(te)))
    return tr, va, te   
    
def main():
    parser = argparse.ArgumentParser(description='fasta2numeric conversion Ready for Keras')
    parser.add_argument( '-t', help = "Training dataset in FASTA",dest='TRAIN')
    parser.add_argument( '-v', help = "Validation dataset in FASTA",dest='VALID')
    parser.add_argument( '--test', help = "Test dataset in FASTA",dest='TEST')
    parser.add_argument( '-N', help = "N taxa", type=int, dest='Ntaxa')
    args = parser.parse_args()
    
    print("Reading input")
    train_data1, valid_data1, test_data1 = tv_parse(args.TRAIN,args.VALID,args.TEST,args.Ntaxa)    
    #Reshape for Keras
    train_data1=train_data1.reshape(train_data1.shape[0],train_data1.shape[1],train_data1.shape[2],1)
    valid_data1=valid_data1.reshape(valid_data1.shape[0],valid_data1.shape[1],valid_data1.shape[2],1)
    test_data1=test_data1.reshape(test_data1.shape[0],test_data1.shape[1],test_data1.shape[2],1)
    np.save("TRAIN", train_data1)
    np.save("VALID", valid_data1)
    np.save("TEST", test_data1)
     
if __name__ == "__main__":
    main()  
    
    
    
    
    
    
    
    
    
