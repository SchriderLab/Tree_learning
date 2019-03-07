# Tree topology inference from multiple sequence alignments using Deep Learning

This repository contains R (>=3.5.0) and Python (>=3.6) scripts that were used in the project ["Accurate inference of tree topologies from multiple sequence alignments using deep learning"](https://www.biorxiv.org/content/10.1101/559054v1) 

**R scripts** can be found in [INDELible directory](https://github.com/SchriderLab/Tree_learning/tree/master/INDELible). They will generate various control files for [INDElible](http://abacus.gene.ucl.ac.uk/software/indelible/) program that simulates MSAs under given tree topology, branch lengths and different substitution as well as indel model parameters.  

Required CRAN R packages:  
[phangorn](https://cran.r-project.org/web/packages/phangorn/index.html)   
[MCMCpack](https://cran.r-project.org/web/packages/MCMCpack/index.html)  
[dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)  
[scales](https://cran.r-project.org/web/packages/scales/index.html)  

1) [indelible_controlgen_INDEL001.R](https://github.com/SchriderLab/Tree_learning/blob/master/INDELible/indelible_controlgen_INDEL001.R) and [indelible_controlgen_NOINDEL.R](https://github.com/SchriderLab/Tree_learning/blob/master/INDELible/indelible_controlgen_NOINDEL.R)    
These scripts generate control files for MSA simulation with (INDEL001) and without (NOINDEL) indels/gaps. The control files will be stored in three directories (topo1, topo2 and topo3) that correspond to three topologies. These scripts are used to generate MSAs for generating TRAINING, VALIDATION and TEST data sets.  
Example: ```Rscript indelible_controlgen_NOINDEL.R 4 1000 500```(generates 1000 MSAs of length 500 per topology) 

2) [indelible_controlgen_REGIONS_INDEL001.R](https://github.com/SchriderLab/Tree_learning/blob/master/INDELible/indelible_controlgen_REGIONS_INDEL001.R) and [indelible_controlgen_REGIONS_NOINDEL0.R](https://github.com/SchriderLab/Tree_learning/blob/master/INDELible/indelible_controlgen_REGIONS_NOINDEL0.R)  
These scripts generate control files for MSA simulation with (INDEL001) and without (NOINDEL) indels/gaps. The scripts will generate FA, FE, LONG, LONGOUT, SHORT, SHORTINT, SHORTOUT directories each with topo1, topo2 and topo3 subdirectories. These correspond to heterogeneous branch length regions, namely Extended Farris zone (FA), Extended Felsenstein zone (FE), Long branches (LONG), Single long branch (LONGOUT), Short branches (SHORT), Short internal branch (SHORTINT) and Single short branch (SHORTOUT). These MSAs were used to test performance of different tree inference methods.  
Example: ```Rscript indelible_controlgen_REGIONS_INDEL001.R 4 1000 500```(generates 1000 MSAs of length 500 per topology for each region)      

3) [indelible_controlgen_REGIONS_INDEL001_FATFAFE.R](https://github.com/SchriderLab/Tree_learning/blob/master/INDELible/indelible_controlgen_REGIONS_INDEL001_FATFAFE.R) and [indelible_controlgen_REGIONS_NOINDEL_FATFAFE.R](https://github.com/SchriderLab/Tree_learning/blob/master/INDELible/indelible_controlgen_REGIONS_NOINDEL_FATFAFE.R)  
These scripts generate control files for MSA simulation with (INDEL001) and without (NOINDEL) indels/gaps. The scripts will generate FA, FAT, FE, and SHORTINT directories each with topo1, topo2 and topo3 subdirectories. These correspond to heterogeneous branch length regions, namely Farris zone (FA), Felsenstein zone (FE), ["Twisted" Farris zone](https://www.sciencedirect.com/science/article/pii/S1055790315002316?via%3Dihub) (FAT) and Short internal branch (SHORTINT). These MSAs were used to test performance of different tree inference methods.  
Example: ```Rscript indelible_controlgen_REGIONS_INDEL001_FATFAFE.R 4 1000 500```(generates 1000 MSAs of length 500 per topology for each region) 

4) [indelible_controlgen_INDEL001_WARNOW.R](https://github.com/SchriderLab/Tree_learning/blob/master/INDELible/indelible_controlgen_INDEL001_WARNOW.R)   
This script generates control files for MSA simulation with no substitutions, only indels. This is the scenario under which maximum likelihood (ML) tree inference has been shown to be statistically inconsistent ([Warnow, 2012](http://currents.plos.org/treeoflife/index.html%3Fp=1609.html)). These MSAs were used to test performance of different tree inference methods.   
Example: ```Rscript indelible_controlgen_INDEL001_WARNOW.R 4 1000 500``` (generates 1000 MSAs of length 500 per topology)   

**Python scripts** can be found in [KERAS directory](https://github.com/SchriderLab/Tree_learning/tree/master/KERAS). They are used for building, training, validating and testing Convolutional Neuronal Networks (CNNs). These scripts are optimized to run on GPUs.  

Required Python dependencies:  
[Tensorflow](https://www.tensorflow.org/install)  
[Keras API](https://keras.io/)      
[SciPy](https://www.scipy.org/)    
[pandas](https://pandas.pydata.org/) 

1) [keras_CNN_TOPO.py](https://github.com/SchriderLab/Tree_learning/blob/master/KERAS/keras_CNN_TOPO.py)   
This script builds, trains, validates and tests CNN. As an input it takes TRAINING, VALIDATION and TESTING MSAs generated by INDELible and saved in .npy array using [fasta2numeric.py](https://github.com/SchriderLab/Tree_learning/tree/master/Utils) script. As an input this utility script takes TRAINING, VALIDATION and TESTING datasets produced by concatinating MSAs. E.g. ```cat topo1/* topo2/* topo3/* > TRAINING``` The keras model (keras.h5) and optimal CNN weights (best_weights_clas) will be outputted by the script after testing is completed.   
Example: ```keras_CNN_TOPO.py -t TRAIN.npy -v VALID.npy  --test TEST.npy  -N 4``` (tested only on 4-taxon MSA cases i.e. ```-N 4```)  
  
 ```
 Options:
  -h, --help   
  -t Training dataset in .npy
  -v Validation dataset in .npy
  --test Test dataset in .npy
  -N N taxa 
  ```

2) [keras_CNN_apply.py](https://github.com/SchriderLab/Tree_learning/blob/master/KERAS/keras_CNN_apply.py)  
This script infers a tree from an MSA. It requires keras model and weights files produced by [keras_CNN_TOPO.py](https://github.com/SchriderLab/Tree_learning/blob/master/KERAS/keras_CNN_TOPO.py), a data set in FASTA format.  
Example: ```keras_CNN_apply.py -t TEST.fasta -w best_weights_clas -k keras_model.h5 -N 4```
```
Options:
  -h, --help 
  -t Evaluation dataset in FASTA
  -w Weights file
  -k Keras model JASON
  -N N taxa
```
3) [keras_CNN_BOOT.py](https://github.com/SchriderLab/Tree_learning/blob/master/KERAS/keras_CNN_BOOT.py)  
This script performs MSA nonparametric bootstrapping. It requires keras model and weights files produced by [keras_CNN_TOPO.py](https://github.com/SchriderLab/Tree_learning/blob/master/KERAS/keras_CNN_TOPO.py), a data set in FASTA format and labeles for  data set.   
Example: ```keras_CNN_BOOT.py --test TEST.fasta --lab labels.txt -w best_weights_clas -k keras_model.h5 -b 100 -N 4```   
```
Options:
  -h, --help
  --test Test dataset in FASTA
  --lab Labels of TEST dataset
  -w Weights file
  -k Keras model
  -b N bootstrap replicates
  -N N taxa
```
