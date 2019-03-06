# Tree topology inference from multiple sequence alignments using Deep Learning

This repository contains R and Phython scripts that were used in the project ["Accurate inference of tree topologies from multiple sequence alignments using deep learning"](https://www.biorxiv.org/content/10.1101/559054v1) 

R scripts can be found in INDELible directory. They will generate various control files for [INDElible](http://abacus.gene.ucl.ac.uk/software/indelible/) program that simulates MSAs under given tree topology, branch lengths and different substitution as well as indel model parameters.  

Required R packages:  
[phangorn](https://cran.r-project.org/web/packages/phangorn/index.html)   
[MCMCpack](https://cran.r-project.org/web/packages/MCMCpack/index.html)  
[dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)  
[scales](https://cran.r-project.org/web/packages/scales/index.html)  

1) [indelible_controlgen_INDEL001.R](https://github.com/SchriderLab/Tree_learning/blob/master/INDELible/indelible_controlgen_INDEL001.R) and [indelible_controlgen_NOINDEL.R](https://github.com/SchriderLab/Tree_learning/blob/master/INDELible/indelible_controlgen_NOINDEL.R)    
These scripts generate control files for MSA simulation with (INDEL001) and without (NOINDEL) indels/gaps. The control files will be stored in three directories (topo1, topo2 and topo3) that correspond to three topologies. These scripts are used to generate MSAs for generating TRAINING, VVALIDATION and test data sets. 

2) 
