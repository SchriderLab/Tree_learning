#!/usr/bin/env Rscript

### ARGUMENTS: N taxa N_sims Aln_length

args = commandArgs(trailingOnly=TRUE)
library('phangorn')
library('MCMCpack')
library('dplyr')
library('scales')
options(scipen=999)
#Model block generating function
model_gen=function(modelset,file)
{
  model_id=1
  models_selected=c()
  for (model in modelset)
  {
    model_orig=model
    #Invariant sites Unif
    I=runif(1,0,1)
    A=runif(1,0,5)
    #Nucl proportions DIRICHLET 
    options(digits=5)
    Pi=format(rdirichlet(1, alpha=c(5,5,5,5)))
    models_selected=c(models_selected,paste(model,'Model',model_id,sep = ''))
    write(paste('\n[MODEL] ',model,'Model',model_id,sep = ''),file,append=T)
    options(digits=2)
    if (model %in% c('HKY','K80')){
      model=paste(c(model,' ',format(runif(1,0,3))),sep = '')
    } else if (model == 'TrN'){
      model=paste(c(model,' ',format(runif(2,0,3))),sep = '')
    } else if (model %in% c('TIM' ,'TIMef')){
      model=paste(c(model,' ',format(runif(3,0,3))),sep = '')
    } else if (model == 'TVM'){
      model=paste(c(model,' ',format(runif(4,0,3))),sep = '')
    } else if (model %in% c('SYM','GTR')){
      model=paste(c(model,' ',format(runif(5,0,3))),sep = '')
    } else if (model == 'UNREST'){
      model=paste(c(model,' ',format(runif(11,0,3))),sep = '')
    } else {
      model=model
    }
    model_id=model_id+1
    write(paste(' [submodel] ',paste(model,collapse=' '),'\n [rates] ',I,' ',A,' 0'),file,append=T)
    if (model_orig %in% c('F81','HKY','TrN','TIM','TVM','GTR'))
    {
      write(paste(' [statefreq]',paste(Pi,collapse=' ')),file,append=T)
    }
  }
  return(models_selected)
}


tree_genFE_SHORTULTRA=function(tr,n_sim)
{
  nbranch=length(tr$edge[,1])
  tr$edge.lengths=rep(0,nbranch)
  tr_newick=unlist(strsplit(write.tree(tr),""))
  boot1=matrix(runif(nbranch*n_sim,c(0.1,0,0,0.1,0,0.1,0,0,0,0.1,0,0.1,0,0.1,0,0,0.1,0,0,0.1),c(0.5,0.05,0.05,0.5,0.05,0.5,0.05,0.05,0.05,0.5,0.05,0.5,0.05,0.5,0.05,0.05,0.5,0.05,0.05,0.5)),ncol=5,byrow=T)
  boot2=matrix(runif(nbranch*n_sim,c(0,0,0,0,0),c(0.01,0.01,0.01,0.01,0.01)),ncol=5,byrow=T)
  boot=rbind(boot1,boot2)
  boot=boot[sample(1:nrow(boot),n_sim),]
  pos=which(tr_newick==0)
  trees=data.frame(t(tr_newick))
  trees=trees[rep(1,n_sim),]
  trees[,pos]=boot[,1:5]
  tree_list=as.vector(apply(trees,1,paste,collapse=""))
  return(tree_list)
}




indelib_gen=function(n_taxa,n_sim,aln_length) # n_sim = number of simulations per topology
{
  print(paste("I am simulating",n_sim,"alignements per topology of Length =",aln_length,"N taxa =",n_taxa))
  taxa=c('A','B','C','D','E','F','G','H','I','J','K')
  all_topo=allTrees(n_taxa, rooted = FALSE, tip.label = taxa[1:n_taxa])
  iter=0
  for (tr in all_topo)
  {
    iter=iter+1
    write(paste('[TYPE] NUCLEOTIDE 2\n[SETTINGS]\n [output] FASTA\n [randomseed] ',round(runif(1,1,100000))),paste(iter,'control.txt',sep=""))
    n_datasets=n_sim
    #Set MODEL block
    modelset=sample(c('JC','TIM','TIMef','GTR','UNREST'),n_datasets,replace=T)
    MODEL=model_gen(modelset,paste(iter,'control.txt',sep=""))
    #Set TREE block
    ID_TREE=paste("t",rep(iter,n_sim),rep("_sim",times=n_datasets),1:n_datasets,sep="")
    print(iter)
    print("Newick")
    NEWICK=tree_genFE_SHORTULTRA(all_topo[[iter]],n_sim)
    print("Done newick")
    write.table(data.frame('[TREE]',ID_TREE,NEWICK),paste(iter,'control.txt',sep=""),append=T,quote=F,row.names=F,col.names =F)
    #Set PARTITIONS block
    PNAME=paste("p",1:n_datasets,sep="")
    write.table(data.frame('[PARTITIONS]',PNAME,"[",ID_TREE,MODEL,aln_length,"]"),paste(iter,'control.txt',sep=""),append=T,quote=F,row.names=F,col.names =F)
    #Set EVOLVE block
    write('[EVOLVE]',paste(iter,'control.txt',sep=""),append=T)
    write.table(data.frame(PNAME,1,apply(data.frame(ID_TREE,"_",MODEL),1,paste,collapse="")),paste(iter,'control.txt',sep=""),append=T,quote=F,row.names=F,col.names =F)
  }
}

indelib_gen(as.numeric(args[1]),as.numeric(args[2]),as.numeric(args[3]))