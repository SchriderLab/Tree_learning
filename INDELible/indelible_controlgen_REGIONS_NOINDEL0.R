#!/usr/bin/env Rscript

### ARGUMENTS: N taxa N_sims Aln_length Parameter of BETA distr

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
#TREE generating function 
tree_gen=function(tr,n_sim,parameter)
{
  nbranch=length(tr$edge[,1])
  tr$edge.lengths=rep(0,nbranch)
  tr_newick=unlist(strsplit(write.tree(tr),""))
  boot=matrix(rbeta(n_sim,parameter,parameter),ncol=nbranch)
  pos=which(tr_newick==0)
  trees=data.frame(t(tr_newick))
  trees=trees[rep(1,n_sim),]
  trees[,pos]=boot[,1:5]
  tree_list=as.vector(apply(trees,1,paste,collapse=""))
  return(tree_list)
}  


tree_genFA=function(tr,n_sim)
{
  nbranch=length(tr$edge[,1])
  tr$edge.lengths=rep(0,nbranch)
  tr_newick=unlist(strsplit(write.tree(tr),""))
  boot=matrix(runif(nbranch*n_sim,c(0.5,0.5,0,0,0,0,0,0,0.5,0.5),c(1,1,1,0.05,0.05,0.05,0.05,1,1,1)),ncol=5,byrow=T)
  pos=which(tr_newick==0)
  trees=data.frame(t(tr_newick))
  trees=trees[rep(1,n_sim),]
  trees[,pos]=boot[,1:5]
  tree_list=as.vector(apply(trees,1,paste,collapse=""))
  return(tree_list)
}
tree_genFE=function(tr,n_sim)
{
  nbranch=length(tr$edge[,1])
  tr$edge.lengths=rep(0,nbranch)
  tr_newick=unlist(strsplit(write.tree(tr),""))
  boot=matrix(runif(nbranch*n_sim,c(0.5,0,0,0.5,0,0.5,0,0,0,0.5,0,0.5,0,0.5,0,0,0.5,0,0,0.5),c(1,0.05,1,1,0.05,1,0.05,1,0.05,1,0.05,1,1,1,0.05,0.05,1,1,0.05,1)),ncol=5,byrow=T)
  pos=which(tr_newick==0)
  trees=data.frame(t(tr_newick))
  trees=trees[rep(1,n_sim),]
  trees[,pos]=boot[,1:5]
  tree_list=as.vector(apply(trees,1,paste,collapse=""))
  return(tree_list)
}

tree_genSHORT=function(tr,n_sim)
{
  nbranch=length(tr$edge[,1])
  tr$edge.lengths=rep(0,nbranch)
  tr_newick=unlist(strsplit(write.tree(tr),""))
  boot=matrix(runif(nbranch*n_sim,c(0,0,0,0,0),c(0.05,0.05,1,0.05,0.05)),ncol=5,byrow=T)
  pos=which(tr_newick==0)
  trees=data.frame(t(tr_newick))
  trees=trees[rep(1,n_sim),]
  trees[,pos]=boot[,1:5]
  tree_list=as.vector(apply(trees,1,paste,collapse=""))
  return(tree_list)
}

tree_genLONG=function(tr,n_sim)
{
  nbranch=length(tr$edge[,1])
  tr$edge.lengths=rep(0,nbranch)
  tr_newick=unlist(strsplit(write.tree(tr),""))
  boot=matrix(runif(nbranch*n_sim,c(0.5,0.5,0,0.5,0.5),c(1,1,1,1,1)),ncol=5,byrow=T)
  pos=which(tr_newick==0)
  trees=data.frame(t(tr_newick))
  trees=trees[rep(1,n_sim),]
  trees[,pos]=boot[,1:5]
  tree_list=as.vector(apply(trees,1,paste,collapse=""))
  return(tree_list)
}

tree_genLONGOUT=function(tr,n_sim)
{
  nbranch=length(tr$edge[,1])
  tr$edge.lengths=rep(0,nbranch)
  tr_newick=unlist(strsplit(write.tree(tr),""))
  boot=matrix(runif(nbranch*n_sim,c(0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0.5,0,0,0,0,0,0.5),c(1,0.05,1,0.05,0.05,0.05,1,1,0.05,0.05,0.05,0.05,1,1,0.05,0.05,0.05,1,0.05,1)),ncol=5,byrow=T)
  pos=which(tr_newick==0)
  trees=data.frame(t(tr_newick))
  trees=trees[rep(1,n_sim),]
  trees[,pos]=boot[,1:5]
  tree_list=as.vector(apply(trees,1,paste,collapse=""))
  return(tree_list)
}

tree_genSHORTOUT=function(tr,n_sim)
{
  nbranch=length(tr$edge[,1])
  tr$edge.lengths=rep(0,nbranch)
  tr_newick=unlist(strsplit(write.tree(tr),""))
  boot=matrix(runif(nbranch*n_sim,c(0,0.5,0,0.5,0.5,0.5,0,0,0.5,0.5,0.5,0.5,0,0,0.5,0.5,0.5,0,0.5,0),c(0.05,1,1,1,1,1,0.05,1,1,1,1,1,1,0.05,1,1,1,1,1,0.05)),ncol=5,byrow=T)
  pos=which(tr_newick==0)
  trees=data.frame(t(tr_newick))
  trees=trees[rep(1,n_sim),]
  trees[,pos]=boot[,1:5]
  tree_list=as.vector(apply(trees,1,paste,collapse=""))
  return(tree_list)
}

tree_genSHORTINT=function(tr,n_sim)
{
  nbranch=length(tr$edge[,1])
  tr$edge.lengths=rep(0,nbranch)
  tr_newick=unlist(strsplit(write.tree(tr),""))
  boot=matrix(runif(nbranch*n_sim,c(0,0,0,0,0),c(1,1,0.05,1,1)),ncol=5,byrow=T)
  pos=which(tr_newick==0)
  trees=data.frame(t(tr_newick))
  trees=trees[rep(1,n_sim),]
  trees[,pos]=boot[,1:5]
  tree_list=as.vector(apply(trees,1,paste,collapse=""))
  return(tree_list)
}


indelib_gen=function(n_taxa,n_sim,aln_length,parameter,region) # n_sim = number of simulations per topology
{
  dir.create(region, showWarnings = FALSE)
  print(paste("I am simulating",n_sim,"alignements per topology of Length =",aln_length,"N taxa =",n_taxa))
  taxa=c('A','B','C','D','E','F','G','H','I','J','K')
  all_topo=allTrees(n_taxa, rooted = FALSE, tip.label = taxa[1:n_taxa])
  iter=0
  for (tr in all_topo)
  {
    iter=iter+1
    dir.create(paste(region,"/topo",iter,sep=""), showWarnings = FALSE)
    write(paste('[TYPE] NUCLEOTIDE 2\n[SETTINGS]\n [output] FASTA\n [randomseed] ',round(runif(1,1,100000))),paste(region,"/topo",iter,'/control.txt',sep=""))
    n_datasets=n_sim
    #Set MODEL block
    modelset=sample(c('JC','TIM','TIMef','GTR','UNREST'),n_datasets,replace=T)
    MODEL=model_gen(modelset,paste(region,"/topo",iter,'/control.txt',sep=""))
    #Set TREE block
    ID_TREE=paste("t",rep(iter,n_sim),rep("_sim",times=n_datasets),1:n_datasets,sep="")
    print(iter)
    print("Newick")
    if (region=="BETA")
    {
      NEWICK=tree_gen(all_topo[[iter]],n_sim,parameter)
    } else if (region == "FA"){
      NEWICK=tree_genFA(all_topo[[iter]],n_sim)
    } else if (region == "FE"){
      NEWICK=tree_genFE(all_topo[[iter]],n_sim)
    } else if (region == "SHORT"){
      NEWICK=tree_genSHORT(all_topo[[iter]],n_sim)
    } else if (region == "LONG"){
      NEWICK=tree_genLONG(all_topo[[iter]],n_sim)
    } else if (region == "LONGOUT"){
      NEWICK=tree_genLONGOUT(all_topo[[iter]],n_sim)
    } else if (region == "SHORTINT"){
      NEWICK=tree_genSHORTINT(all_topo[[iter]],n_sim)
    } else {
      NEWICK=tree_genSHORTOUT(all_topo[[iter]],n_sim)
    }
    print("Done newick")
    write.table(data.frame('[TREE]',ID_TREE,NEWICK),paste(region,"/topo",iter,'/control.txt',sep=""),append=T,quote=F,row.names=F,col.names =F)
    #Set PARTITIONS block
    PNAME=paste("p",1:n_datasets,sep="")
    write.table(data.frame('[PARTITIONS]',PNAME,"[",ID_TREE,MODEL,aln_length,"]"),paste(region,"/topo",iter,'/control.txt',sep=""),append=T,quote=F,row.names=F,col.names =F)
    #Set EVOLVE block
    write('[EVOLVE]',paste(region,"/topo",iter,'/control.txt',sep=""),append=T)
    write.table(data.frame(PNAME,1,apply(data.frame(ID_TREE,"_",MODEL),1,paste,collapse="")),paste(region,"/topo",iter,'/control.txt',sep=""),append=T,quote=F,row.names=F,col.names =F)
  }
}
for (r in c("FA","FE","SHORT","LONG","LONGOUT","SHORTINT","SHORTOUT"))
{  
  indelib_gen(as.numeric(args[1]),as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]),r)
}
