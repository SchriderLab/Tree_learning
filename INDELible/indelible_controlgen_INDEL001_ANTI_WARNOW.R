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
    models_selected=c(models_selected,paste(model,'Model',model_id,sep = ''))
    write(paste('\n[MODEL] ',model,'Model',model_id,sep = ''),file,append=T)
    model=model
    model_id=model_id+1
    write(paste(' [submodel] ',paste(model,collapse=' '),'\n [rates] ',0,' ',0,' 0','\n [indelmodel] POW 1.5 50\n [indelrate] 0.01'),file,append=T)
    
  }
  return(models_selected)
}
#TREE generating function 
tree_gen=function(tr,n_sim)
{
  nbranch=length(tr$edge[,1])
  tr$edge.lengths=rep(0,nbranch)
  tr_newick=unlist(strsplit(write.tree(tr),""))
  boot=matrix(runif(nbranch*n_sim,0,0.5),ncol=nbranch)
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
    modelset=rep('JC',n_datasets)
    MODEL=model_gen(modelset,paste(region,"/topo",iter,'/control.txt',sep=""))
    #Set TREE block
    ID_TREE=paste("t",rep(iter,n_sim),rep("_sim",times=n_datasets),1:n_datasets,sep="")
    print(iter)
    print("Newick")
    NEWICK=tree_gen(all_topo[[iter]],n_sim)
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
for (r in c("WARNOW"))
{  
  indelib_gen(as.numeric(args[1]),as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]),r)
}
