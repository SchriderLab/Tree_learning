library("MASS")
library("ape")
library("phangorn")
library("ggplot2")
library("gtools")
library('pals')
library('vioplot')
library('dplyr')
library('scales')
library('geometry')
library('gplots')
library("corrplot")
library("ggtern")
library("ggpubr")

plot_topo=function(brl,namezone,xlabz,ylabz)
{
  plot(1,1,xlim=c(0.4,1.6),ylim=c(0.6,1.4),col="white",xaxt='n',yaxt='n',main=namezone,xlab="",ylab="")
  brl=ifelse(brl==0,0.15,0.45)
  ext=brl
  b1=ext[1];b2=ext[2];b=ext[3];b3=ext[4];b4=ext[5];as=1;lns=1
  l=b/2
  ax=as-l-sqrt(b1^2/2)
  bx=as-l
  cx=as+l
  dx=as+l+sqrt(b3^2/2)
  ex=as+l+sqrt(b4^2/2)
  fx=as-l-sqrt(b2^2/2)
  ay=sqrt(b1^2/2)+lns
  by=lns
  cy=lns
  dy=sqrt(b3^2/2)+lns
  ey=lns-sqrt(b4^2/2)
  fy=lns-sqrt(b2^2/2) 
  lines(c(ax,bx,cx,dx,cx,ex,cx,bx,fx),c(ay,by,cy,dy,cy,ey,cy,by,fy),type="l",lwd=1.5)
  lines(ax,ay,pch=21,cex=2,type='p',bg="white")
  lines(fx,fy,pch=21,cex=2,type='p',bg="white")
  lines(dx,dy,pch=21,cex=2,type='p',bg="white")
  lines(ex,ey,pch=21,cex=2,type='p',bg="white")
  lines(1,1,pch=21,cex=2,type='p',bg="white")
  text(ax,ay,1,cex=1)
  text(fx,fy,2,cex=1)
  text(dx,dy,3,cex=1)
  text(ex,ey,4,cex=1)
  text(1,1,5,cex=1)
  mtext(xlabz,side=1,line=0.7,cex=1)
  mtext(ylabz,side=2,line=0.5,cex=1)
}

boot_viol=function(my_table)
{    
    bootreg=c()
    for (b in 1:1000)
    {
      tabboot=my_table[,c("MP","NJ","ML","BI",'CNN')]
      boot=apply(apply(tabboot,2,sample,replace=T),2,mean)
      bootreg=rbind(bootreg,boot)
    } 
    return(bootreg)
}



tt=read.table("/Users/anton/Downloads/gapregions_1000.table")
names(tt)=c("ID","true_T","model","I","G","pars","nj","ml","bi","CNN_class")
tt$zone=rep(c( "EXP", "FA", "FAE" ,"FAT", "FE", "FEE","LONG","LONGOUT","LONGULTRA","SHORT","SHORTINT","SHORTOUT","SHORTULTRA"),each=3000)

#Read newicks
master=read.tree(text=as.character(tt$true_T))
pars=read.tree(text=as.character(tt$pars))
nj=read.tree(text=as.character(tt$nj))
ml=read.tree(text=as.character(tt$ml))
ba=read.tree(text=as.character(tt$bi))

#Calculate RF
RF=c()
for (i in 1:length(master))
{
  treerf=RF.dist(c(pars[[i]],nj[[i]],ml[[i]],ba[[i]]),master[[i]])
  RF=rbind(RF,treerf)
  print(i)
} 

tt=cbind(tt,data.frame(RF)+1)
tt$X5=as.numeric(tt$CNN_class==rep(rep(c(0,1,2),times=c(1000,1000,1000)),13))
tt[,c("X1","X2","X3","X4",'X5')]=ifelse(tt[,c("X1","X2","X3","X4",'X5')]!=1,0,1)

###Extract branch lengths
mm=unlist(strsplit(gsub("[:|,|[A-Z]|\\(|\\)|;| "," ",tt$true_T),split=" "))
br_l=data.frame(matrix(as.numeric(mm[mm!=""]),ncol=5,byrow=T),stringsAsFactors=F)
names(br_l)=c("A","B","inter","C","D")
tt=cbind(tt,br_l)
names(tt)[12:16]=c("MP","NJ","ML","BI","CNN")
tt$AB=tt$A+tt$B
tt$CD=tt$C+tt$D
tt$AC=tt$A+tt$C
tt$AD=tt$A+tt$D
tt$BC=tt$B+tt$C
tt$BC=tt$B+tt$C
tt$BD=tt$B+tt$D
tt$ABCD=tt$AB+tt$CD


###Accuracy
z_list=list()
for(z in unique(tt$zone))
{
    print(z)
    zacc=apply(tt[tt$zone==z,c("MP","NJ","ML","BI","CNN")],2,mean)
    zboot=boot_viol(tt[tt$zone==z,c("MP","NJ","ML","BI","CNN")])
    z_list[[z]]=list(z_acc=zacc,z_boot=zboot)   
}


#FA
bootreg=z_list[["FA"]][["z_boot"]]
vioplot(bootreg[,1],bootreg[,2],bootreg[,3],bootreg[,4],bootreg[,5],col="grey",pchMed=16,names=NA,ylim=c(0.53,1.02))
mtext("Accuracy",2,padj=-3,cex=0.7)
mtext(c("MP","NJ","ML","BI","CNN"),1,at=seq(1,5),padj=1,cex=0.6)
text(1:5,apply(bootreg,2,max)+0.017,z_list[["FA"]][["z_acc"]],digits=3),cex=0.7,xpd = TRUE,col="red")