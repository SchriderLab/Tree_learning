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


color.bar <- function(lut, min, max=-min, nticks=3, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  #quartz()
  lines(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1,cex.axis=0.6)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
color.bar(inferno(1000),0,1)
quartz.save("heatlegend.pdf", type = "pdf",antialias=F,bg="white",dpi=800,pointsize=12)
dev.off()














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

plot_dens=function(zonetab,zone,y,xlim,n,title)
{
    ztab=zonetab[zonetab$zone==zone,]
    ztab$y_axis=y
    for (col in c("MP","NJ","ML","BI",'CNN'))
    {
      bind1=ztab[ztab[,col]!=1,c("inter","y_axis")]
      bind2=as.matrix(ztab[,c("inter","y_axis")])
      k1=kde2d(bind1[,1],bind1[,2], n=200,lims = c(c(0,xlim),c(min(bind2[,2]),max(bind2[,2]))))
      k2=kde2d(bind2[,1],bind2[,2], n=200,lims = c(c(0,xlim),c(min(bind2[,2]),max(bind2[,2]))))
      k1$z=k1$z/k2$z
      if (title == 1)
      {
          image(k1, col=inferno(2000),yaxs="i",main="",cex.axis=0.8)
          title(col,adj=0,line=0.5)
      }else{
          image(k1, col=inferno(2000),yaxs="i",main="",cex.axis=0.8)
      }    
      lines(bind1,type="p",pch=21,cex=0.3,bg="gray",lwd=0.3)

    }
}

plot_viol=function(z_list,zone,ylim1,ylim2,laby,topviol)
{    
    bootreg=z_list[[zone]][["z_boot"]]
    vioplot(bootreg[,1],bootreg[,2],bootreg[,3],bootreg[,4],bootreg[,5],col="grey",pchMed=16,names=NA,ylim=c(ylim1,ylim2))
    mtext(laby,2,padj=-3,cex=0.7)
    mtext(c("MP","NJ","ML","BI","CNN"),1,at=seq(1,5),padj=1,cex=0.6)
    text(1:5,apply(bootreg,2,max)+topviol,round(z_list[[zone]][["z_acc"]],digits=3),cex=0.7,xpd = TRUE,col="red")
}









table_prep=function(path_tab)
{    
    tt=read.table(path_tab)
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
    return(tt)
}

###Accuracy
acc_get=function(tt)
{
    z_list=list()
    for(z in unique(tt$zone))
    {
        print(z)
        zacc=apply(tt[tt$zone==z,c("MP","NJ","ML","BI","CNN")],2,mean)
        zboot=boot_viol(tt[tt$zone==z,c("MP","NJ","ML","BI","CNN")])
        z_list[[z]]=list(z_acc=zacc,z_boot=zboot)   
    }
    return(z_list)
}
"/Users/anton/Downloads/gapregions_1000.table"
###MAIN FIGS

###Zones GAPS
quartz(width=7.7, height=11)
par(mfcol=c(7,4),mar=c(2,3,2,1))
#FA
plot_topo(c(1,1,0,0,0),"a) Farris zone",expression('B'[5]),expression('B'[1+2]))
plot_viol(z_list,"FA",0.53,1.02,"Accuracy",0.023)
y=as.numeric(apply(tt[tt$zone=="FA",c("AB","CD","AC","AD","BC","BD")],1,max))
plot_dens(tt,"FA",y,0.05,200,1)
#FAT 
plot_topo(c(1,1,0,0,0),"b) Twisted Farris zone",expression('B'[5]),expression('B'[1+2]))
plot_viol(z_list,"FAT",0.52,1.02,"",0.023)
y=as.numeric(apply(tt[tt$zone=="FAT",c("AB","CD","AC","AD","BC","BD")],1,max))
plot_dens(tt,"FAT",y,0.05,200,0)
#FE
plot_topo(c(1,0,0,1,0),"c) Felsenstein zone",expression('B'[5]),expression('B'[1+3]))
plot_viol(z_list,"FE",0.1,0.8,"",0.03)
y=as.numeric(apply(tt[tt$zone=="FE",c("AB","CD","AC","AD","BC","BD")],1,max))
plot_dens(tt,"FE",y,0.05,200,0)
#SHORTINT
plot_topo(c(1,1,0,1,1),"d) Short internal branch",expression('B'[5]),expression('B'[1+2+3+4]))
plot_viol(z_list,"SHORTINT",0.45,0.71,"",0.013)
y=as.numeric(tt[tt$zone=="SHORTINT","ABCD"])
plot_dens(tt,"SHORTINT",y,0.05,200,0)
quartz.save("Bias_gap.jpeg", type = "jpeg",antialias=F,bg="white",dpi=400,pointsize=12)
dev.off()



###Zones NOGAPS
quartz(width=7.7, height=11)
par(mfcol=c(7,4),mar=c(2,3,2,1))
#FA
plot_topo(c(1,1,0,0,0),"a) Farris zone",expression('B'[5]),expression('B'[1+2]))
plot_viol(z_list,"FA",0.53,1.02,"Accuracy",0.023)
y=as.numeric(apply(tt[tt$zone=="FA",c("AB","CD","AC","AD","BC","BD")],1,max))
plot_dens(tt,"FA",y,0.05,200,1)
#FAT 
plot_topo(c(1,1,0,0,0),"b) Twisted Farris zone",expression('B'[5]),expression('B'[1+2]))
plot_viol(z_list,"FAT",0.52,1.02,"",0.023)
y=as.numeric(apply(tt[tt$zone=="FAT",c("AB","CD","AC","AD","BC","BD")],1,max))
plot_dens(tt,"FAT",y,0.05,200,0)
#FE
plot_topo(c(1,0,0,1,0),"c) Felsenstein zone",expression('B'[5]),expression('B'[1+3]))
plot_viol(z_list,"FE",0.1,0.8,"",0.03)
y=as.numeric(apply(tt[tt$zone=="FE",c("AB","CD","AC","AD","BC","BD")],1,max))
plot_dens(tt,"FE",y,0.05,200,0)
#SHORTINT
plot_topo(c(1,1,0,1,1),"d) Short internal branch",expression('B'[5]),expression('B'[1+2+3+4]))
plot_viol(z_list,"SHORTINT",0.45,0.71,"",0.013)
y=as.numeric(tt[tt$zone=="SHORTINT","ABCD"])
plot_dens(tt,"SHORTINT",y,0.05,200,0)
quartz.save("Bias_nogap.jpeg", type = "jpeg",antialias=F,bg="white",dpi=400,pointsize=12)
dev.off()



##SUPPL FIGS GAP / NOGAP
quartz(width=10.6, height=10.1)
par(mfrow=c(8,7),mar=c(2,2,2,2))
#FAE
plot_topo(c(1,1,1,0,0),"a) Extended \nFarris zone",expression('B'[5]),expression('B'[1+2]))
plot_viol(z_list,"FAE",0.8,1.02,"Accuracy",0.023)
y=as.numeric(apply(tt[tt$zone=="FAE",c("AB","CD","AC","AD","BC","BD")],1,max))
plot_dens(tt,"FAE",y,0.5,200,1)

#FEE
plot_topo(c(1,0,1,1,0),"b) Extended \nFelsenstein zone",expression('B'[5]),expression('B'[1+3]))
plot_viol(z_list,"FEE",0.7,1.02,"Accuracy",0.023)
y=as.numeric(apply(tt[tt$zone=="FEE",c("AB","CD","AC","AD","BC","BD")],1,max))
plot_dens(tt,"FEE",y,0.5,200,1)

#LONG
plot_topo(c(1,1,1,1,1),"c) Long branches",expression('B'[5]),expression('B'[1+2+3+4]))
plot_viol(z_list,"LONG",0.7,1.02,"Accuracy",0.023)
y=as.numeric(tt[tt$zone=="LONG","ABCD"])
plot_dens(tt,"LONG",y,0.5,200,1)

#LONGULTRA
plot_topo(c(1,1,1,1,1),"d) Extra-long branches",expression('B'[5]),expression('B'[1+2+3+4]))
plot_viol(z_list,"LONGULTRA",0.6,1,"Accuracy",0.023)
y=as.numeric(tt[tt$zone=="LONGULTRA","ABCD"])
plot_dens(tt,"LONGULTRA",y,1,200,1)

#LONGOUT
plot_topo(c(1,0,1,0,0),"e) Single long branch",expression('B'[5]),expression('B'[1]))
plot_viol(z_list,"LONGOUT",0.8,1.02,"Accuracy",0.023)
y=as.numeric(apply(tt[tt$zone=="LONGOUT",c("A","B","C","D")],1,max))
plot_dens(tt,"LONGOUT",y,0.5,200,1)

#SHORT
plot_topo(c(0,0,1,0,0),"f) Short branches",expression('B'[5]),expression('B'[1+2+3+4]))
plot_viol(z_list,"SHORT",0.96,1,"Accuracy",0.004)
y=as.numeric(tt[tt$zone=="SHORT","ABCD"])
plot_dens(tt,"SHORT",y,0.5,200,1)

#SHORTULTRA
plot_topo(c(0,0,0,0,0),"g) Extra-short branches",expression('B'[5]),expression('B'[1+2+3+4]))
plot_viol(z_list,"SHORTULTRA",0.6,1,"Accuracy",0.023)
y=as.numeric(tt[tt$zone=="SHORTULTRA","ABCD"])
plot_dens(tt,"SHORTULTRA",y,0.01,200,1)

#SHORTOUT
plot_topo(c(0,1,1,1,1),"e) Single short branch",expression('B'[5]),expression('B'[1]))
plot_viol(z_list,"SHORTOUT",0.8,1.02,"Accuracy",0.023)
y=as.numeric(apply(tt[tt$zone=="SHORTOUT",c("A","B","C","D","BC")],1,min))
plot_dens(tt,"SHORTOUT",y,0.5,200,1)
quartz.save("Suppl_Bias_gap.jpeg", type = "jpeg",antialias=F,bg="white",dpi=400,pointsize=12)
dev.off()
