library(data.table) ## version 1.17.8
library(FactoMineR) ## version 2.16
library(mvcor)
library(scales) ## version
source("changeFunctionRV.R")
source("changeFunctionRVmod.R")


difmat<-vector("list",2)
difmat[[1]]<-matrix(NA,nrow=48,ncol=6) ## bounds and Z
difmat[[2]]<-matrix(NA,nrow=48,ncol=6) ## bounds and Z


difmat_m<-vector("list",2)
difmat_m[[1]]<-matrix(NA,nrow=48,ncol=6) ## bounds and Z
difmat_m[[2]]<-matrix(NA,nrow=48,ncol=6) ## bounds and Z


for(Z in 1:6){
cat("Started Z = ",Z,"\n")

############## green genome ##########

## read in the data
G<-fread("hsub_mlpntest_tchum_green.geno",header=FALSE,sep=",")
gmat<-as.matrix(G[,-c(1:3)])

hsubSNPs_green<-fread("hsubSNPs_green.txt",header=FALSE)
snp_green<-cbind(as.numeric(unlist(hsubSNPs_green[,1])),as.numeric(unlist(hsubSNPs_green[,2])))
ch<-snp_green[,1]
pos<-snp_green[,2]

ch8<-which(ch==4)
gmat_ch8<-gmat[ch8,]
pos_ch8<-pos[ch8]

trt<-read.table("sort_green_phfull_tchum2019.txt",header=TRUE,comment.char="%",sep=",")

regions<-read.csv("ChumashColorRegions-green_genome.csv")

nr<-dim(regions)[1]

## drop missing
nas<-which(is.na(trt$Survival)==TRUE)

res_green<-matrix(NA,nrow=nr,ncol=14)
res_green_m<-matrix(NA,nrow=nr,ncol=14)
for(rr in 1:nr){
	cat("Workgin on ",rr,"\n",sep="")
	rmat<-matrix(regions[rr,6:(6+regions$num_windows[rr]*2-1)],byrow=TRUE,ncol=2)
	res_green[rr,]<-LDChangeRV(regions=rmat,G=gmat_ch8[,-nas],pos=pos_ch8,
				   trt=trt[-nas,],npcs=Z)
	res_green_m[rr,]<-LDChangeRVm(regions=rmat,G=gmat_ch8[,-nas],pos=pos_ch8,
				   trt=trt[-nas,],npcs=Z)
}
cols<-c(3,5,7,10,12,14)

########### red genome #############

## read in the data
G<-fread("hsub_mlpntest_tchum_red.geno",header=FALSE,sep=",")
gmat<-as.matrix(G[,-c(1:3)])

hsubSNPs_red<-fread("hsubSNPs_red.txt",header=FALSE)
snp_red<-cbind(as.numeric(unlist(hsubSNPs_red[,1])),as.numeric(unlist(hsubSNPs_red[,2])))

ch<-snp_red[,1]
pos<-snp_red[,2]

ch8<-which(ch==4)
gmat_ch8<-gmat[ch8,]
pos_ch8<-pos[ch8]

trt<-read.table("sort_red_phfull_tchum2019.txt",header=TRUE,comment.char="%",sep=",")

regions<-read.csv("ChumashColorRegions-red_genome.csv")

nr<-dim(regions)[1]

## drop missing
nas<-which(is.na(trt$Survival)==TRUE)

res_red<-matrix(NA,nrow=nr,ncol=14)
res_red_m<-matrix(NA,nrow=nr,ncol=14)
for(rr in 1:nr){
	cat("Working on ",rr,"\n",sep="")
	rmat<-matrix(regions[rr,6:(6+regions$num_windows[rr]*2-1)],byrow=TRUE,ncol=2)
	res_red[rr,]<-LDChangeRV(regions=rmat,G=gmat_ch8[,-nas],pos=pos_ch8,
				 trt=trt[-nas,],npcs=Z)
	res_red_m[rr,]<-LDChangeRVm(regions=rmat,G=gmat_ch8[,-nas],pos=pos_ch8,
				 trt=trt[-nas,],npcs=Z)
}
cols<-c(3,5,7,10,12,14)



of<-paste("LDpvals_RV",Z,".pdf",sep="")

## combined plot
pdf(of,width=8,height=5)
par(mar=c(5,5,3,1))
cl<-1.4;ca<-1.1;ct<-1.4
res_comb<-rbind(res_green[,cols],res_red[,cols])
cs<-rep(c(alpha("forestgreen",.5),alpha("firebrick2",.5)),each=24)
boxplot(res_comb,outline=FALSE,ylab="P-value",xlab="Null model",
	names=c("random","linear RG","linear GB","random","linear RG","linear GB"),
	col=NULL,border="gray30",cex.lab=cl,cex.axis=ca)
points(jitter(rep(1:6,each=48)),res_comb,pch=19,col=rep(cs,6))
abline(v=3.5)
abline(h=.05,lty=3)
mtext("A/C",side=3,line=1,adj=.25,cex=ct)
mtext("MM",side=3,line=1,adj=.75,cex=ct)

res_comb<-rbind(res_green_m[,cols],res_red_m[,cols])
cs<-rep(c(alpha("forestgreen",.5),alpha("firebrick2",.5)),each=24)
boxplot(res_comb,outline=FALSE,ylab="P-value",xlab="Null model",
	names=c("random","linear RG","linear GB","random","linear RG","linear GB"),
	col=NULL,border="gray30",cex.lab=cl,cex.axis=ca)
points(jitter(rep(1:6,each=48)),res_comb,pch=19,col=rep(cs,6))
abline(v=3.5)
abline(h=.05,lty=3)
mtext("A/C",side=3,line=1,adj=.25,cex=ct)
mtext("MM",side=3,line=1,adj=.75,cex=ct)
dev.off()

difmat[[1]][,Z]<-c(res_green[,1]-res_green[,2],res_red[,1]-res_red[,2])
difmat[[2]][,Z]<-c(res_green[,8]-res_green[,9],res_red[,8]-res_red[,9])

difmat_m[[1]][,Z]<-c(res_green_m[,1]-res_green_m[,2],res_red_m[,1]-res_red_m[,2])
difmat_m[[2]][,Z]<-c(res_green_m[,8]-res_green_m[,9],res_red_m[,8]-res_red_m[,9])

rm("G","gmat")
save(list=ls(),file=paste("LDChange_RV",Z,".rdat",sep=""))
}



## combined plot
pdf("sfig_LDRV.pdf",width=10.5,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,4.5,2.5,1))
cl<-1.55;ca<-1.1;ct<-1.2
for(k in 1:6){
load(paste("LDChange_RV",k,".rdat",sep=""))
res_comb<-rbind(cbind(res_green[,1]-res_green[,2],res_green[,1]-res_green[,4],
		      res_green[,1]-res_green[,6],res_green[,8]-res_green[,9],
		      res_green[,8]-res_green[,11],res_green[,8]-res_green[,13]),
		cbind(res_red[,1]-res_red[,2],res_red[,1]-res_red[,4],
                      res_red[,1]-res_red[,6],res_red[,8]-res_red[,9],
                      res_red[,8]-res_red[,11],res_red[,8]-res_red[,13])
		)
ub<-max(as.numeric(res_comb))+0.01
lb<-min(as.numeric(res_comb))-0.01
res_p<-rbind(res_green[,cols],res_red[,cols])
cat(apply(res_p<0.05,2,sum),"\n")
cx<-c(21,19)
cs<-rep(c(alpha("forestgreen",.5),alpha("firebrick2",.5)),each=24)
boxplot(res_comb,outline=FALSE,ylab="Null-centered RV",xlab="Null model",
	names=c("random","linear RG","linear GB","random","linear RG","linear GB"),
	col=NULL,border="gray30",cex.lab=cl,cex.axis=ca,ylim=c(lb,ub))
points(jitter(rep(1:6,each=48)),res_comb,pch=cx[1+as.numeric(res_p <0.05)],col=rep(cs,6))
abline(v=3.5)
abline(h=0,lty=3)
mtext("A/C",side=3,line=-1.75,adj=.25,cex=ct)
mtext("MM",side=3,line=-1.75,adj=.75,cex=ct)
tout<-paste("(",letters[k],") RV correlations for k = ",k,sep="")
title(main=tout,cex.main=cl)
if(k==1){legend(.3,-.005,c("24_0163 (green genome)","24_0159 (red genome)"),pch=19,cex=1.2,col=c(alpha("forestgreen",.5),alpha("firebrick2",.5)),bty='n')}
}
dev.off()
