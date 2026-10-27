library(data.table)
library(scales)
source("changeFunction_between.R")

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
for(rr in 1:nr){
	cat("Workgin on ",rr,"\n",sep="")
	rmat<-matrix(regions[rr,6:(6+regions$num_windows[rr]*2-1)],byrow=TRUE,ncol=2)
	res_green[rr,]<-LDChangePC_between(regions=rmat,G=gmat_ch8[,-nas],pos=pos_ch8,trt=trt[-nas,])
}
cols<-c(3,5,7,10,12,14)

plot(jitter(rep(1:6,each=24)),res_green[,cols],pch=19)
abline(h=.05,lty=3)

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
for(rr in 1:nr){
	cat("Working on ",rr,"\n",sep="")
	rmat<-matrix(regions[rr,6:(6+regions$num_windows[rr]*2-1)],byrow=TRUE,ncol=2)
	res_red[rr,]<-LDChangePC_between(regions=rmat,G=gmat_ch8[,-nas],pos=pos_ch8,trt=trt[-nas,])
}
cols<-c(3,5,7,10,12,14)

plot(jitter(rep(1:6,each=24)),res_red[,cols],pch=19)
abline(h=.05,lty=3)

rm("G","gmat")
save(list=ls(),file="LDChange_between_pcs.rdat")

## combined plot
pdf("LDpvals_betweenpcs.pdf",width=8,height=5)
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
dev.off()
