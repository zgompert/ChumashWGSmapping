## looking at LD across GWA signal
library(data.table)
library(scales)

## GR8 experimental population
datGR8<-fread("hsub_mlpntest_GR8_tchum.geno",header=FALSE,sep=",")
dim(datGR8)
#[1] 12498367      460
gGR8<-as.matrix(datGR8[,-c(1:3)])
snpsGR8<-read.table("SNPsRG_GR8.txt",header=FALSE)

## ch8
ch8<-which(snpsGR8[,1]==8)
g_ch8<-gGR8[ch8,]
snps_ch8<-snpsGR8[ch8,]
p_ch8<-apply(g_ch8,1,mean)/2

## common snps > 50 mbp < 70 mbp
set<-which(p_ch8 > .05 & p_ch8 < .95 & snps_ch8[,2] > 50000000 & snps_ch8[,2] < 70000000)
## grab len=2000 equally spaced snps from this set
len<-2000
sset<-set[seq(1,length(set),length.out=len)]

r2<-cor(t(g_ch8[sset,]))^2
r2[upper.tri(r2)]<-NA
r2[diag(r2)]<-NA

image(r2,breaks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)

## now look at subsets by phenotype
ph_GR8<-read.table("sub_ph_tchumGR8.txt",header=FALSE)
nas<-which(is.na(ph_GR8[,1])==TRUE)
grp<-rep(NA,dim(ph_GR8)[1])
grp[-nas]<- (ph_GR8[-nas,1] > .05 & ph_GR8[-nas,2] < .21) +1
plot(ph_GR8[,1],ph_GR8[,2],col=c("green","brown")[grp],pch=19)

r2a<-cor(t(g_ch8[sset,which(grp==1)]))^2
r2a[upper.tri(r2a)]<-NA
r2a[diag(r2a)]<-NA

image(r2a,breaks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)

r2b<-cor(t(g_ch8[sset,which(grp==2)]))^2
r2b[upper.tri(r2b)]<-NA
r2b[diag(r2b)]<-NA

image(r2b,breaks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)

pdf("LD_tchum_GR8.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(1,3,3,1))

image(r2,breaks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)
title(main="Combined",cex.main=1.3)
image(r2a,breaks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)
title(main="Green",cex.main=1.3)
image(r2b,breaks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)
title(main="Melanic",cex.main=1.3)
dev.off()

## looking at putative SV vs phenotype
sv<-which(snps_ch8[,2] > 55000000 & snps_ch8[,2] < 63000000)
pc<-prcomp(t(g_ch8[sv,]),center=TRUE,scale=FALSE)
#Importance of components:
#                            PC1      PC2      PC3     PC4     PC5     PC6
#Standard deviation     16.09235 14.01399 12.82640 9.07003 8.89636 8.69101
#Proportion of Variance  0.09514  0.07216  0.06044 0.03022 0.02908 0.02775
#Cumulative Proportion   0.09514  0.16730  0.22774 0.25797 0.28705 0.31480

plot(pc$x[,1],pc$x[,2],pch=19)
plot(ph_GR8[,1],ph_GR8[,2])
cs<-rep(1,length(pc$x[,1]))
cs[ph_GR8[,1] > .05 & ph_GR8[,2] < .21]<-2
plot(ph_GR8[,1],ph_GR8[,2],col=c("green","brown")[cs])

plot(pc$x[,1],pc$x[,2],pch=19,col=c("green","brown")[cs])
