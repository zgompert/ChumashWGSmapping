## looking at LD across GWA signal
library(data.table)
library(scales)

## 2019 experimental population
dat2019<-fread("hsub_mlpntest_EXP_tchum.geno",header=FALSE,sep=",")
dim(dat2019)
#[1] 5997988     432
g2019<-as.matrix(dat2019[,-c(1:3)])
snps2019<-read.table("SNPsRG_sub.txt",header=FALSE)

## ch8
ch8<-which(snps2019[,1]==8)
g_ch8<-g2019[ch8,]
snps_ch8<-snps2019[ch8,]
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
ph_2019<-read.table("ph_tchum2019.txt",header=FALSE)
nas<-which(is.na(ph_2019[,1])==TRUE)
ko<-kmeans(x=ph_2019[-nas,1:2],center=2,iter.max=20,nstart=50)
grp<-rep(NA,dim(ph_2019)[1])
grp[-nas]<-ko$cluster

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

pdf("LD_tchum_2019.pdf",width=9,height=9)
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
image(r2b,bre## ph with details
phdat<-read.table("sort_phfull_tchum2019.txt",header=TRUE,comment.char="?")

sv<-which(snps_ch8[,2] > 55000000 & snps_ch8[,2] < 63000000)
pc<-prcomp(t(g_ch8[sv,]),center=TRUE,scale=FALSE)
aks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)
title(main="Melanic",cex.main=1.3)
dev.off()

## ph with details
phdat<-read.table("sort_phfull_tchum2019.txt",header=TRUE,comment.char="?")

sv<-which(snps_ch8[,2] > 55000000 & snps_ch8[,2] < 63000000)
pc<-prcomp(t(g_ch8[sv,]),center=TRUE,scale=FALSE)
# summary(pc)
#Importance of components:
#                           PC1      PC2     PC3     PC4     PC5     PC6     PC7
#Standard deviation     14.6505 10.30433 8.91830 6.45832 6.37013 5.87209 5.27724
#Proportion of Variance  0.1521  0.07522 0.05635 0.02955 0.02875 0.02443 0.01973
#Cumulative Proportion   0.1521  0.22728 0.28362 0.31317 0.34192 0.36634 0.38607
par(mfrow=c(1,2))
plot(pc$x[,1],pc$x[,2],pch=19,col=phdat$hex,xlab="PC1 (15.2%)",ylab="PC2 (7.5%)",cex.lab=1.3)
plot(pc$x[,2],pc$x[,3],pch=19,col=phdat$hex,xlab="PC1 (7.5%)",ylab="PC3 (5.6%)",cex.lab=1.3)

ko<-kmeans(x=pc$x[,1:2],centers=6,nstart=50,iter.max=100)
plot(pc$x[,1],pc$x[,2],pch=19,col=phdat$hex,xlab="PC1 (15.2%)",ylab="PC2 (7.5%)",cex.lab=1.3)
text(pc$x[,1],pc$x[,2],ko$cluster)

## cluster based on a single run, homoz = 3, 5, 4
cgen<-matrix(0,ncol=3,nrow=429)
## "A" allele = melanic
cgen[ko$cluster==3,1]<-2
cgen[ko$cluster==6 | ko$cluster==1,1]<-1
## "B" allele = top of cluster triangle
cgen[ko$cluster==4,2]<-2
cgen[ko$cluster==6 | ko$cluster==2,2]<-1
## "C" allele = green
cgen[ko$cluster==5,3]<-2
cgen[ko$cluster==1 | ko$cluster==2,3]<-1
write.table(file="svGenotypes.txt",t(cgen),row.names=FALSE,col.names=FALSE,quote=FALSE)
save(list=ls(),file="ldk.rdat")


r2_1<-cor(t(g_ch8[sset,which(ko$cluster==1)]))^2
r2_1[upper.tri(r2_1)]<-NA
r2_1[diag(r2_1)]<-NA

image(r2_1,breaks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)

r2_2<-cor(t(g_ch8[sset,which(ko$cluster==2)]))^2
r2_2[upper.tri(r2_2)]<-NA
r2_2[diag(r2_2)]<-NA

image(r2_2,breaks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)

r2_3<-cor(t(g_ch8[sset,which(ko$cluster==3)]))^2
r2_3[upper.tri(r2_3)]<-NA
r2_3[diag(r2_3)]<-NA

image(r2_3,breaks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)

r2_4<-cor(t(g_ch8[sset,which(ko$cluster==4)]))^2
r2_4[upper.tri(r2_4)]<-NA
r2_4[diag(r2_4)]<-NA

## check out lack of variation, NA, for this one
image(r2_4,breaks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)

r2_5<-cor(t(g_ch8[sset,which(ko$cluster==5)]))^2
r2_5[upper.tri(r2_5)]<-NA
r2_5[diag(r2_5)]<-NA

## and this one
image(r2_5,breaks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)

r2_6<-cor(t(g_ch8[sset,which(ko$cluster==6)]))^2
r2_6[upper.tri(r2_6)]<-NA
r2_6[diag(r2_6)]<-NA

## sort of this one, bet less so
image(r2_6,breaks=c(-.1,.01,.05,.1,.5,1.1),col=alpha("firebrick",c(.04,.1,.3,.6,1)),axes=FALSE)
ub<-max(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lb<-min(which(snps_ch8[sset,2] > 55000000 & snps_ch8[sset,2] < 63000000))
lines(c(lb,ub)/len,c(lb,lb)/len)
lines(c(ub,ub)/len,c(lb,ub)/len)
