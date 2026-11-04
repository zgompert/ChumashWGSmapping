library(lostruct)
packageVersion("lostruct")
#[1] ‘0.0.0.9000’

library(data.table)
setDTthreads(8)

### Green
G<-as.matrix(fread("ch8_mlpntest_tchum_green.csv",header=FALSE,sep=","))
dim(G)
#[1] 3278110     455

## full run
pcs<-eigen_windows(G,win=1000,k=2)
pcd<-pc_dist(pcs,npc=2)
fit2d <- cmdscale(pcd, eig=TRUE, k=2 )
plot( fit2d$points, xlab="Coordinate 1", ylab="Coordinate 2", col=rainbow(1.2*nrow(pcd)) )

## I think axis 1 is the big inversion and axis 2 is the color-associated region
plot(fit2d$points[,1])
plot(fit2d$points[,2])

## read in snp positions
spos<-scan("ch8_snppos_tchum_green.txt")

L<-3278*1000 ## windows by size
pmat<-matrix(spos[1:L],nrow=3278,ncol=1000,byrow=TRUE)
mpos<-apply(pmat,1,median)
plot(mpos,fit2d$points[,1])
plot(mpos,fit2d$points[,1],pch=19)
plot(mpos,fit2d$points[,2],pch=19)
plot(mpos,fit2d$points[,2],pch=19,xlim=c(5.8e7,6.6e7))

## just looking at LD
gwar<-which(spos >= 5.8e7 & spos <=6.6e7)
## then thin it, too many snps
gwar<-gwar[seq(1,length(gwar),by=10)]
## PCA by snps not inds
pcsnp<-prcomp(G[gwar,],center=TRUE,scale=FALSE)

outs<-which(abs(fit2d$points[,2]) > mean(fit2d$points[,2])+3*sd(fit2d$points[,2]))
cs<-rep(1,3278)
cs[outs]<-2
table(cs)
## maybe a tad conservative, mostly contiguous
plot(mpos,fit2d$points[,2],pch=19,col=c("gray","firebrick")[cs])
min(pmat[2429,])
#[1] 59926636
max(pmat[2470,])
#[1] 64870040



save(list=ls(),file="locpca_green.rdat")

rm(list=ls())
### Red
G<-as.matrix(fread("ch8_mlpntest_tchum_red.csv",header=FALSE,sep=","))
dim(G)
#[1] 3293533      455

## full run
pcs<-eigen_windows(G,win=1000,k=2)
pcd<-pc_dist(pcs,npc=2)
fit2d <- cmdscale(pcd, eig=TRUE, k=2 )
plot( fit2d$points, xlab="Coordinate 1", ylab="Coordinate 2", col=rainbow(1.2*nrow(pcd)) )

## I think axis 1 is the big inversion and axis 2 is the color-associated region
plot(fit2d$points[,1])
plot(fit2d$points[,2])

cl<-1.3;cm<-1.3
pdf("winpca_green.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4.5,4.5,2,1))
plot(mpos,-1*fit2d$points[,2],pch=19,ylim=c(-.6,.8),type='n',xlab="Position",ylab="Coordinate",cex.lab=cl)
polygon(c(6e7,6.5e7,6.5e7,6e7),c(-1,-1,1,1),border=NA,col=alpha("brown",.2))
points(mpos,-1*fit2d$points[,2],pch=19,col=alpha("firebrick",.7))
points(mpos,fit2d$points[,1],pch=19,col=alpha("gray20",.7))
legend(.1e7,.7,c("coordinate 1","coordinate 2"),pch=19,col=alpha(c("gray20","firebrick"),.7),bty='n')
title("(a) Chromosome 8 for 24_0163",cex.main=cm)

plot(mpos,-1*fit2d$points[,2],pch=19,ylim=c(-.6,.8),type='n',xlim=c(5.5e7,7e7),,xlab="Position",ylab="Coordinate",cex.lab=cl)
polygon(c(6e7,6.5e7,6.5e7,6e7),c(-1,-1,1,1),border=NA,col=alpha("brown",.2))
points(mpos,-1*fit2d$points[,2],pch=19,col=alpha("firebrick",.7))
points(mpos,fit2d$points[,1],pch=19,col=alpha("gray20",.7))
title("(b) Zoomed in view for 24_0163",cex.main=cm)
dev.off()

## read in snp positions
spos<-scan("ch8_snppos_tchum_red.txt")

L<-3293*1000 ## windows by size
pmat<-matrix(spos[1:L],nrow=3293,ncol=1000,byrow=TRUE)
mpos<-apply(pmat,1,median)
plot(mpos,fit2d$points[,1])
plot(mpos,fit2d$points[,1],pch=19)
plot(mpos,fit2d$points[,2],pch=19)
plot(mpos,fit2d$points[,2],pch=19,xlim=c(5.8e7,6.6e7))
plot(mpos,fit2d$points[,1],pch=19)



## define local pca region
outs<-which(abs(fit2d$points[,2]) > mean(fit2d$points[,2])+3*sd(fit2d$points[,2]))
cs<-rep(1,3293)
cs[outs]<-2
table(cs)
## maybe a tad conservative, but okay, use continuous set not extra peak
plot(mpos,fit2d$points[,2],pch=19,col=c("gray","firebrick")[cs])
min(pmat[2438,])
#[1] 60656171
max(pmat[2481,])
#[1] 63854024

cl<-1.3;cm<-1.3
pdf("winpca_red.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4.5,4.5,2,1))
plot(mpos,-1*fit2d$points[,2],pch=19,ylim=c(-.6,.8),type='n',xlab="Position",ylab="Coordinate",cex.lab=cl)
polygon(c(6.05e7,6.4e7,6.4e7,6.05e7),c(-1,-1,1,1),border=NA,col=alpha("brown",.2))
points(mpos,-1*fit2d$points[,2],pch=19,col=alpha("firebrick",.7))
points(mpos,fit2d$points[,1],pch=19,col=alpha("gray20",.7))
legend(.1e7,.7,c("coordinate 1","coordinate 2"),pch=19,col=alpha(c("gray20","firebrick"),.7),bty='n')
title("(c) Chromosome 8 for 24_0159",cex.main=cm)

plot(mpos,-1*fit2d$points[,2],pch=19,ylim=c(-.6,.8),type='n',xlim=c(5.5e7,7e7),,xlab="Position",ylab="Coordinate",cex.lab=cl)
polygon(c(6.05e7,6.4e7,6.4e7,6.05e7),c(-1,-1,1,1),border=NA,col=alpha("brown",.2))
points(mpos,-1*fit2d$points[,2],pch=19,col=alpha("firebrick",.7))
points(mpos,fit2d$points[,1],pch=19,col=alpha("gray20",.7))
title("(d) Zoomed in view for 24_0159",cex.main=cm)
dev.off()

save(list=ls(),file="locpca_red.rdat")
