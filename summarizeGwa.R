## summarize single SNP GWA
library(data.table)
library(scales)

## graphical parameters
cl<-1.4;ca<-1.1

## get het dat
het<-as.matrix(fread("HetAndP_color.txt",header=TRUE))
## SNPs with het were already dropped in the GWA, drop those first
hetx<-het[het[,1] != 1,]

## read gwa results for RG
resRG<-fread("output/o_RG.assoc.txt",header=TRUE)

## chrom and pos for RG
snpsRG<-as.matrix(fread("SNPsRG.txt",header=FALSE))

## -log10 P from LRG
lpRG<--log10(resRG$p_lrt)
L<-length(lpRG)
sigRG<-which(lpRG > 2 & hetx[,3] > -.2 & hetx[,3] < .2 & hetx[,2] > .025 & hetx[,2] < .975) ## get rid of extreme hets and rare

bnds<-c(which(snpsRG[-L,1] != snpsRG[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=snpsRG[,1],mean)

gws<- -log10(0.05/L)

pdf("GWA2019_RG_SSNP.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(sigRG,lpRG[sigRG],ylim=c(0,27),type='n',axes=FALSE,xlab="Chromosome",ylab="RG -Log10(P)",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,28,28),border=NA,col=alpha("gray20",.4))
}
points(sigRG,lpRG[sigRG],pch=19)
abline(h=gws,lty=3)
axis(1,at=mids,unique(snpsRG[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
box()
dev.off()

## now for Chr8 only
sigRG8<-which(lpRG > 2 & snpsRG[,1]==8)

pdf("GWA2019_RG_SSNP_Ch8.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(snpsRG[sigRG8,2],lpRG[sigRG8],ylim=c(0,27),pch=19,xlab="Position (bp)",ylab="RG -Log10(P)",cex.lab=cl)
dev.off()

## read gwa results for GB
resGB<-fread("output/o_GB.assoc.txt",header=TRUE)

## chrom and pos for RG
snpsGB<-as.matrix(fread("SNPsGB.txt",header=FALSE))

## -log10 P from LRG
lpGB<--log10(resGB$p_lrt)
L<-length(lpGB)
sigGB<-which(lpGB > 2 & hetx[,3] > -.2 & hetx[,3] < .2 & hetx[,2] > .025 & hetx[,2] < .975) ## get rid of extreme hets and rare

bnds<-c(which(snpsGB[-L,1] != snpsGB[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=snpsGB[,1],mean)

gws<- -log10(0.05/L)

pdf("GWA2019_GB_SSNP.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(sigGB,lpGB[sigGB],ylim=c(0,33.5),type='n',axes=FALSE,xlab="Chromosome",ylab="GB -Log10(P)",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,34,34),border=NA,col=alpha("gray20",.4))
}
points(sigGB,lpGB[sigGB],pch=19)
abline(h=gws,lty=3)
axis(1,at=mids,unique(snpsGB[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
box()
dev.off()

## now for Chr8 only
sigGB8<-which(lpGB > 2 & snpsGB[,1]==8)

pdf("GWA2019_GB_SSNP_Ch8.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(snpsGB[sigGB8,2],lpGB[sigGB8],ylim=c(0,33.5),pch=19,xlab="Position (bp)",ylab="GB -Log10(P)",cex.lab=cl)
dev.off()


## read gwa results for survival, both hosts and as binary
resSurv<-fread("output/o_SURV.assoc.txt",header=TRUE)

## chrom and pos for RG
snpsSurv<-as.matrix(fread("SNPsSurv.txt",header=FALSE))

## -log10 P from LRG
lpSurv<--log10(resSurv$p_lrt)
L<-length(lpSurv)
sigSurv<-which(lpSurv > 2) ## no filter yet

bnds<-c(which(snpsSurv[-L,1] != snpsSurv[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=snpsSurv[,1],mean)

gws<- -log10(0.05/L)

pdf("GWA2019_Surv_SSNP.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(sigSurv,lpSurv[sigSurv],ylim=c(0,10),type='n',axes=FALSE,xlab="Chromosome",ylab="Survival -Log10(P)",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,28,28),border=NA,col=alpha("gray20",.4))
}
points(sigSurv,lpSurv[sigSurv],pch=19)
abline(h=gws,lty=3)
axis(1,at=mids,unique(snpsSurv[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
box()
dev.off()

