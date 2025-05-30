## summarize single SNP GWA
library(data.table)
library(scales)

## graphical parameters
cl<-1.4;ca<-1.1


## read gwa results for RG
resRG<-fread("output/o_GR8_RG.assoc.txt",header=TRUE)

## chrom and pos for RG, note extreme het removed
snpsRG<-as.matrix(fread("SNPsRG_GR8.txt",header=FALSE))

## -log10 P from LRG
lpRG<--log10(resRG$p_lrt)
L<-length(lpRG)
sigRG<-which(lpRG > 2)

bnds<-c(which(snpsRG[-L,1] != snpsRG[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=snpsRG[,1],mean)

gws<- -log10(0.05/L)

pdf("GWAGR8_RG_SSNP.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(sigRG,lpRG[sigRG],ylim=c(0,33),type='n',axes=FALSE,xlab="Chromosome",ylab="RG -Log10(P)",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,35,35),border=NA,col=alpha("gray20",.4))
}
points(sigRG,lpRG[sigRG],pch=19)
abline(h=gws,lty=3)
axis(1,at=mids,unique(snpsRG[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
box()
dev.off()

## now for Chr8 only
sigRG8<-which(lpRG > 2 & snpsRG[,1]==8)

pdf("GWAGR_RG_SSNP_Ch8.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(snpsRG[sigRG8,2],lpRG[sigRG8],ylim=c(0,33),pch=19,xlab="Position (bp)",ylab="RG -Log10(P)",cex.lab=cl)
dev.off()

## read gwa results for GB
resGB<-fread("output/o_GR8_GB.assoc.txt",header=TRUE)

## chrom and pos for GB, same as RG
snpsGB<-snpsRG

## -log10 P from LGB
lpGB<--log10(resGB$p_lrt)
L<-length(lpGB)
sigGB<-which(lpGB > 2)

bnds<-c(which(snpsGB[-L,1] != snpsGB[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=snpsGB[,1],mean)

gws<- -log10(0.05/L)

pdf("GWAGR8_GB_SSNP.pdf",width=8,height=4.5)
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

pdf("GWAGR8_GB_SSNP_Ch8.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(snpsGB[sigGB8,2],lpGB[sigGB8],ylim=c(0,33.5),pch=19,xlab="Position (bp)",ylab="GB -Log10(P)",cex.lab=cl)
dev.off()

