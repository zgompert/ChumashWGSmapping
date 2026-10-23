## summarize single SNP GWA
library(data.table)
library(scales)

## graphical parameters
cl<-1.5;ca<-1.1;cm<-1.5

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

## now for Chr8 only
sigRG8<-which(lpRG > 2 & snpsRG[,1]==8)

## read gwa results for GB
resGB<-fread("output/o_GB.assoc.txt",header=TRUE)

## chrom and pos for RG
snpsGB<-as.matrix(fread("SNPsGB.txt",header=FALSE))

## -log10 P from LRG
lpGB<--log10(resGB$p_lrt)
L<-length(lpGB)
sigGB<-which(lpGB > 2 & hetx[,3] > -.2 & hetx[,3] < .2 & hetx[,2] > .025 & hetx[,2] < .975) ## get rid of extreme hets and rare

########## GR8
## get het dat
het_GR8<-as.matrix(fread("HetAndP_color_GR8.txt",header=TRUE))
## SNPs with het devs dropped already
hetx_GR8<-het_GR8[het_GR8[,3] > -.2 & het_GR8[,3] < .2,]

resRG_GR8<-fread("output/o_GR8_RG.assoc.txt",header=TRUE)
## chrom and pos for RG, note extreme het removed
snpsRG_GR8<-as.matrix(fread("SNPsRG_GR8.txt",header=FALSE))

maxSNPs_GR8<-tapply(INDEX=snpsRG_GR8[,1],X=snpsRG_GR8[,2],max)

lpRG_GR8<--log10(resRG_GR8$p_lrt)
L_GR8<-length(lpRG_GR8)
sigRG_GR8<-which(lpRG_GR8 > 2 & hetx_GR8[,3] > -.2 & hetx_GR8[,3] < .2 & hetx_GR8[,2] > .025
 & hetx_GR8[,2] < .975)

resGB_GR8<-fread("output/o_GR8_GB.assoc.txt",header=TRUE)
## chrom and pos for RG, note extreme het removed
snpsGB_GR8<-as.matrix(fread("SNPsRG_GR8.txt",header=FALSE))

lpGB_GR8<--log10(resGB_GR8$p_lrt)
L_GR8<-length(lpGB_GR8)
sigGB_GR8<-which(lpGB_GR8 > 2 & hetx_GR8[,3] > -.2 & hetx_GR8[,3] < .2 & hetx_GR8[,2] > .025 & hetx_GR8[,2] < .975)


### figure 3 for manuscript #####
pdf("fig_gwaZoom_SSNP.pdf",width=8,height=9.5)
layout(matrix(1:5,nrow=5,ncol=1),widths=8,heights=c(1.5,2,2,2,2))

par(mar=c(4.5,5.5,2.5,1.5))

ub<-max(snpsRG[snpsRG[,1]==8,2])
plot(c(1,ub),c(.3,.3),type='l',lwd=5,axes=FALSE,xlab="Position (bp)",ylab="",ylim=c(0,1),cex.lab=cl)
polygon(c(5.2e7,6.6e7,6.6e7,5.2e7),c(.0,.0,.7,.7),col=alpha("tan",.2))
lines(c(0.2e7,1.2e7),c(.9,.9),lwd=1.8)
lines(c(0.2e7,0.2e7),c(.85,.95),lwd=1.8)
lines(c(1.2e7,1.2e7),c(.85,.95),lwd=1.8)
text(0.75e7,.7,"10 Mb")
title(main="(a) Chromosome 8",cex.main=cm)

zRG<-which(hetx[,3] > -.2 & hetx[,3] < .2 & hetx[,2] > .025 & hetx[,2] < .975 & snpsRG[,1]==8 & snpsRG[,2]>52000000 & snpsRG[,2]<66000000)

plot(snpsRG[zRG,2],lpRG[zRG],pch=21,xlab="Position (bp)",ylab="-log10(P)",cex.lab=cl,col=NA,bg=alpha("firebrick",.5))
title(main="(b) Chromosome 8 GWA RG for HF5",cex.main=cm)
lines(c(5.2e7,5.3e7),c(20,20),lwd=1.8)
lines(c(5.2e7,5.2e7),c(19,21),lwd=1.8)
lines(c(5.3e7,5.3e7),c(19,21),lwd=1.8)
text(5.25e7,16,"1 Mb")


zRG_GR8<-which(hetx_GR8[,3] > -.2 & hetx_GR8[,3] < .2 & hetx_GR8[,2] > .025 & hetx_GR8[,2] < .975 & snpsRG_GR8[,1]==8 & snpsRG_GR8[,2]>52000000 & snpsRG_GR8[,2]<66000000)

plot(snpsRG_GR8[zRG_GR8,2],lpRG_GR8[zRG_GR8],pch=21,xlab="Position (bp)",ylab="-log10(P)",cex.lab=cl,col=NA,bg=alpha("firebrick",.5))
title(main="(c) Chromosome 8 GWA RG for GR8",cex.main=cm)

zGB<-which(hetx[,3] > -.2 & hetx[,3] < .2 & hetx[,2] > .025 & hetx[,2] < .975 & snpsGB[,1]==8 & snpsGB[,2]>52000000 & snpsGB[,2]<66000000)

plot(snpsGB[zGB,2],lpGB[zGB],pch=21,xlab="Position (bp)",ylab="-log10(P)",cex.lab=cl,col=NA,bg=alpha("royalblue3",.5))
title(main="(d) Chromosome 8 GWA GB for HF5",cex.main=cm)

zGB_GR8<-which(hetx_GR8[,3] > -.2 & hetx_GR8[,3] < .2 & hetx_GR8[,2] > .025 & hetx_GR8[,2] < .975 & snpsGB_GR8[,1]==8 & snpsGB_GR8[,2]>52000000 & snpsGB_GR8[,2]<66000000)

plot(snpsGB_GR8[zGB_GR8,2],lpGB_GR8[zGB_GR8],pch=21,xlab="Position (bp)",ylab="-log10(P)",cex.lab=cl,col=NA,bg=alpha("royalblue3",.5))
title(main="(e) Chromosome 8 GWA GB for GR8",cex.main=cm)

dev.off()


