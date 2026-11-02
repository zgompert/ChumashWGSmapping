## summarize BSLMM GWA, green genome
library(data.table)
library(scales)

## graphical parameters
cl<-1.4;ca<-1.1

## read gwa results for RG
resRG<-fread("output/pips_o_poly_bslmm_green_ph1_ch0.param.txt",header=FALSE)


pips<-as.numeric(unlist(resRG[,5]))
L<-length(pips)
sigRG<-which(pips > 0.00001) 

bnds<-c(which(resRG[-L,1] != resRG[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=resRG[,1],mean)

## by scaffold number 4 = ch8
pdf("GWA_green_RG_poly.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(sigRG,pips[sigRG],ylim=c(0,1),type='n',axes=FALSE,xlab="Chromosome",ylab="Posterior inclusion probability",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,28,28),border=NA,col=alpha("gray20",.4))
}
points(sigRG,pips[sigRG],pch=19)
axis(1,at=mids,unique(unlist(resRG[,1])),cex.axis=ca)
axis(2,cex.axis=ca)
box()
dev.off()

## now for Chr8 only
sigRG8<-which(pips > 0.00001 & resRG[,1]==4)

pdf("GWA_green_RG_poly_Ch8.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(unlist(resRG[sigRG8,2]),pips[sigRG8],ylim=c(0,1),pch=19,xlab="Position (bp)",ylab="Posterior inclusion probability",cex.lab=cl)
dev.off()


sum(pips[resRG[,1]==4])
#[1] 3.24537
1-prod(1-pips[resRG[,1]==4])
#[1] 0.9746839
sum(pips)
#[1]  11.45243
sum(pips[pips > 0.05])
#[1] 2.63675
sum(pips[resRG[,1]==4 & pips > 0.05])
#[1]1.845746

## regions from single res
reg_green<-read.csv("ChumashColorRegions-green_genome.csv")
r1lb<-reg_green[1,6]
r1ub<-reg_green[1,7]
r2lb<-reg_green[1,8]
r2ub<-reg_green[1,9]
r3lb<-reg_green[1,10]
r3ub<-reg_green[1,11]

sum(pips[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r1ub])
#[1] 1.648558 
sum(pips[resRG[,1]==4 & resRG[,2] > r2lb & resRG[,2] < r2ub])
#[1] 0.3680095
sum(pips[resRG[,1]==4 & resRG[,2] > r3lb & resRG[,2] < r3ub])
#[1] 0.08634567
sum(pips[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r3ub])
#[1] 2.103527
pipsRG<-pips

nr<-dim(reg_green)[1]
pipRG_region<-rep(NA,nr)
for(k in 1:nr){
	sp<-rep(NA,reg_green$num_windows[k])
	for(a in 1:reg_green$num_windows[k]){
		b<-(a-1)*2
		lb<-reg_green[k,6+b]
		ub<-reg_green[k,7+b]
		sp[a]<-sum(pips[resRG[,1]==4 & resRG[,2] > lb & resRG[,2] < ub])
	}
	pipRG_region[k]<-sum(sp)
}

## read gwa results for GB
resGB<-fread("output/pips_o_poly_bslmm_green_ph2_ch0.param.txt",header=FALSE)

pips<-as.numeric(unlist(resGB[,5]))
L<-length(pips)
sigGB<-which(pips > 0.00001) 
bnds<-c(which(resGB[-L,1] != resGB[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=resGB[,1],mean)


pdf("GWA_green_GB_poly.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(sigGB,pips[sigGB],ylim=c(0,1),type='n',axes=FALSE,xlab="Chromosome",ylab="Posterior inclusion probability",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),border=NA,col=alpha("gray20",.4))
}
points(sigGB,pips[sigGB],pch=19)
axis(1,at=mids,unique(unlist(resGB[,1])),cex.axis=ca)
axis(2,cex.axis=ca)
box()
dev.off()

## now for Chr8 only
sigGB8<-which(pips > 0.00001 & resGB[,1]==4)

pdf("GWA_green_GB_poly_Ch8.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(unlist(resGB[sigGB8,2]),pips[sigGB8],ylim=c(0,1),pch=19,xlab="Position (bp)",ylab="Posterior inclusion probability",cex.lab=cl)
dev.off()

sum(pips[resGB[,1]==4])
#[1]  3.084562
1-prod(1-pips[resGB[,1]==4])
#[1] 0.9674149
sum(pips)
#[1] 13.13639
sum(pips[pips > 0.05])
#[1] 1.411071
sum(pips[resGB[,1]==4 & pips > 0.05])
#[1] 1.411071

## regions from single res
reg_green<-read.csv("ChumashColorRegions-green_genome.csv")
r1lb<-reg_green[1,6]
r1ub<-reg_green[1,7]
r2lb<-reg_green[1,8]
r2ub<-reg_green[1,9]
r3lb<-reg_green[1,10]
r3ub<-reg_green[1,11]

sum(pips[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r1ub])
#[1] 1.23329
sum(pips[resGB[,1]==4 & resGB[,2] > r2lb & resGB[,2] < r2ub])
#[1] 0.8682238
sum(pips[resGB[,1]==4 & resGB[,2] > r3lb & resGB[,2] < r3ub])
#[1] 0.06072583

sum(pips[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r3ub])
pipsGB<-pips

nr<-dim(reg_green)[1]
pipGB_region<-rep(NA,nr)
for(k in 1:nr){
        sp<-rep(NA,reg_green$num_windows[k])
        for(a in 1:reg_green$num_windows[k]){
                b<-(a-1)*2
                lb<-reg_green[k,6+b]
                ub<-reg_green[k,7+b]
                sp[a]<-sum(pips[resGB[,1]==4 & resGB[,2] > lb & resGB[,2] < ub])
        }
        pipGB_region[k]<-sum(sp)
}



## number of QTL total
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r1ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r1ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r1ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r1ub])
#[1] 2.881533 
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r2lb & resGB[,2] < r2ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r2lb & resRG[,2] < r2ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r2lb & resGB[,2] < r2ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r2lb & resRG[,2] < r2ub])
#[1] 1.235952
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r3lb & resGB[,2] < r3ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r3lb & resRG[,2] < r3ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r3lb & resGB[,2] < r3ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r3lb & resRG[,2] < r3ub])
#[1] 0.1470668
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r3ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r3ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r3ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r3ub])
#[1] 4.273493


summary(pipRG_region)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.174   1.865   2.042   1.936   2.102   2.103 
summary(pipGB_region)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.428   2.084   2.152   2.073   2.162   2.162

pdf("sfig_qtn_green.pdf",width=5,height=5)
par(mar=c(4.5,5.5,2.5,1.5))
plot(pipRG_region,pipGB_region,pch=19,col=alpha("gray30",.5),xlab="Number of QTN RG",ylab="Number of QTN GB",cex.lab=1.3,cex.axis=1)
title(main="(a) Green genome (24_0163H1)",cex.main=1.3)
dev.off()
