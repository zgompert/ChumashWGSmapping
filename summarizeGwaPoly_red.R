## summarize BSLMM GWA, red genome
library(data.table)
library(scales)

## graphical parameters
cl<-1.4;ca<-1.1

## read gwa results for RG
resRG<-fread("output/pips_o_poly_bslmm_red_ph1_ch0.param.txt",header=FALSE)


pips<-as.numeric(unlist(resRG[,5]))
L<-length(pips)
sigRG<-which(pips > 0.00001) 

bnds<-c(which(resRG[-L,1] != resRG[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=resRG[,1],mean)

## by scaffold number 4 = ch8
pdf("GWA_red_RG_poly.pdf",width=8,height=4.5)
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

pdf("GWA_red_RG_poly_Ch8.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(unlist(resRG[sigRG8,2]),pips[sigRG8],ylim=c(0,1),pch=19,xlab="Position (bp)",ylab="Posterior inclusion probability",cex.lab=cl)
dev.off()


sum(pips[resRG[,1]==4])
#[1] 2.317189 
1-prod(1-pips[resRG[,1]==4])
#[1] 0.9085576
sum(pips)
#[1] 9.944713
sum(pips[pips > 0.05])
#[1] 1.680172
sum(pips[resRG[,1]==4 & pips > 0.05])
#[1] 1.014856

## regions from single res
reg_red<-read.csv("ChumashColorRegions-red_genome.csv")
r1lb<-reg_red[1,6]
r1ub<-reg_red[1,7]
r2lb<-reg_red[1,8]
r2ub<-reg_red[1,9]

sum(pips[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r1ub])
#[1] 0.6506288
sum(pips[resRG[,1]==4 & resRG[,2] > r2lb & resRG[,2] < r2ub])
#[1] 0.4121167
sum(pips[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r2ub])
#[1] 1.062892
pipsRG<-pips


## read gwa results for GB
resGB<-fread("output/pips_o_poly_bslmm_red_ph2_ch0.param.txt",header=FALSE)

pips<-as.numeric(unlist(resGB[,5]))
L<-length(pips)
sigGB<-which(pips > 0.00001) 
bnds<-c(which(resGB[-L,1] != resGB[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=resGB[,1],mean)


pdf("GWA_red_GB_poly.pdf",width=8,height=4.5)
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

pdf("GWA_red_GB_poly_Ch8.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(unlist(resGB[sigGB8,2]),pips[sigGB8],ylim=c(0,1),pch=19,xlab="Position (bp)",ylab="Posterior inclusion probability",cex.lab=cl)
dev.off()

sum(pips[resGB[,1]==4])
#[1] 3.160622 
1-prod(1-pips[resGB[,1]==4])
#[1] 0.9695394
sum(pips)
#[1] 11.0378
sum(pips[pips > 0.05])
#[1] 1.501023
sum(pips[resGB[,1]==4 & pips > 0.05])
#[1] 1.501023

## regions from single res
reg_red<-read.csv("ChumashColorRegions-red_genome.csv")
r1lb<-reg_red[1,6]
r1ub<-reg_red[1,7]
r2lb<-reg_red[1,8]
r2ub<-reg_red[1,9]

sum(pips[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r1ub])
#[1] 1.531005
sum(pips[resGB[,1]==4 & resGB[,2] > r2lb & resGB[,2] < r2ub])
#[1] 0.5784108
sum(pips[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r2ub])
#2.109506
pipsGB<-pips

## number of QTL total
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r1ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r1ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r1ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r1ub])
#[1] 2.181588 
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r2lb & resGB[,2] < r2ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r2lb & resRG[,2] < r2ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r2lb & resGB[,2] < r2ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r2lb & resRG[,2] < r2ub])
#[1] 0.9905177
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] <r2ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r2ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r2ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r2ub])
#[1] 3.172344
