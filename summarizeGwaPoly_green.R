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

## number of QTL total
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r1ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r1ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r1ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r1ub])
#[1] 2.881533 
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r2lb & resGB[,2] < r2ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r2lb & resRG[,2] < r2ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r2lb & resGB[,2] < r2ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r2lb & resRG[,2] < r2ub])
#[1] 1.235952
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r3lb & resGB[,2] < r3ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r3lb & resRG[,2] < r3ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r3lb & resGB[,2] < r3ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r3lb & resRG[,2] < r3ub])
#[1] 0.1470668
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r3ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r3ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r3ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r3ub])
#[1] 4.273493

### ANNOTATION NOT DONE #####

## check annotation for SNPs with pips > X
gtf<-read.table("/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/t_chumash/braker/braker.gtf",header=FALSE,fill=TRUE)

genes<-which(gtf$V3=="gene")
CDS<-which(gtf$V3=="CDS")

chromtab<-read.table("ChromLookUp.txt",header=FALSE)


## RG
resRG<-as.matrix(fread("SNPsRG_sub.txt",header=FALSE))
pips<-as.numeric(resRG$V2)
pip10<-which(pips > .1)
resRG[pip10,]
#      V1        V2
# [1,]  1 731095677
# [2,]  1 732566323
# [3,]  1 771077484
# [4,]  1 781206805
# [5,]  7  70200596
# [6,]  7  75600304
# [7,]  8  59330221
# [8,]  8  59587922
# [9,]  8  59630464
#[10,]  8 136007084
#[11,] 13  32771453


RGinGene<-rep(NA,length(pip10))
for(i in 1:length(pip10)){
	ch<-chromtab[chromtab[,1]==resRG[pip10][i],2]
	xx<-which(gtf[genes,1]==ch & resRG[pip10,2][i] > gtf[genes,4] & resRG[pip10,2][i] < gtf[genes,5])
	if(length(xx)>=1){
		RGinGene[i]<-gtf[genes[xx],9]
	}
}
## SNP 2 (1, 732566323) = g15738 = Proteasome-substrate-size regulator, mid region
## SNP 11 (13, 32771453) = g3312 = PH domain

## GB
resGB<-resRG
pips<-as.numeric(resGB$V2)
pip10<-which(pips > .1)
resGB[pip10,]
#      V1        V2
#[1,]  8 56786125
#[2,]  8 59330149
#[3,]  8 59592931
#[4,]  8 59695850


GBinGene<-rep(NA,length(pip10))
for(i in 1:length(pip10)){
	ch<-chromtab[chromtab[,1]==resGB[pip10][i],2]
	xx<-which(gtf[genes,1]==ch & resGB[pip10,2][i] > gtf[genes,4] & resGB[pip10,2][i] < gtf[genes,5])
	if(length(xx)>=1){
		GBinGene[i]<-gtf[genes[xx],9]
	}
}
## SNP 1 (8, 56786125) in g407 and g408 = RNA polymerage III and UBX domain



