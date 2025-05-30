## summariz## chrom and pos for RG, has extreme het removed already
snpsRG<-as.matrix(fread("SNPsRG_GR8.txt",header=FALSE))

e BSLMM GWA for GR*
library(data.table)
library(scales)

## graphical parameters
cl<-1.4;ca<-1.1

## read gwa results for RG
resRG<-fread("output/pips_o_GR8_poly_bslmm_ph1_ch1.param.txt",header=FALSE)

## chrom and pos for RG, has extreme het removed already
snpsRG<-as.matrix(fread("SNPsRG_GR8.txt",header=FALSE))

pips<-as.numeric(resRG$V2)
L<-length(pips)
sigRG<-which(pips > 0.00001) 

bnds<-c(which(snpsRG[-L,1] != snpsRG[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=snpsRG[,1],mean)


pdf("GWAGR8_RG_poly.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(sigRG,pips[sigRG],ylim=c(0,1),type='n',axes=FALSE,xlab="Chromosome",ylab="Posterior inclusion probability",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,28,28),border=NA,col=alpha("gray20",.4))
}
points(sigRG,pips[sigRG],pch=19)
axis(1,at=mids,unique(snpsRG[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
box()
dev.off()

sum(pips[snpsRG[,1]==8])
#[1] 3.273542
1-prod(1-pips[snpsRG[,1]==8])
#[1] 0.999203
sum(pips)
#[1] 11.11169
sum(pips[pips > 0.01])
#[1] 6.415333
sum(pips[snpsRG[,1]==8 & pips > 0.01])
#[1] 2.913792


## now for Chr8 only
sigRG8<-which(pips > 0.00001 & snpsRG[,1]==8)

pdf("GWAGR8_RG_poly_Ch8.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(snpsRG[sigRG8,2],pips[sigRG8],ylim=c(0,1),pch=19,xlab="Position (bp)",ylab="Posterior inclusion probability",cex.lab=cl)
dev.off()


## read gwa results for GB
resGB<-fread("output/pips_o_GR8_poly_bslmm_ph2_ch1.param.txt",header=FALSE)

## chrom and pos for GB
snpsGB<-snpsRG ## same

pips<-as.numeric(resGB$V2)
L<-length(pips)
sigGB<-which(pips > 0.00001) 

bnds<-c(which(snpsGB[-L,1] != snpsGB[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=snpsGB[,1],mean)

pdf("GWAGR8_GB_poly.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(sigGB,pips[sigGB],ylim=c(0,1),type='n',axes=FALSE,xlab="Chromosome",ylab="Posterior inclusion probability",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),border=NA,col=alpha("gray20",.4))
}
points(sigGB,pips[sigGB],pch=19)
axis(1,at=mids,unique(snpsGB[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
box()
dev.off()

## now for Chr8 only
sigGB8<-which(pips > 0.00001 & snpsGB[,1]==8)

pdf("GWAGR8_GB_poly_Ch8.pdf",width=8,height=4.5)
par(mar=c(5,5,1,1))
plot(snpsGB[sigGB8,2],pips[sigGB8],ylim=c(0,1),pch=19,xlab="Position (bp)",ylab="Posterior inclusion probability",cex.lab=cl)
dev.off()

sum(pips[snpsGB[,1]==8])
#[1] 3.417886
1-prod(1-pips[snpsGB[,1]==8])
#[1] 0.9788123
sum(pips)
#[1] 14.59631
sum(pips[pips > 0.01])
#[1] 2.930607
sum(pips[snpsGB[,1]==8 & pips > 0.01])
#[1] 2.498587
as.numeric(resGB$V2)
