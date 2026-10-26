## summarize single SNP GWA for GR8
library(data.table)
library(scales)
setDTthreads(threads=10)


## graphical parameters
cl<-1.4;ca<-1.1

ph<-read.csv("2015_Tchumash_color.csv")
phGR8<-ph[which(ph$location=="GR8.06" & is.na(ph$latRG)==FALSE),]


## gwa ##

## chromosome ids
## last two are 159 h1 and 163 h1
chtab<-matrix(c(1,43,1,1, ## 1 vs 4 vs 9 (dropped) arbitrary
        2,1392,6,7,
        4,43,2,2,
        5,56,5,5,
        6,1469,7,6,
        7,1510,10,10,
        8,113,4,4,
        10,1213,8,8,
        11,48,9,9,
        12,1403,11,11,
        13,1308,3,3),nrow=11,ncol=4,byrow=TRUE)


## read gwa results for 
res_green_RG<-fread("output/o_GR8_RG.assoc.txt",header=TRUE)
res_green_GB<-fread("output/o_GR8_GB.assoc.txt",header=TRUE)

snp_green<-as.matrix(fread("snps_gr8_green.txt",header=FALSE))

## renumber chromosomes
newchrom<-rep(NA,dim(snp_green)[1])
for(k in 1:dim(chtab)[1]){
	sc<-chtab[k,4]
	ch<-chtab[k,1]
	a<-which(snp_green[,1]==sc)
	newchrom[a]<-ch
}
snp_green[,1]<-newchrom

sort_snp_green<-snp_green[order(snp_green[,1]),]
maxSNPs<-tapply(INDEX=sort_snp_green[,1],X=sort_snp_green[,2],max)

## -log10 P from LRG
lp_green_RG<--log10(res_green_RG$p_lrt[order(snp_green[,1])])
L<-length(lp_green_RG)
sig_green_RG<-which(lp_green_RG > 2) ## get rid of extreme hets and rare


## -log10 P from LGB
lp_green_GB<--log10(res_green_GB$p_lrt[order(snp_green[,1])])
L<-length(lp_green_GB)
sig_green_GB<-which(lp_green_GB > 2)

relPos<-rep(NA,L)
for(i in 1:length(maxSNPs)){
        a<-which(sort_snp_green[,1]==names(maxSNPs)[i])
        if(i==1){
                relPos[a]<-as.numeric(sort_snp_green[a,2])/sum(maxSNPs)
        } else{
                relPos[a]<-(as.numeric(sort_snp_green[a,2])+sum(maxSNPs[1:(i-1)]))/sum(maxSNPs)
        }
}
bnds<-relPos[c(which(sort_snp_green[-L,1] != sort_snp_green[-1,1])+.5,L+.5)]
mids<-relPos[tapply(X=1:L,INDEX=sort_snp_green[,1],mean)]

pdf("sfig_gwa_SSNP_GR8.pdf",width=8,height=6)
cm<-1.3;cl<-1.3;ca<-1.0
par(mfrow=c(2,1))

par(mar=c(4.5,4.5,2.5,1))
plot(relPos[sig_green_RG],lp_green_RG[sig_green_RG],ylim=c(2,34),type='n',axes=FALSE,xlab="Chromosome",ylab=expression(-log[10]~(P)),cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,35,35),border=NA,col=alpha("gray20",.4))
}
points(relPos[sig_green_RG],lp_green_RG[sig_green_RG],pch=21,col=NA,bg=alpha("firebrick",.5))
axis(1,at=mids,unique(sort_snp_green[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
title(main="(a) GR8 GWA RG for genome 24_0163H1 (green genome assembly)",cex.main=cm)
box()

#####
plot(relPos[sig_green_GB],lp_green_GB[sig_green_GB],ylim=c(2,39),type='n',axes=FALSE,xlab="Chromosome",ylab=expression(-log[10]~(P)),cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,40,40),border=NA,col=alpha("gray20",.4))
}
points(relPos[sig_green_GB],lp_green_GB[sig_green_GB],pch=21,col=NA,bg=alpha("royalblue3",.5))
axis(1,at=mids,unique(sort_snp_green[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
title(main="(b) GR8 GWA GB for genome 24_0163H1 (green genome assembly)",cex.main=cm)
box()

dev.off()

## ch8 figure
sigRG8<-which(lp_green_RG > 2 & sort_snp_green[,1]==8) 
sigGB8<-which(lp_green_GB > 2 & sort_snp_green[,1]==8) 

pdf("sfig_gwa_zoom_GR8.pdf",width=8,height=6)
cm<-1.3;cl<-1.3;ca<-1.0
par(mfrow=c(2,1))

par(mar=c(4.5,4.5,2.5,1))
plot(sort_snp_green[sigRG8,2],lp_green_RG[sigRG8],ylim=c(2,34),xlim=c(5.5e7,7.0e7),xlab="Position (bp)",ylab=expression(-log[10]~(P)),cex.lab=cl,pch=21,col=NA,bg=alpha("firebrick",.5))
title(main="(a) GR8 chromosome 8 GWA RG for genome 24_0163H1",cex.main=cm)

plot(sort_snp_green[sigGB8,2],lp_green_GB[sigGB8],ylim=c(2,39),xlim=c(5.5e7,7.0e7),xlab="Position (bp)",ylab=expression(-log[10]~(P)),cex.lab=cl,pch=21,col=NA,bg=alpha("royalblue3",.5))

title(main="(b) GR8 chromosome 8 GWA GB for genome 24_0163H1",cex.main=cm)
dev.off()
