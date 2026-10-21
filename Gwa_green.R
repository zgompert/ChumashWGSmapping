## summarize single SNP GWA
library(data.table)
library(scales)

## graphical parameters
cl<-1.4;ca<-1.1

## green genome ##

## get SNP data, need to drop
## full set from VCF
fullSNPs_green<-fread("../gendat_rev/FullSNPTable_green.txt",header=FALSE) 

## drop for geno
keepGreen<-scan("../gendat_rev/KeepSnps_green.txt")
subSNPs_green<-fullSNPs_green[keepGreen==1,]
dim(subSNPs_green)
#[1] 39442472        4
het<-as.matrix(fread("HetAndP_green.txt",header=TRUE))
## drop for het
keep<-which(het[,3] > -.2 & het[,3] < .2)
hsubSNPs_green<-subSNPs_green[keep,]
dim(hsubSNPs_green)
#[1] 39184098        4
write.table(hsubSNPs_green,file="hsubSNPs_green.txt",row.names=FALSE,col.names=FALSE,quote=FALSE) ## remove hap1_chr external
hsubSNPs_green<-fread("hsubSNPs_green.txt",header=FALSE)
snp_green<-cbind(as.numeric(unlist(hsubSNPs_green[,1])),as.numeric(unlist(hsubSNPs_green[,2])))

## read gwa results for green  RG
res_green_RG<-fread("output/o_green_RG.assoc.txt",header=TRUE)

maxSNPs<-tapply(INDEX=snp_green[,1],X=snp_green[,2],max)

## -log10 P from LRG
lp_green_RG<--log10(res_green_RG$p_lrt)
L<-length(lp_green_RG)
#sigRG<-which(lpRG > 2 & hetx[,3] > -.2 & hetx[,3] < .2 & hetx[,2] > .025 & hetx[,2] < .975) ## get rid of extreme hets and rare
sig_green_RG<-which(lp_green_RG > 2) ## get rid of extreme hets and rare
sig_green_RG_ch8<-which(lp_green_RG > 2 & snp_green[,1]==4)
all_green_RG_ch8<-which(snp_green[,1]==4)


## read gwa results for GB
res_green_GB<-fread("output/o_green_GB.assoc.txt",header=TRUE)

## -log10 P from LGB
lp_green_GB<--log10(res_green_GB$p_lrt)
L<-length(lp_green_GB)
sig_green_GB<-which(lp_green_GB > 2)
sig_green_GB_ch8<-which(lp_green_GB > 2 & snp_green[,1]==4)
all_green_GB_ch8<-which(snp_green[,1]==4)



relPos<-rep(NA,L)
for(i in 1:length(maxSNPs)){
	a<-which(snp_green[,1]==names(maxSNPs)[i])
	if(i==1){
		relPos[a]<-as.numeric(snp_green[a,2])/sum(maxSNPs)
	} else{
		relPos[a]<-(as.numeric(snp_green[a,2])+sum(maxSNPs[1:(i-1)]))/sum(maxSNPs)
	}
}
bnds<-relPos[c(which(snp_green[-L,1] != snp_green[-1,1])+.5,L+.5)]
mids<-relPos[tapply(X=1:L,INDEX=snp_green[,1],mean)]


pdf("fig_gwa_green_SSNP.pdf",width=8,height=18)
cm<-1.5;cl<-1.5;ca<-1.1
par(mfrow=c(6,1))
par(mar=c(4.5,5.5,2.5,1.5))
plot(relPos[sig_green_RG],lp_green_RG[sig_green_RG],ylim=c(2,30),type='n',axes=FALSE,xlab="Chromosome",ylab="-log10(P)",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,35,35),border=NA,col=alpha("gray20",.4))
}
points(relPos[sig_green_RG],lp_green_RG[sig_green_RG],pch=21,col=NA,bg=alpha("firebrick",.5))
axis(1,at=mids,unique(snp_green[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
title(main="(A) GWA RG for HF5, green genome",cex.main=cm)
box()

plot(relPos[sig_green_GB],lp_green_GB[sig_green_GB],ylim=c(2,34),type='n',axes=FALSE,xlab="Chromosome",ylab="-log10(P)",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,35,35),border=NA,col=alpha("gray20",.4))
}
points(relPos[sig_green_GB],lp_green_GB[sig_green_GB],pch=21,col=NA,bg=alpha("royalblue3",.5))
axis(1,at=mids,unique(snp_green[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
title(main="(B) GWA GB for HF5, green genome",cex.main=cm)
box()

plot(snp_green[sig_green_RG_ch8,2],lp_green_RG[sig_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(0,1e8))
title(main="(C) GWA Ch8 RG for HF5, green genome",cex.main=cm)
plot(snp_green[sig_green_GB_ch8,2],lp_green_GB[sig_green_GB_ch8],pch=21,col=NA,bg=alpha("royalblue3",.5),xlim=c(0,1e8))
title(main="(D) GWA Ch8 GB for HF5, green genome",cex.main=cm)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
title(main="(E) GWA Ch8 zoom RG for HF5, green genome",cex.main=cm)
plot(snp_green[all_green_GB_ch8,2],lp_green_GB[all_green_GB_ch8],pch=21,col=NA,bg=alpha("royalblue3",.5),xlim=c(5.8e7,6.7e7))
title(main="(F) GWA Ch8 zoom GB for HF5, green genome",cex.main=cm)


dev.off()

## plot from segmentation analysis
segs<-read.csv("ChumashColorRegions-green_genome.csv")

pdf("seg_regions_green.pdf",width=8,height=5)
par(mar=c(5,5,1,1))
ch8_snp_comb<-c(snp_green[all_green_GB_ch8,2],snp_green[all_green_RG_ch8,2])
ch8_lp_comb<-c(lp_green_GB[all_green_GB_ch8],lp_green_RG[all_green_RG_ch8])
xx<-order(ch8_snp_comb)
plot(ch8_snp_comb,ch8_lp_comb,pch=19,col=alpha(.15,colour="gray40"),xlim=c(5.8e7,6.6e7),xlab="Position (bp)",ylab="log10 P-value",cex.lab=1.4,cex.axis=1.1)
y<-seq(from=10,to=30,length.out=dim(segs)[1])
for(i in 1:dim(segs)[1]){
	for(a in seq(6,12,2)){
		b<-a+1
		lines(c(segs[i,a],segs[i,b]),c(y[i],y[i]),lwd=2)	
	}
}
dev.off()

## segmentation analyes
## x = position in bp, y = -log10 p
smoothSignal <- function(X, Y, winsize = 50000, winstep = 10000) {

    maxpos <- floor(max(X) / winsize) * winsize - 1

    winlb <- seq(from = 1, to = maxpos, by = winstep)
    winub <- winlb + winsize

    # First and last indices of X falling in each window [winlb, winub)
    first <- findInterval(winlb - 1, X) + 1
    last  <- findInterval(winub - 1, X)

    ss <- vapply(seq_along(first), function(k) {
        if (first[k] <= last[k]) {
            median(Y[first[k]:last[k]], na.rm = TRUE)
        } else {
            NA_real_
        }
    }, numeric(1))

    out<-cbind((winlb+winub)/2,ss)
    return(out)
}

## ss = smoothed signal, two columns, position and metric
## thresh = units above mean, minSz = minimum size of a peak
## minSep = minimum seperation between peaks
peakFinder<-function(ss=NA,thresh=0.5,minSz=200000,minSep=200000){
	N<-dim(ss)[1]
	## identify windows above the threshold
	pt<-as.numeric(ss[,2] >= median(ss[,2],na.rm=TRUE)+thresh)
	## get bounds
	lb<-which(pt[-N]==0 & pt[-1]==1)
	ub<-which(pt[-N]==1 & pt[-1]==0)
	sz<-ss[ub,1]-ss[lb,1]
	sep<-ss[lb[-1],1]-ss[ub[-length(ub)],1]
	## combine cases where sep is less than minSep
	while(min(sep) < minSep){
		msi<-which.min(sep)
		ub[msi]<-ub[msi+1]
		ub<-ub[-(msi+1)]
		lb<-lb[-(msi+1)]
		sz<-ss[ub,1]-ss[lb,1]
		sep<-ss[lb[-1],1]-ss[ub[-length(ub)],1]
	}
	## drop small ones
	keep<-which(sz >= minSz)
	o<-cbind(ss[lb[keep],1],ss[ub[keep],1],sz[keep])
	return(o)
}


## smooth over 200kbp
ch8_snp_comb<-c(snp_green[all_green_GB_ch8,2],snp_green[all_green_RG_ch8,2])
ch8_lp_comb<-c(lp_green_GB[all_green_GB_ch8],lp_green_RG[all_green_RG_ch8])
xx<-order(ch8_snp_comb)
o<-smoothSignal(X=ch8_snp_comb[xx],Y=ch8_lp_comb[xx],winsize=200000)
plot(o[,1],o[,2],type='l',xlim=c(5.8e7,6.6e7))


po<-peakFinder(ss=o,thresh=0.5,minSz=200000,minSep=200000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po## alla re lower, upper, size
#         [,1]     [,2]    [,3]
#[1,] 60040001 61640001 1600000
#[2,] 62480001 63860001 1380000
#[3,] 64490001 64890001  400000

po<-peakFinder(ss=o,thresh=0.75,minSz=200000,minSep=200000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60090001 61640001 1550000
#[2,] 62480001 63230001  750000
#[3,] 64490001 64800001  310000

po<-peakFinder(ss=o,thresh=1,minSz=200000,minSep=200000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60160001 60970001 810000
#[2,] 62480001 63150001 670000
#[3,] 64490001 64760001 270000

po<-peakFinder(ss=o,thresh=0.5,minSz=400000,minSep=300000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
#         [,1]     [,2]    [,3]
#[1,] 60040001 61830001 1790000
#[2,] 62450001 63880001 1430000
#[3,] 64490001 64890001  400000

po<-peakFinder(ss=o,thresh=.75,minSz=400000,minSep=300000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60090001 61640001 1550000
#[2,] 62480001 63460001  980000


po<-peakFinder(ss=o,thresh=1,minSz=400000,minSep=300000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60160001 60970001 810000
#[2,] 62480001 63460001 980000

po<-peakFinder(ss=o,thresh=0.5,minSz=300000,minSep=300000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60040001 61830001 1790000
#[2,] 62450001 63880001 1430000
#[3,] 64490001 64890001  400000

po<-peakFinder(ss=o,thresh=0.75,minSz=300000,minSep=300000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60090001 61640001 1550000
#[2,] 62480001 63460001  980000
#[3,] 64490001 64800001  310000

po<-peakFinder(ss=o,thresh=1,minSz=300000,minSep=300000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[0,] 60160001 60970001 810000
#[2,] 62480001 63460001 980000

po<-peakFinder(ss=o,thresh=0.5,minSz=100000,minSep=100000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
         [,1]     [,2]    [,3]
#[1,] 60040001 61640001 1600000
#[2,] 62450001 63280001  830000
#[3,] 63440001 63880001  440000
#[4,] 64490001 64890001  400000


po<-peakFinder(ss=o,thresh=0.75,minSz=100000,minSep=100000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60090001 61330001 1240000
#[2,] 61440001 61640001  200000
#[3,] 62480001 63230001  750000
#[4,] 64490001 64800001  310000


po<-peakFinder(ss=o,thresh=1,minSz=100000,minSep=100000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60160001 60320001 160000
#[2,] 60420001 60700001 280000
#[3,] 62480001 63150001 670000
#[4,] 64490001 64760001 270000


## smooth over 300kbp
ch8_snp_comb<-c(snp_green[all_green_GB_ch8,2],snp_green[all_green_RG_ch8,2])
ch8_lp_comb<-c(lp_green_GB[all_green_GB_ch8],lp_green_RG[all_green_RG_ch8])
xx<-order(ch8_snp_comb)
o<-smoothSignal(X=ch8_snp_comb[xx],Y=ch8_lp_comb[xx],winsize=300000)


po<-peakFinder(ss=o,thresh=0.5,minSz=200000,minSep=200000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po## all are lower, upper, size
#         [,1]     [,2]    [,3]
#[1,] 60050001 61690001 1640000
#[2,] 62430001 63880001 1450000
#[3,] 64440001 64860001  420000

po<-peakFinder(ss=o,thresh=0.75,minSz=200000,minSep=200000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60100001 61600001 1500000
#[2,] 62430001 63280001  850000
#[3,] 64440001 64760001  320000


po<-peakFinder(ss=o,thresh=1,minSz=200000,minSep=200000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60250001 61020001 770000
#[2,] 62430001 63160001 730000
#[3,] 64510001 64710001 200000

po<-peakFinder(ss=o,thresh=0.5,minSz=400000,minSep=300000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60050001 61690001 1640000
#[2,] 62430001 63880001 1450000
#[3,] 64440001 64860001  420000

po<-peakFinder(ss=o,thresh=.75,minSz=400000,minSep=300000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60100001 61600001 1500000
#[2,] 62430001 63280001  850000


po<-peakFinder(ss=o,thresh=1,minSz=400000,minSep=300000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60250001 61020001 770000
#[2,] 62430001 63160001 730000

po<-peakFinder(ss=o,thresh=0.5,minSz=300000,minSep=300000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60050001 61690001 1640000
#[2,] 62430001 63880001 1450000
#[3,] 64440001 64860001  420000

po<-peakFinder(ss=o,thresh=0.75,minSz=300000,minSep=300000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60100001 61600001 1500000
#[2,] 62430001 63280001  850000
#[3,] 64440001 64760001  320000

po<-peakFinder(ss=o,thresh=1,minSz=300000,minSep=300000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60250001 61020001 770000
#[2,] 62430001 63160001 730000

po<-peakFinder(ss=o,thresh=0.5,minSz=100000,minSep=100000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60050001 61690001 1640000
#[2,] 62430001 63330001  900000
#[3,] 63630001 63880001  250000
#[4,] 64440001 64860001  420000


po<-peakFinder(ss=o,thresh=0.75,minSz=100000,minSep=100000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60100001 61200001 1100000
#[2,] 61300001 61600001  300000
#[3,] 62430001 63280001  850000
#[4,] 64440001 64760001  320000

po<-peakFinder(ss=o,thresh=1,minSz=100000,minSep=100000)
plot(snp_green[all_green_RG_ch8,2],lp_green_RG[all_green_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60250001 60810001 560000
#[2,] 62430001 63160001 730000
#[3,] 64510001 64710001 200000

save(list=ls(),file="gwa_green.rdat")
