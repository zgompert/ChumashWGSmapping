## summarize single SNP GWA
library(data.table)
library(scales)

## graphical parameters
cl<-1.4;ca<-1.1

## red genome ##

## get SNP data, need to drop
## full set from VCF
fullSNPs_red<-fread("../gendat_rev/FullSNPTable_red.txt",header=FALSE) 

## drop for geno
keepGreen<-scan("../gendat_rev/KeepSnps_red.txt")
subSNPs_red<-fullSNPs_red[keepGreen==1,]
dim(subSNPs_red)
#[1] 39216576         4
het<-as.matrix(fread("HetAndP_red.txt",header=TRUE))
## drop for het
keep<-which(het[,3] > -.2 & het[,3] < .2)
hsubSNPs_red<-subSNPs_red[keep,]
dim(hsubSNPs_red)
#[1] 38967857        4
write.table(hsubSNPs_red,file="hsubSNPs_red.txt",row.names=FALSE,col.names=FALSE,quote=FALSE) ## remove hap1_chr external
hsubSNPs_red<-fread("hsubSNPs_red.txt",header=FALSE)
snp_red<-cbind(as.numeric(unlist(hsubSNPs_red[,1])),as.numeric(unlist(hsubSNPs_red[,2])))

## read gwa results for red  RG
res_red_RG<-fread("output/o_red_RG.assoc.txt",header=TRUE)

maxSNPs<-tapply(INDEX=snp_red[,1],X=snp_red[,2],max)

## -log10 P from LRG
lp_red_RG<--log10(res_red_RG$p_lrt)
L<-length(lp_red_RG)
#sigRG<-which(lpRG > 2 & hetx[,3] > -.2 & hetx[,3] < .2 & hetx[,2] > .025 & hetx[,2] < .975) ## get rid of extreme hets and rare
sig_red_RG<-which(lp_red_RG > 2) ## get rid of extreme hets and rare
sig_red_RG_ch8<-which(lp_red_RG > 2 & snp_red[,1]==4)
all_red_RG_ch8<-which(snp_red[,1]==4)


## read gwa results for GB
res_red_GB<-fread("output/o_red_GB.assoc.txt",header=TRUE)

## -log10 P from LGB
lp_red_GB<--log10(res_red_GB$p_lrt)
L<-length(lp_red_GB)
sig_red_GB<-which(lp_red_GB > 2)
sig_red_GB_ch8<-which(lp_red_GB > 2 & snp_red[,1]==4)
all_red_GB_ch8<-which(snp_red[,1]==4)



relPos<-rep(NA,L)
for(i in 1:length(maxSNPs)){
	a<-which(snp_red[,1]==names(maxSNPs)[i])
	if(i==1){
		relPos[a]<-as.numeric(snp_red[a,2])/sum(maxSNPs)
	} else{
		relPos[a]<-(as.numeric(snp_red[a,2])+sum(maxSNPs[1:(i-1)]))/sum(maxSNPs)
	}
}
bnds<-relPos[c(which(snp_red[-L,1] != snp_red[-1,1])+.5,L+.5)]
mids<-relPos[tapply(X=1:L,INDEX=snp_red[,1],mean)]


pdf("fig_gwa_red_SSNP.pdf",width=8,height=18)
cm<-1.5;cl<-1.5;ca<-1.1
par(mfrow=c(6,1))
par(mar=c(4.5,5.5,2.5,1.5))
plot(relPos[sig_red_RG],lp_red_RG[sig_red_RG],ylim=c(2,30),type='n',axes=FALSE,xlab="Chromosome",ylab="-log10(P)",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,35,35),border=NA,col=alpha("gray20",.4))
}
points(relPos[sig_red_RG],lp_red_RG[sig_red_RG],pch=21,col=NA,bg=alpha("firebrick",.5))
axis(1,at=mids,unique(snp_red[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
title(main="(A) GWA RG for HF5, red genome",cex.main=cm)
box()

plot(relPos[sig_red_GB],lp_red_GB[sig_red_GB],ylim=c(2,34),type='n',axes=FALSE,xlab="Chromosome",ylab="-log10(P)",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,35,35),border=NA,col=alpha("gray20",.4))
}
points(relPos[sig_red_GB],lp_red_GB[sig_red_GB],pch=21,col=NA,bg=alpha("royalblue3",.5))
axis(1,at=mids,unique(snp_red[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
title(main="(B) GWA GB for HF5, red genome",cex.main=cm)
box()

plot(snp_red[sig_red_RG_ch8,2],lp_red_RG[sig_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(0,1e8))
title(main="(C) GWA Ch8 RG for HF5, red genome",cex.main=cm)
plot(snp_red[sig_red_GB_ch8,2],lp_red_GB[sig_red_GB_ch8],pch=21,col=NA,bg=alpha("royalblue3",.5),xlim=c(0,1e8))
title(main="(D) GWA Ch8 GB for HF5, red genome",cex.main=cm)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
title(main="(E) GWA Ch8 zoom RG for HF5, red genome",cex.main=cm)
plot(snp_red[all_red_GB_ch8,2],lp_red_GB[all_red_GB_ch8],pch=21,col=NA,bg=alpha("royalblue3",.5),xlim=c(5.8e7,6.7e7))
title(main="(F) GWA Ch8 zoom GB for HF5, red genome",cex.main=cm)


dev.off()

## plot from segmentation analysis
segs<-read.csv("ChumashColorRegions-red_genome.csv")

pdf("seg_regions_red.pdf",width=8,height=5)
par(mar=c(5,5,1,1))
ch8_snp_comb<-c(snp_red[all_red_GB_ch8,2],snp_red[all_red_RG_ch8,2])
ch8_lp_comb<-c(lp_red_GB[all_red_GB_ch8],lp_red_RG[all_red_RG_ch8])
xx<-order(ch8_snp_comb)
plot(ch8_snp_comb,ch8_lp_comb,pch=19,col=alpha(.15,colour="gray40"),xlim=c(5.8e7,6.6e7),xlab="Position (bp)",ylab="log10 P-value",cex.lab=1.4,cex.axis=1.1)
y<-seq(from=10,to=30,length.out=dim(segs)[1])
for(i in 1:dim(segs)[1]){
        for(a in seq(6,10,2)){
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
ch8_snp_comb<-c(snp_red[all_red_GB_ch8,2],snp_red[all_red_RG_ch8,2])
ch8_lp_comb<-c(lp_red_GB[all_red_GB_ch8],lp_red_RG[all_red_RG_ch8])
xx<-order(ch8_snp_comb)
o<-smoothSignal(X=ch8_snp_comb[xx],Y=ch8_lp_comb[xx],winsize=200000)


po<-peakFinder(ss=o,thresh=0.5,minSz=200000,minSep=200000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po## alla re lower, upper, size
#         [,1]     [,2]    [,3]
#[1,] 60690001 62040001 1350000
#[2,] 62820001 63890001 1070000

po<-peakFinder(ss=o,thresh=0.75,minSz=200000,minSep=200000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60720001 62010001 1290000
#[2,] 62830001 63790001  960000

po<-peakFinder(ss=o,thresh=1,minSz=200000,minSep=200000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60780001 61740001 960000
#[2,] 62850001 63160001 310000
#[3,] 63380001 63760001 380000

po<-peakFinder(ss=o,thresh=0.5,minSz=400000,minSep=300000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
#         [,1]     [,2]    [,3]
#[1,] 60690001 62390001 1700000
#[2,] 62820001 63890001 1070000

po<-peakFinder(ss=o,thresh=.75,minSz=400000,minSep=300000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60720001 62260001 1540000
#[2,] 62830001 63790001  960000

po<-peakFinder(ss=o,thresh=1,minSz=400000,minSep=300000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60780001 61740001 960000
#[2,] 62850001 63760001 910000

po<-peakFinder(ss=o,thresh=0.5,minSz=300000,minSep=300000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60690001 62390001 1700000
#[2,] 62820001 63890001 1070000

po<-peakFinder(ss=o,thresh=0.75,minSz=300000,minSep=300000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60720001 62260001 1540000
#[2,] 62830001 63790001  960000

po<-peakFinder(ss=o,thresh=1,minSz=300000,minSep=300000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60780001 61740001 960000
#[2,] 62850001 63760001 910000

po<-peakFinder(ss=o,thresh=0.5,minSz=100000,minSep=100000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60690001 62040001 1350000
#[2,] 62240001 62390001  150000
#[3,] 62820001 63890001 1070000

po<-peakFinder(ss=o,thresh=0.75,minSz=100000,minSep=100000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60720001 62010001 1290000
#[2,] 62830001 63790001  960000

po<-peakFinder(ss=o,thresh=1,minSz=100000,minSep=100000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60780001 61740001 960000
#[2,] 62850001 63020001 170000
#[3,] 63380001 63760001 380000

## smooth over 300kbp
ch8_snp_comb<-c(snp_red[all_red_GB_ch8,2],snp_red[all_red_RG_ch8,2])
ch8_lp_comb<-c(lp_red_GB[all_red_GB_ch8],lp_red_RG[all_red_RG_ch8])
xx<-order(ch8_snp_comb)
o<-smoothSignal(X=ch8_snp_comb[xx],Y=ch8_lp_comb[xx],winsize=300000)

po<-peakFinder(ss=o,thresh=0.5,minSz=200000,minSep=200000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po## all are lower, upper, size
#         [,1]     [,2]    [,3]
#[1,] 60680001 62090001 1410000
#[2,] 62780001 63860001 1080000

po<-peakFinder(ss=o,thresh=0.75,minSz=200000,minSep=200000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60740001 61960001 1220000
#[2,] 62790001 63760001  970000

po<-peakFinder(ss=o,thresh=1,minSz=200000,minSep=200000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60790001 61780001 990000
#[2,] 62820001 63100001 280000
#[3,] 63330001 63730001 400000

po<-peakFinder(ss=o,thresh=0.5,minSz=400000,minSep=300000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60680001 62340001 1660000
#[2,] 62780001 63860001 1080000

po<-peakFinder(ss=o,thresh=.75,minSz=400000,minSep=300000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60740001 61960001 1220000
#[2,] 62790001 63760001  970000

po<-peakFinder(ss=o,thresh=1,minSz=400000,minSep=300000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60790001 61780001 990000
#[2,] 62820001 63730001 910000

po<-peakFinder(ss=o,thresh=0.5,minSz=300000,minSep=300000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60680001 62340001 1660000
#[2,] 62780001 63860001 1080000

po<-peakFinder(ss=o,thresh=0.75,minSz=300000,minSep=300000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60740001 61960001 1220000
#[2,] 62790001 63760001  970000

po<-peakFinder(ss=o,thresh=1,minSz=300000,minSep=300000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60790001 61780001 990000
#[2,] 62820001 63730001 910000

po<-peakFinder(ss=o,thresh=0.5,minSz=100000,minSep=100000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60680001 62090001 1410000
#[2,] 62780001 63860001 1080000

po<-peakFinder(ss=o,thresh=0.75,minSz=100000,minSep=100000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]    [,3]
#[1,] 60740001 61960001 1220000
#[2,] 62790001 63760001  970000

po<-peakFinder(ss=o,thresh=1,minSz=100000,minSep=100000)
plot(snp_red[all_red_RG_ch8,2],lp_red_RG[all_red_RG_ch8],pch=21,col=NA,bg=alpha("firebrick",.5),xlim=c(5.8e7,6.7e7))
abline(v=po[,1],lty=1)
abline(v=po[,2],lty=2)
po
#         [,1]     [,2]   [,3]
#[1,] 60790001 61780001 990000
#[2,] 62820001 63100001 280000
#[3,] 63330001 63730001 400000
save(list=ls(),file="gwa_red.rdat")
