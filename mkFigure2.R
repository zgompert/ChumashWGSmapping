## summarize single SNP GWA
library(data.table)
library(scales)

## graphical parameters
cl<-1.4;ca<-1.1

ph2019<-read.csv("/uufs/chpc.utah.edu/common/home/u6000989/../gompert-group4/projects/timema_genarch/color_exp_2019/2019_Tchumash_transplant_table.csv")

ph<-read.csv("2015_Tchumash_color.csv")
phGR8<-ph[which(ph$location=="GR8.06" & is.na(ph$latRG)==FALSE),]
colGR8<-rep(NA,dim(phGR8)[1])
## approximate hex
for(i in 1:length(colGR8)){
	colGR8[i]<-ph2019$hex[which.min(abs(phGR8$latRG[i]-ph2019$RG)+abs(phGR8$latGB[i]-ph2019$GB))]
}


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
#write.table(hsubSNPs_green,file="hsubSNPs_green.txt",row.names=FALSE,col.names=FALSE,quote=FALSE) ## remove hap1_chr external
hsubSNPs_green<-fread("hsubSNPs_green.txt",header=FALSE)
snp_green<-cbind(as.numeric(unlist(hsubSNPs_green[,1])),as.numeric(unlist(hsubSNPs_green[,2])))

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

## read gwa results for green  RG
res_green_RG<-fread("output/o_green_RG.assoc.txt",header=TRUE)


## -log10 P from LRG
lp_green_RG<--log10(res_green_RG$p_lrt[order(snp_green[,1])])
L<-length(lp_green_RG)
#sigRG<-which(lpRG > 2 & hetx[,3] > -.2 & hetx[,3] < .2 & hetx[,2] > .025 & hetx[,2] < .975) ## get rid of extreme hets and rare
sig_green_RG<-which(lp_green_RG > 2) ## get rid of extreme hets and rare

## read gwa results for GB
res_green_GB<-fread("output/o_green_GB.assoc.txt",header=TRUE)

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

## red genome ##

## get SNP data, need to drop
## full set from VCF
fullSNPs_red<-fread("../gendat_rev/FullSNPTable_red.txt",header=FALSE)

## drop for geno
keepRed<-scan("../gendat_rev/KeepSnps_red.txt")
subSNPs_red<-fullSNPs_red[keepRed==1,]
dim(subSNPs_red)
#[1] 39216576        4
het<-as.matrix(fread("HetAndP_red.txt",header=TRUE))
## drop for het
keep<-which(het[,3] > -.2 & het[,3] < .2)
hsubSNPs_red<-subSNPs_red[keep,]
dim(hsubSNPs_red)
#[1] 39184098        4
#write.table(hsubSNPs_red,file="hsubSNPs_red.txt",row.names=FALSE,col.names=FALSE,quote=FALSE) ## remove hap1_chr external
hsubSNPs_red<-fread("hsubSNPs_red.txt",header=FALSE)
snp_red<-cbind(as.numeric(unlist(hsubSNPs_red[,1])),as.numeric(unlist(hsubSNPs_red[,2])))

## renumber chromosomes
newchrom<-rep(NA,dim(snp_red)[1])
for(k in 1:dim(chtab)[1]){
	sc<-chtab[k,3]
	ch<-chtab[k,1]
	a<-which(snp_red[,1]==sc)
	newchrom[a]<-ch
}
snp_red[,1]<-newchrom
sort_snp_red<-snp_red[order(snp_red[,1]),]
maxSNPs_red<-tapply(INDEX=sort_snp_red[,1],X=sort_snp_red[,2],max)

## read gwa results for red  RG
res_red_RG<-fread("output/o_red_RG.assoc.txt",header=TRUE)


## -log10 P from LRG
lp_red_RG<--log10(res_red_RG$p_lrt[order(snp_red[,1])])
L<-length(lp_red_RG)
sig_red_RG<-which(lp_red_RG > 2) ## get rid of extreme hets and rare

## read gwa results for GB
res_red_GB<-fread("output/o_red_GB.assoc.txt",header=TRUE)

## -log10 P from LGB
lp_red_GB<--log10(res_red_GB$p_lrt[order(snp_red[,1])])
L<-length(lp_red_GB)
sig_red_GB<-which(lp_red_GB > 2)


relPos_red<-rep(NA,L)
for(i in 1:length(maxSNPs_red)){
        a<-which(sort_snp_red[,1]==names(maxSNPs_red)[i])
        if(i==1){
                relPos_red[a]<-as.numeric(sort_snp_red[a,2])/sum(maxSNPs_red)
        } else{
                relPos_red[a]<-(as.numeric(sort_snp_red[a,2])+sum(maxSNPs_red[1:(i-1)]))/sum(maxSNPs_red)
        }
}

bnds_red<-relPos_red[c(which(sort_snp_red[-L,1] != sort_snp_red[-1,1])+.5,L+.5)]
mids_red<-relPos_red[tapply(X=1:L,INDEX=sort_snp_red[,1],mean)]



### figure 2 for manuscript #####
pdf("fig_gwa_SSNP.pdf",width=8,height=12)
cm<-1.5;cl<-1.5;ca<-1.1
layout(matrix(c(1:8,rep(9:12,each=4)),nrow=6,ncol=4,byrow=TRUE),width=c(3,1,3,1),
       heights=c(1,3,2,2,2,2))

par(mar=c(0,4.5,2.5,0))
hist(ph2019$RG,axes=FALSE,xlab="",ylab="",main="",col=alpha("firebrick",.5),border=NA)
title(main="(a) Colour for HF5",cex.main=cm)
plot(c(1,1),type='n',xlab="",ylab="",axes=FALSE)
hist(phGR8$latRG,axes=FALSE,xlab="",ylab="",main="",col=alpha("firebrick",.5),border=NA)
title(main="(b) Colour for GR8",cex.main=cm)
plot(c(1,1),type='n',xlab="",ylab="",axes=FALSE)

par(mar=c(4.5,4.5,0,0))
plot(ph2019$RG,ph2019$GB,pch=19,col=ph2019$hex,xlab="RG",ylab="GB",cex.lab=cl,cex.axis=ca)
par(mar=c(4.5,0,0,1.5))
hh<-hist(ph2019$GB,plot=FALSE)
barplot(hh$counts,horiz=TRUE,names.arg="",xlab="",ylab="",space=0,axes=FALSE,col=alpha("royalblue3",.5),border=NA)
par(mar=c(4.5,4.5,0,0))
plot(phGR8$latRG,phGR8$latGB,pch=19,col=colGR8,xlab="RG",ylab="GB",cex.lab=cl,cex.axis=ca)
par(mar=c(4.5,0,0,1.5))
hh<-hist(phGR8$latGB,plot=FALSE)
barplot(hh$counts,horiz=TRUE,names.arg="",xlab="",ylab="",space=0,axes=FALSE,col=alpha("royalblue3",.5),border=NA)

par(mar=c(4.5,5.5,2.5,1.5))
plot(relPos[sig_green_RG],lp_green_RG[sig_green_RG],ylim=c(2,30),type='n',axes=FALSE,xlab="Chromosome",ylab="-log10(P)",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,31,31),border=NA,col=alpha("gray20",.4))
}
points(relPos[sig_green_RG],lp_green_RG[sig_green_RG],pch=21,col=NA,bg=alpha("firebrick",.5))
axis(1,at=mids,unique(sort_snp_green[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
title(main="(c) GWA RG for genome 24_0163H1 (green genome assembly)",cex.main=cm)
box()
####
plot(relPos_red[sig_red_RG],lp_red_RG[sig_red_RG],ylim=c(2,30),type='n',axes=FALSE,xlab="Chromosome",ylab="-log10(P)",cex.lab=cl)
N<-length(bnds_red)
for(i in seq(1,N-1,2)){
        polygon(c(bnds_red[i],bnds_red[i+1],bnds_red[i+1],bnds_red[i]),c(-1,-1,31,31),border=NA,col=alpha("gray20",.4))
}
points(relPos_red[sig_red_RG],lp_red_RG[sig_red_RG],pch=21,col=NA,bg=alpha("firebrick",.5))
axis(1,at=mids_red,unique(sort_snp_red[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
title(main="(d) GWA RG for genome 24_0159H1 (red genome assembly)",cex.main=cm)
box()
####

#####
plot(relPos[sig_green_GB],lp_green_GB[sig_green_GB],ylim=c(2,34),type='n',axes=FALSE,xlab="Chromosome",ylab="-log10(P)",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-1,-1,35,35),border=NA,col=alpha("gray20",.4))
}
points(relPos[sig_green_GB],lp_green_GB[sig_green_GB],pch=21,col=NA,bg=alpha("royalblue3",.5))
axis(1,at=mids,unique(sort_snp_green[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
title(main="(e) GWA GB for genome 24_0163H1 (green genome assembly)",cex.main=cm)
box()

#####
plot(relPos_red[sig_red_GB],lp_red_GB[sig_red_GB],ylim=c(2,35.5),type='n',axes=FALSE,xlab="Chromosome",ylab="-log10(P)",cex.lab=cl)
N<-length(bnds)
for(i in seq(1,N-1,2)){
        polygon(c(bnds_red[i],bnds_red[i+1],bnds_red[i+1],bnds_red[i]),c(-1,-1,36,36),border=NA,col=alpha("gray20",.4))
}
points(relPos[sig_red_GB],lp_red_GB[sig_red_GB],pch=21,col=NA,bg=alpha("royalblue3",.5))
axis(1,at=mids_red,unique(sort_snp_red[,1]),cex.axis=ca)
axis(2,cex.axis=ca)
title(main="(f) GWA GB for genome 24_0159H1 (red genome assembly)",cex.main=cm)
box()

dev.off()


