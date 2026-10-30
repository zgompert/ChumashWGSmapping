library(data.table)
library(scales)

############## green genome ##########

## read in the data
G<-fread("hsub_mlpntest_tchum_green.geno",header=FALSE,sep=",")
gmat<-as.matrix(G[,-c(1:3)])

hsubSNPs_green<-fread("hsubSNPs_green.txt",header=FALSE)
snp_green<-cbind(as.numeric(unlist(hsubSNPs_green[,1])),as.numeric(unlist(hsubSNPs_green[,2])))
ch<-snp_green[,1]
pos<-snp_green[,2]

ch8<-which(ch==4)
gmat_ch8<-gmat[ch8,]
pos_ch8<-pos[ch8]

regions<-read.csv("ChumashColorRegions-green_genome.csv")

bnds<-seq(from=0,to=100000000,by=100000)
sz<-1000000
step<-100000
L<-length(bnds)-10
snpCnts<-rep(NA,L)
for(k in 1:L){
	a<-sum(pos_ch8 >= bnds[k] & pos_ch8 < bnds[k+10])
	snpCnts[k]<-a
}

plot(bnds[1:L]+500000,snpCnts,type='l')

anno<-read.table("~/../gompert-group4/data/timema/hic_genomes/Annotation/t_chumash_24_0163_g_h1/braker.gff3",header=FALSE)

anno_ch8<-anno[anno[,1]=="hap1_Chr04",]
cds_ch8<-anno_ch8[anno_ch8[,3]=="CDS",]
cds_ch8_lb<-cds_ch8[,4]
cds_ch8_ub<-cds_ch8[,5]

sz<-1000000
step<-100000
L<-length(bnds)-10
geneCnts<-rep(NA,L)
for(k in 1:L){
	gs<-rep(0,sz)
	a<-which(cds_ch8_lb >= bnds[k] & cds_ch8_ub <= bnds[k+10])
	lb<-cds_ch8_lb[a]%%sz;ub<-cds_ch8_ub[a]%%sz
	if(length(a) > 0){
		for(z in 1:length(a)){
			gs[lb[z]:ub[z]]<-1
		}
	}
	a<-which(cds_ch8_lb >= bnds[k] & cds_ch8_lb <= bnds[k+10] & cds_ch8_ub > bnds[k+10])
	lb<-cds_ch8_lb[a]%%sz;ub<-rep(bnds[k+1],length(a))%%sz
	if(length(a) > 0){
		for(z in 1:length(a)){
			gs[lb[z]:ub[z]]<-1
		}
	}
	a<-which(cds_ch8_lb <= bnds[k] & cds_ch8_ub >= bnds[k] & cds_ch8_ub <= bnds[k+10])
	if(length(a) > 0){
		lb<-rep(bnds[k+1],length(a))%%sz;ub<-cds_ch8_ub[a]%%sz
		for(z in 1:length(a)){
			gs[lb[z]:ub[z]]<-1
		}
	}

	geneCnts[k]<-mean(gs)
}

library(scales)
pdf("sfig_green_density.pdf",width=8,height=9)
par(mfrow=c(2,1))
par(mar=c(5,5,3,1))
plot(bnds[1:L]+500000,snpCnts,type='n',xlab="Position (bp) on chromosome 8",ylab="SNPs per Mbp",cex.lab=1.3,cex.axis=1)
polygon(c(regions[1,c(6,11)],rev(regions[1,c(6,11)])),c(-1e4,-1e4,1e5,1e5),border=NA,col=alpha("gray",.5))
lines(bnds[1:L]+500000,snpCnts)

title(main="(a) SNP density",cex.main=1.3)
plot(bnds[1:L]+500000,geneCnts,type='n',xlab="Position (bp) on chromosome 8",ylab="% CDS per Mbp",cex.lab=1.3,cex.axis=1,ylim=c(0,.2))
polygon(c(regions[1,c(6,11)],rev(regions[1,c(6,11)])),c(-1,-1,1,1),border=NA,col=alpha("gray",.5))
lines(bnds[1:L]+500000,geneCnts)
title(main="(b) Gene density",cex.main=1.3)
dev.off()
