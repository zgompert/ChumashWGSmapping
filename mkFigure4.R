## make figure 4, initial LD figure and PCAs 
library(data.table)
library(scales)

## graphical parameters
cl<-1.4;ca<-1.1;cm<-1.5


##################################################################
G<-fread("hsub_mlpntest_tchum_green.geno",header=FALSE,sep=",")
#ph<-read.table("ph_tchum_green.txt",header=FALSE)
ph<-read.table("sort_green_phfull_tchum2019.txt",header=TRUE,comment.char="%",sep=",")

dim(G)
#39184098      458
gmat<-as.matrix(G[,-c(1:3)])
hsubSNPs_green<-fread("hsubSNPs_green.txt",header=FALSE)
snp_green<-cbind(as.numeric(unlist(hsubSNPs_green[,1])),as.numeric(unlist(hsubSNPs_green[,2])))

ch8<-which(snp_green==4)## scaffold 4 is chromosome 8
pos<-snp_green[ch8,2]
g_ch8<-gmat[ch8,]
p_ch8<-apply(g_ch8,1,mean)/2

dim(g_ch8)
#[1] 3278110     455

rm("G") ## for space

## read gwa results for green  RG
res_green_RG<-fread("output/o_green_RG.assoc.txt",header=TRUE)

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

segs<-read.table("ChumashColorRegions-green_genome.csv",header=TRUE,sep=",")
ch8_snp_comb<-c(snp_green[all_green_GB_ch8,2],snp_green[all_green_RG_ch8,2])
ch8_lp_comb<-c(lp_green_GB[all_green_GB_ch8],lp_green_RG[all_green_RG_ch8])


### figure 4 manuscript #####
pdf("fig_LD_PCA.pdf",width=8.4,height=10)
layout(matrix(c(1,1,2,2,3,3,3,3,4,5,5,6),nrow=3,ncol=4,byrow=TRUE),widths=c(2,1,1,2),heights=c(3,3,2))
par(mar=c(4.5,5.5,2.5,1.5))
##########
set<-which(p_ch8 > .05 & p_ch8 < .95 &  pos > 58500000 & pos < 66500000)
pset<-pos[set]
wins<-seq(from=58.5e6,to=66.6e6,by=10000)
cnts<-rep(NA,length(wins)-1)
for(k in 1:(length(wins)-1)){
    cnts[k]<-sum(pset >= wins[k] & pset < wins[k+1])
    }
sset<-rep(NA,sum(cnts > 0))  
a<-1  
for(k in 1:(length(wins)-1)){
    if(cnts[k] > 0){
        if(cnts[k] > 1){sset[a]<-as.numeric(sample(x=set[pset >= wins[k] & pset < wins[k+1]],1))} else{
            sset[a]<-as.numeric(set[pset >= wins[k] & pset < wins[k+1]])}
            a<-a+1
        }

        }
        

r2<-cor(t(g_ch8[sset,]))^2
r2[upper.tri(r2)]<-NA
diag(r2)<-NA
image(r2,breaks=c(-.1,.01,.025,.05,.1,1.1),col=heat.colors(n=5,rev=TRUE),axes=FALSE)
lb<-min(which(pos[sset] > 60000000))/length(sset)
ub<-max(which(pos[sset] < 65000000))/length(sset)


lines(c(lb,ub),c(lb,lb))
lines(c(ub,ub),c(lb,ub))
legend(0,.95,c("<0.01","0.01-0.025","0.025-0.05","0.05-0.10",">0.10"),fill=heat.colors(n=5,rev=TRUE),bty='n',title="LD",cex=1.3)
title(main="(a) LD chromosome 8",cex.main=cm)

##########
regs<-read.table("ChumashColorRegions-green_genome.csv",header=TRUE,sep=",")
r1lb<-regs$window1_lb[1]
r2lb<-regs$window2_lb[1]
r3lb<-regs$window3_lb[1]
r1ub<-regs$window1_ub[1]
r2ub<-regs$window2_ub[1]
r3ub<-regs$window3_ub[1]
r1<-which(pos >= r1lb & pos <= r1ub)
r2<-which(pos >= r2lb & pos <= r2ub)
r3<-which(pos >= r3lb & pos <= r3ub)
o1<-prcomp(t(g_ch8[r1,]),center=TRUE,scale=FALSE)
o2<-prcomp(t(g_ch8[r2,]),center=TRUE,scale=FALSE)
o3<-prcomp(t(g_ch8[r3,]),center=TRUE,scale=FALSE)

pcs<-cbind(o1$x[,1:2],o2$x[,1:2],o3$x[,1:2])
rr<-abs(cor(pcs))
rr[upper.tri(rr)]<-NA
diag(rr)<-NA
image(rr,breaks=c(-.1,.2,.4,.6,.8,1.1),col=heat.colors(n=5,rev=TRUE),axes=FALSE)
aaa<-seq(0,1,.2)
text(aaa,aaa,c("PC1r1","PC2r1","PC1r2","PC2r2","PC1r3","PC2r3"),cex=ca)
legend(0,1.05,c("<0.2","0.2-0.4","0.4-0.6","0.6-0.8",">0.8"),fill=heat.colors(n=5,rev=TRUE),bty='n',title="|r|",cex=1.3)
title(main="(b) PC correlations",cex.main=cm)

xx<-order(ch8_snp_comb)
zz<-which(ch8_snp_comb >=5.7e7 & ch8_snp_comb <= 6.8e7)
plot(ch8_snp_comb[zz],ch8_lp_comb[zz],type='n',xlim=c(5.7e7,6.8e7),xlab="Position (bp)",ylab="log10 P-value",cex.lab=1.4,cex.axis=1.1)
y<-seq(from=10,to=30,length.out=dim(segs)[1])
i<-1
for(a in seq(6,12,2)){
	b<-a+1
        polygon(c(segs[i,a],segs[i,b],segs[i,b],segs[i,a]),c(-10,-10,50,50),col=alpha(.2,colour="firebrick"),border=NA)        
}
points(ch8_snp_comb[zz],ch8_lp_comb[zz],pch=19,col=alpha(.15,colour="gray40"))
mns<-c((r1lb+r1ub)/2,(r2lb+r2ub)/2,(r3lb+r3ub)/2)
text(mns,c(33,33,33),c("R1","R2","R3"),cex=1.3)
title(main="(c) Colour-associated regions of chromosome 8",cex.main=cm)



plot(o1$x[,1],o1$x[,2],pch=19,col=ph$hex,cex.lab=cl,cex.axis=ca,xlab="PC1",ylab="PC2")
title(main="(d) PCA region 1",cex.main=cm)

plot(o2$x[,1],o2$x[,2],pch=19,col=ph$hex,cex.lab=cl,cex.axis=ca,xlab="PC1",ylab="PC2")
title(main="(e) PCA region 2",cex.main=cm)

plot(o3$x[,1],o3$x[,2],pch=19,col=ph$hex,cex.lab=cl,cex.axis=ca,xlab="PC1",ylab="PC2")
title(main="(f) PCA region 3",cex.main=cm)

dev.off()
