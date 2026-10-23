## make figure 5, PC change 
library(data.table)
library(scales)
library(RColorBrewer)

## graphical parameters
cl<-1.4;ca<-1.1;cm<-1.5

load("/uufs/chpc.utah.edu/common/home/u6000989/../gompert-group4/projects/timema_genarch/color_exp_2019/bayes_color_2019.rdat")

load("LDChange_pcs.rdat")
# selection coefficients from selection gradient, using alphas, hierarchical level effects
bMM<-est_hcl_mm[1:6,3]
bAC<-est_hcl_ac[1:6,3]
ph<-read.table("ph_tchum_green.txt",header=FALSE)
pmns<-apply(ph,2,mean,na.rm=TRUE)
psds<-apply(ph,2,sd,na.rm=TRUE)

gridsurf<-function(b=NA,xx=seq(-1.7,3.3,.02),yy=seq(-2.1,2.8,.02)){
        gm<-matrix(NA,nrow=length(yy),ncol=length(xx))
        lgm<-gm
        for(i in 1:length(yy)){
                lgm[i,]<-b[6] + b[1]*xx + b[2]*yy[i] + b[3]*xx^2 + b[4]*yy[i]^2 + b[5] * xx * yy[i]
        }
        gm<-1/(1+exp(-1 * lgm))
        gm<-gm/max(as.vector(gm))
        gm
        
}

gridsurf2<-function(b=NA,xx=seq(-1.7,3.3,.02),yy=seq(-2.1,2.8,.02)){
        gm<-matrix(NA,nrow=length(yy),ncol=length(xx))
        lgm<-gm
        for(i in 1:length(yy)){
                lgm[i,]<-b[6] + b[1]*xx + b[2]*yy[i] + b[3]*xx^2 + b[4]*yy[i]^2 + b[5] * xx * yy[i]
        }
        gm<-1/(1+exp(-1 * lgm))
        gm<-gm/max(as.vector(gm))
        gm
        
}

calcfit<-function(b=NA,xx=NA,yy=NA){
        lgm<-b[6] + b[1]*xx + b[2]*yy + b[3]*xx^2 + b[4]*yy^2 + b[5] * xx * yy
        gm<-1/(1+exp(-1 * lgm))
        gm
}

##############################################################
G<-fread("hsub_mlpntest_tchum_green.geno",header=FALSE,sep=",")
#ph<-read.table("ph_tchum_green.txt",header=FALSE)
ph<-read.table("sort_green_phfull_tchum2019.txt",header=TRUE,comment.char="%",sep=",")
phs<-read.table("../gemma/sort_phfull_tchum2019.txt",header=TRUE,comment.char="?")


dim(G)
#39184098      458
gmat<-as.matrix(G[,-c(1:3)])
hsubSNPs_green<-fread("hsubSNPs_green.txt",header=FALSE)
snp_green<-cbind(as.numeric(unlist(hsubSNPs_green[,1])),as.numeric(unlist(hsubSNPs_green[,2])))

ch8<-which(snp_green==4)## scaffold 4 is chromosome 8
pos<-snp_green[ch8,2]
g_ch8<-gmat[ch8,]
p_ch8<-apply(g_ch8,1,mean)/2

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
trt<-read.table("sort_green_phfull_tchum2019.txt",header=TRUE,comment.char="%",sep=",")

nas<-which(is.na(trt$Survival)==TRUE)

rmat<-matrix(regs[1,6:(6+3*2-1)],byrow=TRUE,ncol=2)
pcc<-LDChangePCverbose(regions=rmat,G=g_ch8[,-nas],pos=pos,trt=trt[-nas,])

t_ac<-which(trt$Treatment=="AC")
t_mm<-which(trt$Treatment=="MM")

t_ac_surv<-which(trt$Treatment=="AC" & trt$Survival==1)
t_mm_surv<-which(trt$Treatment=="MM" & trt$Survival==1)

### figure 5 manuscript #####
pdf("fig_PCA_comb.pdf",width=9.6,height=12.8)
layout(matrix(c(1,2,7,8,3,4,9,10,5,6,11,12),nrow=4,ncol=3,byrow=FALSE),widths=c(3,3,3),heights=c(3,3,3,3))
par(mar=c(4.5,5.5,2.5,1.5))
##########
lx<-min(Dac$x1)
ux<-max(Dac$x1)
ly<-min(Dac$x2)
uy<-max(Dac$x2)
sRG<-(c(lx,ux)*psds[1] + pmns[1])
sGB<-(c(ly,uy)*psds[2] + pmns[2])

betas<-est_hcl_ac[1:6,3]
ooo<-gridsurf2(betas,yy=seq(ly,uy,.02),xx=seq(lx,ux,.02))
image(t(ooo),axes=FALSE,col=alpha(rev(heat.colors(n=20)),1),breaks=seq(0,1,length.out=21),xlab="RG",ylab="GB",cex.lab=cl)
#polygon(c(.55,1.01,1.01,.55),c(1.01,1.01,.65,1.01),border=NA,col="white")
#polygon(c(.6,1.01,1.01,.6),c(1.01,1.01,.7,1.01),border=NA,col="white")
axis(1,seq(0,1,length.out=6),round(seq(sRG[1],sRG[2],length.out=6),2))
axis(2,seq(0,1,length.out=6),round(seq(sGB[1],sGB[2],length.out=6),2))
box()
cc<-rev(heat.colors(n=20))
cs<-rep(cc[1],length(Dac$x1))
bks<-seq(0,1,length.out=21)
ff<-calcfit(betas,Dac$x1,Dac$x2)
for(k in 2:20){
	cs[ff > bks[k]]<-cc[k]
	}
pxx<-(Dac$x1-lx)/(ux-lx)	
pyy<-(Dac$x2-ly)/(uy-ly)
#points(pxx,pyy,col="black",pch=21,bg=cs,cex=1.3)
title(main="(a) Fitness surface A/C",cex.main=cm)

betas<-est_hcl_mm[1:6,3]
lx<-min(Dmm$x1)
ux<-max(Dmm$x1)
ly<-min(Dmm$x2)
uy<-max(Dmm$x2)
ooo<-gridsurf2(betas,yy=seq(ly,uy,.02),xx=seq(lx,ux,.02))
image(t(ooo),axes=FALSE,col=alpha(rev(heat.colors(n=20)),1),breaks=seq(0,1,length.out=21),xlab="RG",ylab="GB",cex.lab=cl)
#polygon(c(.55,1.01,1.01,.55),c(1.01,1.01,.65,1.01),border=NA,col="white")
#polygon(c(.6,1.01,1.01,.6),c(1.01,1.01,.7,1.01),border=NA,col="white")
axis(1,seq(0,1,length.out=6),round(seq(sRG[1],sRG[2],length.out=6),2))
axis(2,seq(0,1,length.out=6),round(seq(sGB[1],sGB[2],length.out=6),2))
box()
cc<-rev(heat.colors(n=20))
cs<-rep(cc[1],length(Dmm$x1))
bks<-seq(0,1,length.out=21)
ff<-calcfit(betas,Dmm$x1,Dmm$x2)
for(k in 2:20){
	cs[ff > bks[k]]<-cc[k]
	}
pxx<-(Dmm$x1-lx)/(ux-lx)	
pyy<-(Dmm$x2-ly)/(uy-ly)
#points(pxx,pyy,col="black",pch=21,bg=cs,cex=1.3)
title(main="(d) Fitness surface MM",cex.main=cm)

AC<-which(phs$Treatment=="AC")
css1<-phs$hex[ac]
css2<-phs$hex[ac]
css1[phs$Survival[ac]==1]<-"black"
css2[phs$Survival[ac]==0]<-alpha("white",.3)
plot(phs$RG[ac],phs$GB[ac],pch=21,col=css1,bg=css2,xlab="RG",ylab="GB",cex.lab=cl,cex.axis=ca)
title(main="(b) Colour and survival on AC",cex.main=cm)

MM<-which(phs$Treatment=="MM")
css1<-phs$hex[mm]
css2<-phs$hex[mm]
css1[phs$Survival[mm]==1]<-"black"
css2[phs$Survival[mm]==0]<-alpha("white",.3)
plot(phs$RG[mm],phs$GB[mm],pch=21,col=css1,bg=css2,cex=1.2,xlab="RG",ylab="GB",cex.lab=cl,cex.axis=ca)
title(main="(e) Colour and survival on MM",cex.main=cm)


dd<-density(pcc$null_acs,adj=1.3)

plot(dd$x,dd$y,type='n',xlab="Mean PC correlation",ylab="Density",cex.lab=cl,cex.axis=ca)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),border=NA,col="gray80")
abline(v=pcc$obs_acs,col="firebrick4",lwd=1.5)
res_green[1,3]
#[1] 0.001998002
mtext("P = 0.002",3,line=-2,adj=.1,cex=1.2)
text(mean(pcc$null_acs),6,"null",cex=1.8)
text(pcc$obs_acs*.9375,10,"obs.",cex=1.8)
title(main="(c) Null PC correlations A/C",cex.main=cm)

dd<-density(pcc$null_mms,adj=1.3)

plot(dd$x,dd$y,type='n',xlab="Mean PC correlation",ylab="Density",cex.lab=cl,cex.axis=ca)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),border=NA,col="gray80")
abline(v=pcc$obs_mms,col="firebrick4",lwd=1.5)
res_green[1,10]
#[1] 0.3626374
mtext("P = 0.363",3,line=-2,adj=.1,cex=1.2)
title(main="(f) Null PC correlations MM",cex.main=cm)


##
pcs<-cbind(o1$x[,1],o2$x[,1],o3$x[,1],o1$x[,2],o2$x[,2],o3$x[,2])
rr<-abs(cor(pcs[t_ac,]))
rr[upper.tri(rr)]<-NA
diag(rr)<-NA
image(rr,breaks=c(-.1,1/6,2/6,3/6,4/6,5/6,1.1),col=heat.colors(n=6,rev=TRUE),axes=FALSE)
#image(rr,breaks=c(-.1,.2,.4,.6,.8,1.1),col=heat.colors(n=5,rev=TRUE),axes=FALSE)
polygon(c(.5,1.1,1.1,.5),c(-0.1,-0.1,.5,.5),col=NA)
aaa<-seq(0,1,.2)
text(aaa,aaa+.015,c("R1","R2","R3","R1","R2","R3"))
text(aaa,aaa-.04,c("PC1","PC1","PC1","PC2","PC2","PC2"))
#text(aaa,aaa,c("R1PC1","R2PC1","R3PC1","R1PC2","R2PC2","R3PC2"))
legend(-.13,1.125,c("0.0, 0.2","0.2, 0.3","0.3, 0.5","0.5, 0.7","0.7, 0.8","0.8, 1.0"),fill=heat.colors(n=6,rev=TRUE),bty='n',title="|r|",cex=1.15)
#legend(-.13,1.125,c("0.0, 0.2","0.2, 0.4","0.4, 0.6","0.6, 0.8","0.8, 1.0"),fill=heat.colors(n=5,rev=TRUE),bty='n',title="|r|",cex=1.15)
title(main="(g) Release PC correlations A/C",cex.main=cm)

pcs<-cbind(o1$x[,1],o2$x[,1],o3$x[,1],o1$x[,2],o2$x[,2],o3$x[,2])
rr<-abs(cor(pcs[t_mm,]))
rr[upper.tri(rr)]<-NA
diag(rr)<-NA
image(rr,breaks=c(-.1,1/6,2/6,3/6,4/6,5/6,1.1),col=heat.colors(n=6,rev=TRUE),axes=FALSE)
#image(rr,breaks=c(-.1,.2,.4,.6,.8,1.1),col=heat.colors(n=5,rev=TRUE),axes=FALSE)
polygon(c(.5,1.1,1.1,.5),c(-0.1,-0.1,.5,.5),col=NA)
aaa<-seq(0,1,.2)
text(aaa,aaa+.015,c("R1","R2","R3","R1","R2","R3"))
text(aaa,aaa-.04,c("PC1","PC1","PC1","PC2","PC2","PC2"))
#text(aaa,aaa,c("R1PC1","R2PC1","R3PC1","R1PC2","R2PC2","R3PC2"))
title(main="(j) Release PC correlations MM",cex.main=cm)

pcs<-cbind(o1$x[,1],o2$x[,1],o3$x[,1],o1$x[,2],o2$x[,2],o3$x[,2])
rr<-abs(cor(pcs[t_ac_surv,]))
rr[upper.tri(rr)]<-NA
diag(rr)<-NA
image(rr,breaks=c(-.1,1/6,2/6,3/6,4/6,5/6,1.1),col=heat.colors(n=6,rev=TRUE),axes=FALSE)
#image(rr,breaks=c(-.1,.2,.4,.6,.8,1.1),col=heat.colors(n=5,rev=TRUE),axes=FALSE)
polygon(c(.5,1.1,1.1,.5),c(-0.1,-0.1,.5,.5),col=NA)
aaa<-seq(0,1,.2)
text(aaa,aaa+.015,c("R1","R2","R3","R1","R2","R3"))
text(aaa,aaa-.04,c("PC1","PC1","PC1","PC2","PC2","PC2"))
#text(aaa,aaa,c("R1PC1","R2PC1","R3PC1","R1PC2","R2PC2","R3PC2"))
title(main="(h) Recapture PC correlations A/C",cex.main=cm)

pcs<-cbind(o1$x[,1],o2$x[,1],o3$x[,1],o1$x[,2],o2$x[,2],o3$x[,2])
rr<-abs(cor(pcs[t_mm_surv,]))
rr[upper.tri(rr)]<-NA
diag(rr)<-NA
image(rr,breaks=c(-.1,1/6,2/6,3/6,4/6,5/6,1.1),col=heat.colors(n=6,rev=TRUE),axes=FALSE)
#image(rr,breaks=c(-.1,.2,.4,.6,.8,1.1),col=heat.colors(n=5,rev=TRUE),axes=FALSE)
polygon(c(.5,1.1,1.1,.5),c(-0.1,-0.1,.5,.5),col=NA)
aaa<-seq(0,1,.2)
text(aaa,aaa+.015,c("R1","R2","R3","R1","R2","R3"))
text(aaa,aaa-.04,c("PC1","PC1","PC1","PC2","PC2","PC2"))
#text(aaa,aaa,c("R1PC1","R2PC1","R3PC1","R1PC2","R2PC2","R3PC2"))
title(main="(k) Recapture PC correlations MM",cex.main=cm)

pcs<-cbind(o1$x[,1],o2$x[,1],o3$x[,1],o1$x[,2],o2$x[,2],o3$x[,2])
rr<-abs(cor(pcs[t_ac_surv,]))-abs(cor(pcs[t_ac,]))
rr[upper.tri(rr)]<-NA
diag(rr)<-NA
image(rr,breaks=c(-.3,-.15,-.05,0,.05,.15,.3),col=rev(brewer.pal(7,"RdBu")[-4]),axes=FALSE)
polygon(c(.5,1.1,1.1,.5),c(-0.1,-0.1,.5,.5),col=NA)
aaa<-seq(0,1,.2)
#text(aaa,aaa,c("R1PC1","R2PC1","R3PC1","R1PC2","R2PC2","R3PC2"))
text(aaa,aaa+.015,c("R1","R2","R3","R1","R2","R3"))
text(aaa,aaa-.04,c("PC1","PC1","PC1","PC2","PC2","PC2"))
legend(-.13,1.125,c("-0.30, -0.15","-0.15, -0.05","-0.05, 0.00","0.00, 0.05","0.05, 0.15","0.15, 0.30"),fill=rev(brewer.pal(7,"RdBu")[-4]),bty='n',title="change in |r|",cex=1.3)
title(main="(i) Change PC correlations A/C",cex.main=cm)

pcs<-cbind(o1$x[,1],o2$x[,1],o3$x[,1],o1$x[,2],o2$x[,2],o3$x[,2])
rr<-abs(cor(pcs[t_mm_surv,]))-abs(cor(pcs[t_mm,]))
rr[upper.tri(rr)]<-NA
diag(rr)<-NA
image(rr,breaks=c(-.3,-.15,-.05,0,.05,.15,.3),col=rev(brewer.pal(7,"RdBu")[-4]),axes=FALSE)
polygon(c(.5,1.1,1.1,.5),c(-0.1,-0.1,.5,.5),col=NA)
aaa<-seq(0,1,.2)
#text(aaa,aaa,c("R1PC1","R2PC1","R3PC1","R1PC2","R2PC2","R3PC2"))
text(aaa,aaa+.015,c("R1","R2","R3","R1","R2","R3"))
text(aaa,aaa-.04,c("PC1","PC1","PC1","PC2","PC2","PC2"))
title(main="(l) Change PC correlations MM",cex.main=cm)

dev.off()
