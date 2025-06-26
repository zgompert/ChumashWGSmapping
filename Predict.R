## genomic predictions of phenotypes
library(data.table)
library(rjags)
library(scales)

G<-fread("hsub_mlpntest_EXP_tchum.geno",header=FALSE,sep=",")
ph<-read.table("ph_tchum2019.txt",header=FALSE)
gmat<-as.matrix(G[,-c(1:3)])

## remove genotypic means to match data processing by gemma
center<-function(x=NA){x<-x-mean(x);return(x)}
gcent<-apply(gmat,1,center)

## get effects, RG
rg<-fread("output/mav_o_poly_bslmm_ph1_ch0.param.txt",header=FALSE)
betaRG<-as.numeric(rg$V2)

## compare to observed
pscoreRG<-betaRG %*% t(gcent)
pscoreRG<-pscoreRG+mean(ph[,1],na.rm=TRUE) ## add back mean to make comparable
cor.test(ph[,1],as.numeric(pscoreRG))
#t = 84.912, df = 422, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9661573 0.9767798
#sample estimates:
#      cor 
#0.9719602 

## slope shallower, 0.745, but no bias

## get effects, GB
gb<-fread("output/mav_o_poly_bslmm_ph2_ch0.param.txt",header=FALSE)
betaGB<-as.numeric(gb$V2)

## compare to observed
pscoreGB<-betaGB %*% t(gcent)
pscoreGB<-pscoreGB+mean(ph[,2],na.rm=TRUE) ## add back mean to make comparable
cor.test(ph[,2],as.numeric(pscoreGB))
#data:  ph[, 2] and as.numeric(pscoreGB)
#t = 50.268, df = 422, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9107393 0.9382098
#sample estimates:
#      cor
#0.9256854

## again slope shallower but little bias (intercept not 0 though in this case, still small).

## get the pips
pipRG<-fread("output/pips_o_poly_bslmm_ph1_ch0.param.txt",header=FALSE)
pipGB<-fread("output/pips_o_poly_bslmm_ph2_ch0.param.txt",header=FALSE)

## get pip > .1
hipip<-which(pipRG$V2 > .1 | pipGB$V2 > .1)
length(hipip)
#[1] 15
phipips<-apply(gmat[hipip,],1,mean)/2
N<-dim(gmat)[2]
## loop over radnomization = no LD
rgRan<-matrix(NA,nrow=N,ncol=50)
gbRan<-matrix(NA,nrow=N,ncol=50)
for(i in 1:50){
	gran<-gcent
	for(k in 1:length(hipip)){
		gran[,hipip[k]]<-rbinom(n=N,size=2,prob=phipips[k])-mean(gmat[hipip[k],])
	}
	tempRG<-betaRG %*% t(gran)
	rgRan[,i]<-tempRG+mean(ph[,1],na.rm=TRUE)
	tempGB<-betaGB %*% t(gran)
	gbRan[,i]<-tempGB+mean(ph[,2],na.rm=TRUE)
}	

## maximize LD
## define groups
nona<-which(is.na(ph[,1])==FALSE)
ko<-kmeans(x=ph[nona,1:2],center=2,iter.max=20,nstart=20)
plot(ph$V1[nona],ph$V2[nona],pch=19,col=c("brown","green")[ko$cluster])
grp<-rep(NA,dim(ph)[1])
grp[nona]<-ko$cluster
p1<-apply(gmat[hipip,which(grp==1)],1,mean)/2
p2<-apply(gmat[hipip,which(grp==2)],1,mean)/2
sz1<-mean(grp==1,na.rm=TRUE)
sz2<-mean(grp==2,na.rm=TRUE)
pmax1<-rep(NA,length(hipip))
pmax2<-rep(NA,length(hipip))
for(k in 1:length(hipip)){
	if(p1[k] > p2[k]){
		if(phipips[k] < sz1){ ## all goes in group 1
			pmax1[k]<-phipips[k]/sz1
			pmax2[k]<-0
		} else { ## some leftofer
			pmax1[k]<-1
			pmax2[k]<-(phipips[k]-sz1)/sz2
		} 
	} else{
                if(phipips[k] < sz2){ ## all goes in group 2
                        pmax2[k]<-phipips[k]/sz2
                        pmax1[k]<-0
                } else { ## some leftofer
                        pmax2[k]<-1
                        pmax1[k]<-(phipips[k]-sz2)/sz1
                }
	}
}
mean(abs(p1-p2))
#[1] 0.1613982
mean(abs(pmax1-pmax2))
#[1] 0.3390849


## get pip > .1
hipip<-which(pipRG$V2 > .1 | pipGB$V2 > .1)
length(hipip)
#[1] 15
phipips<-apply(gmat[hipip,],1,mean)/2
N<-dim(gmat)[2]
## loop over radnomization = no LD
rgRan<-matrix(NA,nrow=N,ncol=50)
gbRan<-matrix(NA,nrow=N,ncol=50)
for(i in 1:50){
	gran<-gcent
	for(k in 1:length(hipip)){
		gran[,hipip[k]]<-rbinom(n=N,size=2,prob=phipips[k])-mean(gmat[hipip[k],])
	}
	tempRG<-betaRG %*% t(gran)
	rgRan[,i]<-tempRG+mean(ph[,1],na.rm=TRUE)
	tempGB<-betaGB %*% t(gran)
	gbRan[,i]<-tempGB+mean(ph[,2],na.rm=TRUE)
}	

## loop over anti-radnomization = maxing LD
rgMax<-matrix(NA,nrow=N,ncol=50)
gbMax<-matrix(NA,nrow=N,ncol=50)
for(i in 1:50){
	gmax<-gcent
	for(k in 1:length(hipip)){
		pvec<-rep(phipips[k],N)
		pvec[which(grp==1)]<-pmax1[k]
		pvec[which(grp==2)]<-pmax2[k]
		gmax[,hipip[k]]<-rbinom(n=N,size=2,prob=pvec)-mean(gmat[hipip[k],])
	}
	tempRG<-betaRG %*% t(gmax)
	rgMax[,i]<-tempRG+mean(ph[,1],na.rm=TRUE)
	tempGB<-betaGB %*% t(gmax)
	gbMax[,i]<-tempGB+mean(ph[,2],na.rm=TRUE)
}	

save(list=ls(),file="predict.rdat")

## get selection data
load("/uufs/chpc.utah.edu/common/home/u6000989/../gompert-group4/projects/timema_genarch/color_exp_2019/bayes_color_2019.rdat")
# selection coefficients from selection gradient, using alphas, hierarchical level effects
bMM<-est_hcl_mm[1:6,3]
bAC<-est_hcl_ac[1:6,3]

# b6 + b1 * rg + b2 *gb + b3 * rg^2 + b4 * gb^2 + b5 * rg * gb
## for standardized trait values
pmns<-apply(ph,2,mean,na.rm=TRUE)
psds<-apply(ph,2,sd,na.rm=TRUE)

sRG<-(as.numeric(pscoreRG)-pmns[1])/psds[1]
sGB<-(as.numeric(pscoreGB)-pmns[2])/psds[2]

lpgMM<-bMM[6] + bMM[1] * sRG + bMM[2] * sGB + bMM[3] * sRG^2 + bMM[4] * sGB^2 + bMM[5] * sRG * sGB
pgMM<-1/(1+exp(-1 * lpgMM))

lpgAC<-bAC[6] + bAC[1] * sRG + bAC[2] * sGB + bAC[3] * sRG^2 + bAC[4] * sGB^2 + bAC[5] * sRG * sGB
pgAC<-1/(1+exp(-1 * lpgAC))

sd(pgMM)
#[1] 0.05216237
sd(pgAC)
#[1] 0.1151358
mfAC<-mean(pgAC)
#[1] 0.1358
mfMM<-mean(pgMM)
#[1] 0.162113

## random
ran_mfAC<-rep(NA,50)
ran_mfMM<-rep(NA,50)

for(i in 1:50){
	sRG<-(rgRan[,i]-pmns[1])/psds[1]
	sGB<-(gbRan[,i]-pmns[2])/psds[2]

	lpgMM<-bMM[6] + bMM[1] * sRG + bMM[2] * sGB + bMM[3] * sRG^2 + bMM[4] * sGB^2 + bMM[5] * sRG * sGB
	rpgMM<-1/(1+exp(-1 * lpgMM))

	lpgAC<-bAC[6] + bAC[1] * sRG + bAC[2] * sGB + bAC[3] * sRG^2 + bAC[4] * sGB^2 + bAC[5] * sRG * sGB
	rpgAC<-1/(1+exp(-1 * lpgAC))
	ran_mfMM[i]<-mean(rpgMM)
	ran_mfAC[i]<-mean(rpgAC)
}

## anti-random
max_mfAC<-rep(NA,50)
max_mfMM<-rep(NA,50)

for(i in 1:50){
	sRG<-(rgMax[,i]-pmns[1])/psds[1]
	sGB<-(gbMax[,i]-pmns[2])/psds[2]

	lpgMM<-bMM[6] + bMM[1] * sRG + bMM[2] * sGB + bMM[3] * sRG^2 + bMM[4] * sGB^2 + bMM[5] * sRG * sGB
	rpgMM<-1/(1+exp(-1 * lpgMM))

	lpgAC<-bAC[6] + bAC[1] * sRG + bAC[2] * sGB + bAC[3] * sRG^2 + bAC[4] * sGB^2 + bAC[5] * sRG * sGB
	rpgAC<-1/(1+exp(-1 * lpgAC))
	max_mfMM[i]<-mean(rpgMM)
	max_mfAC[i]<-mean(rpgAC)
}

save(list=ls(),file="predict.rdat")

library(scales)
pdf("relativeFitnessRan.pdf",width=6,height=5)

plot(density(mfMM/ran_mfMM,adj=1.2),xlim=c(.9,1.5),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=1.3)
dmm<-density(mfMM/ran_mfMM,adj=1.2)
polygon(c(dmm$x,rev(dmm$x)),c(rep(0,length(dmm$y)),dmm$y),col=alpha("cadetblue",.5),border=NA)
dac<-density(mfAC/ran_mfAC,adj=1.2)
polygon(c(dac$x,rev(dac$x)),c(rep(0,length(dac$y)),dac$y),col=alpha("orange",.5),border=NA)
abline(v=1,lty=3)
legend(1.4,45,c("MM","A/C"),fill=c(alpha("cadetblue",.5),alpha("orange",.5)),bty='n')
dev.off()

## posterior samples are here
dim(out_hcl_ac[[1]])
#1] 2000  478
length(out_hcl_ac)

## repeat over posterior, 1000 samples total, 20 for each of 50 reps

ran_mfAC<-matrix(NA,nrow=50,ncol=20)
ran_mfMM<-matrix(NA,nrow=50,ncol=20)
max_mfAC<-matrix(NA,nrow=50,ncol=20)
max_mfMM<-matrix(NA,nrow=50,ncol=20)
obs_mfAC<-matrix(NA,nrow=50,ncol=20)
obs_mfMM<-matrix(NA,nrow=50,ncol=20)

mcb_AC<-rbind(out_hcl_ac[[1]][,1:6],out_hcl_ac[[2]][,1:6],out_hcl_ac[[3]][,1:6])
mcb_MM<-rbind(out_hcl_mm[[1]][,1:6],out_hcl_mm[[2]][,1:6],out_hcl_mm[[3]][,1:6])
nMC<-dim(mcb_AC)[1] # 6000

for(i in 1:50){
	for(j in 1:20){
		x<-sample(1:nMC,1)
		## random
		sRG<-(rgRan[,i]-pmns[1])/psds[1]
		sGB<-(gbRan[,i]-pmns[2])/psds[2]

		lpgMM<-mcb_MM[x,6] + mcb_MM[x,1] * sRG + mcb_MM[x,2] * sGB + mcb_MM[x,3] * sRG^2 + 
			mcb_MM[x,4] * sGB^2 + mcb_MM[x,5] * sRG * sGB
		rpgMM<-1/(1+exp(-1 * lpgMM))

		lpgAC<-mcb_AC[x,6] + mcb_AC[x,1] * sRG + mcb_AC[x,2] * sGB + mcb_AC[x,3] * sRG^2 + 
			mcb_AC[x,4] * sGB^2 + mcb_AC[x,5] * sRG * sGB
		rpgAC<-1/(1+exp(-1 * lpgAC))
		ran_mfMM[i,j]<-mean(rpgMM)
		ran_mfAC[i,j]<-mean(rpgAC)

		## antirandom
		sRG<-(rgMax[,i]-pmns[1])/psds[1]
		sGB<-(gbMax[,i]-pmns[2])/psds[2]
	
		lpgMM<-mcb_MM[x,6] + mcb_MM[x,1] * sRG + mcb_MM[x,2] * sGB + mcb_MM[x,3] * sRG^2 +
		       	mcb_MM[x,4] * sGB^2 + mcb_MM[x,5] * sRG * sGB
		rpgMM<-1/(1+exp(-1 * lpgMM))

		lpgAC<-mcb_AC[x,6] + mcb_AC[x,1] * sRG + mcb_AC[x,2] * sGB + mcb_AC[x,3] * sRG^2 + 
			mcb_AC[x,4] * sGB^2 + mcb_AC[x,5] * sRG * sGB
		rpgAC<-1/(1+exp(-1 * lpgAC))
		max_mfMM[i,j]<-mean(rpgMM)
		max_mfAC[i,j]<-mean(rpgAC)

		## observed
		
		sRG<-(as.numeric(pscoreRG)-pmns[1])/psds[1]
		sGB<-(as.numeric(pscoreGB)-pmns[2])/psds[2]

		lpgMM<-mcb_MM[x,6] + mcb_MM[x,1] * sRG + mcb_MM[x,2] * sGB + mcb_MM[x,3] * sRG^2 + 
			mcb_MM[x,4] * sGB^2 + mcb_MM[x,5] * sRG * sGB
		pgMM<-1/(1+exp(-1 * lpgMM))

		lpgAC<-mcb_AC[x,6] + mcb_AC[x,1] * sRG + mcb_AC[x,2] * sGB + mcb_AC[x,3] * sRG^2 + 
			mcb_AC[x,4] * sGB^2 + mcb_AC[x,5] * sRG * sGB
		pgAC<-1/(1+exp(-1 * lpgAC))

		obs_mfMM[i,j]<-mean(pgMM)
		obs_mfAC[i,j]<-mean(pgAC)
	}
}

relf_ran_MM<-as.vector(obs_mfMM/ran_mfMM)
relf_ran_AC<-as.vector(obs_mfAC/ran_mfAC)
relf_max_MM<-as.vector(max_mfMM/obs_mfMM)
relf_max_AC<-as.vector(max_mfAC/obs_mfAC)
relf_ext_MM<-as.vector(max_mfMM/ran_mfMM)
relf_ext_AC<-as.vector(max_mfAC/ran_mfAC)

## prob > 1
mean(relf_ran_MM > 1)
#[1] 0.812
mean(relf_ran_AC > 1)
#[1] 0.942
mean(relf_max_MM > 1)
#[1] 0.815
mean(relf_max_AC > 1)
#[1] 0.945
mean(relf_ext_MM > 1)
#[1] 0.823
mean(relf_ext_AC > 1)
#[1] 0.984

## median, 95ETPI and 90ETIP
quantile(relf_ran_MM,probs=c(.5,.025,.975,.05,.95))
#      50%      2.5%     97.5%        5%       95% 
#1.0749412 0.9136294 1.2994537 0.9393364 1.2558112 
quantile(relf_ran_AC,probs=c(.5,.025,.975,.05,.95))
#      50%      2.5%     97.5%        5%       95% 
#1.2728799 0.9445857 1.6684976 0.9908804 1.5890200 
quantile(relf_max_MM,probs=c(.5,.025,.975,.05,.95))
#      50%      2.5%     97.5%        5%       95% 
#1.0970248 0.8779211 1.3301346 0.9123359 1.2935863 
quantile(relf_max_AC,probs=c(.5,.025,.975,.05,.95))
#      50%      2.5%     97.5%        5%       95% 
#1.2688999 0.9431728 1.6698536 0.9969875 1.6018711 
quantile(relf_ext_MM,probs=c(.5,.025,.975,.05,.95))
#      50%      2.5%     97.5%        5%       95% 
#1.1897078 0.8087499 1.6783626 0.8716734 1.5878912 
quantile(relf_ext_AC,probs=c(.5,.025,.975,.05,.95))
#     50%     2.5%    97.5%       5%      95% 
#1.595446 1.043744 2.416977 1.148429 2.274322 

library(scales)
cl<-1.6;cm<-1.6;ca<-1.1
pdf("relativeFitnessUnc.pdf",width=6,height=11)
par(mfrow=c(3,1))
par(mar=c(5,4.5,2.5,1))
plot(density(relf_ran_MM,adj=1.2),xlim=c(.6,4),ylim=c(0,5),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dmm<-density(relf_ran_MM,adj=1.2)
polygon(c(dmm$x,rev(dmm$x)),c(rep(0,length(dmm$y)),rev(dmm$y)),col=alpha("cadetblue",.5),border=NA)
dac<-density(relf_ran_AC,adj=1.2)
polygon(c(dac$x,rev(dac$x)),c(rep(0,length(dac$y)),rev(dac$y)),col=alpha("orange",.5),border=NA)
abline(v=1,lty=3)
legend(3.5,5,c("MM","A/C"),fill=c(alpha("cadetblue",.5),alpha("orange",.5)),bty='n',cex=1.5)
mtext("MM prob. > 1 = 0.82",1,line=-8,adj=.85,cex=ca)
mtext("A/C prob. > 1 = 0.94",1,line=-6,adj=.85,cex=ca)
title(main="(a) Observed vs minimum LD",cex.main=cm)

plot(density(relf_max_MM,adj=1.2),xlim=c(.6,4),ylim=c(0,4),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl)
dmm<-density(relf_max_MM,adj=1.2)
polygon(c(dmm$x,rev(dmm$x)),c(rep(0,length(dmm$y)),rev(dmm$y)),col=alpha("cadetblue",.5),border=NA)
dac<-density(relf_max_AC,adj=1.2)
polygon(c(dac$x,rev(dac$x)),c(rep(0,length(dac$y)),rev(dac$y)),col=alpha("orange",.5),border=NA)
abline(v=1,lty=3)
mtext("MM prob. > 1 = 0.82",1,line=-8,adj=.85,cex=ca)
mtext("A/C prob. > 1 = 0.95",1,line=-6,adj=.85,cex=ca)
title(main="(b) Maximum vs observed LD",cex.main=cm)

plot(density(relf_ext_MM,adj=1.2),xlim=c(.6,4),ylim=c(0,2),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl)
dmm<-density(relf_ext_MM,adj=1.2)
polygon(c(dmm$x,rev(dmm$x)),c(rep(0,length(dmm$y)),rev(dmm$y)),col=alpha("cadetblue",.5),border=NA)
dac<-density(relf_ext_AC,adj=1.2)
polygon(c(dac$x,rev(dac$x)),c(rep(0,length(dac$y)),rev(dac$y)),col=alpha("orange",.5),border=NA)
abline(v=1,lty=3)
mtext("MM prob. > 1 = 0.82",1,line=-8,adj=.85,cex=ca)
mtext("A/C prob. > 1 = 0.98",1,line=-6,adj=.85,cex=ca)
title(main="(c) Maximum vs minimum LD",cex.main=cm)
dev.off()
save(list=ls(),file="predict.rdat")

