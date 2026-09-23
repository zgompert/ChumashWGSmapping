# genomic predictions of phenotypes
library(data.table)
#library(rjags)
library(scales)

# Read command line arguments, pip threshold
# also used for file names
args <- commandArgs(trailingOnly = TRUE)
# Check if arguments are provided
if (length(args) == 0) {
  stop("No arguments provided")
}

pip<-as.numeric(args[1])
cat(pip,"\n")

## sorts effects to further increase LD 

Gf<-fread("ch8_mlpntest_tchum_green.geno",header=FALSE,sep=",")
ph<-read.table("ph_tchum_green.txt",header=FALSE)
pos<-as.numeric(unlist(Gf[,1]))
keepInd<-which(is.na(ph$V1)==FALSE)
gfmat<-as.matrix(Gf[,-c(1:3)])
gfmat<-gfmat[,keepInd]

## remove genotypic means to match data processing by gemma
center<-function(x=NA){x<-x-mean(x);return(x)}
gcent<-apply(gfmat,1,center)

## get effects, RG
rg<-fread("output/mav_o_ch8_poly_bslmm_green_ph1_ch0.param.txt",header=FALSE)
betaRG<-as.numeric(rg$V2)
epos<-as.numeric(unlist(rg$V1))
keep<-which(pos %in% epos)

gfmat<-gfmat[keep,]
gcent<-gcent[,keep]
L<-length(keep)

ph<-ph[keepInd,]
## compare to observed
pscoreRG<-betaRG %*% t(gcent)
pscoreRG<-pscoreRG+mean(ph[,1],na.rm=TRUE) ## add back mean to make comparable
cor.test(ph[,1],as.numeric(pscoreRG))

summary(lm(ph[,1]~as.numeric(pscoreRG)))

## get effects, GB
gb<-fread("output/mav_o_ch8_poly_bslmm_green_ph2_ch0.param.txt",header=FALSE)
betaGB<-as.numeric(gb$V2)
#
## compare to observed
pscoreGB<-betaGB %*% t(gcent)
pscoreGB<-pscoreGB+mean(ph[,2],na.rm=TRUE) ## add back mean to make comparable
cor.test(ph[,2],as.numeric(pscoreGB))

summary(lm(ph[,2]~as.numeric(pscoreGB)))



phi<-apply(gfmat,1,mean)/2

## get the pips
pipRG<-fread("output/pips_o_ch8_poly_bslmm_green_ph1_ch0.param.txt",header=FALSE)
pipGB<-fread("output/pips_o_ch8_poly_bslmm_green_ph2_ch0.param.txt",header=FALSE)

## what about Va
vaRG<-2*phi*(1-phi)*pipRG$V2^2
vaGB<-2*phi*(1-phi)*pipGB$V2^2


hipipRG<-which(pipRG$V2 > pip)
hipipGB<-which(pipGB$V2 > pip)

## no overlap at pip 0.05, 0.075 or 0.1

hipip<-c(hipipRG,hipipGB)
phipip<-phi[hipip]
hipipBrg<-betaRG[hipipRG]
hipipBgb<-betaGB[hipipGB]
hipipB<-c(hipipBrg,-1*hipipBgb)## use negative so that "positive" is melanic

## denote RG or GB, 1 or 2
waxis<-c(rep(1,length(hipipRG)),rep(2,length(hipipGB)))

N<-dim(gfmat)[2]

## maximize LD
## define groups
L<-length(hipip)
ko<-kmeans(x=ph[,1:2],center=2,iter.max=20,nstart=20)
if(sum(ko$cluster==1) > sum(ko$cluster==2)){ ## flip the assignments
	ko$cluster<-3-ko$cluster
}
plot(ph$V1,ph$V2,pch=19,col=c("brown","green")[ko$cluster])

grp<-rep(NA,dim(ph)[1])
grp<-ko$cluster
p1<-apply(gfmat[hipip,which(grp==1)],1,mean)/2
p2<-apply(gfmat[hipip,which(grp==2)],1,mean)/2
sz1<-mean(grp==1,na.rm=TRUE)
sz2<-mean(grp==2,na.rm=TRUE)
pmax1<-rep(NA,L)
pmax2<-rep(NA,L)
for(k in 1:L){
	if(hipipB[k] > 0){  
		if(phipip[k] < sz1){ ## all goes in group 1
			pmax1[k]<-phipip[k]/sz1
			pmax2[k]<-0
		} else { ## some left over
			pmax1[k]<-1
			pmax2[k]<-(phipip[k]-sz1)/sz2
		} 
	} else{
                if(phipip[k] < sz2){ ## all goes in group 2
                        pmax2[k]<-phipip[k]/sz2
                        pmax1[k]<-0
                } else { ## some left over
                        pmax2[k]<-1
                        pmax1[k]<-(phipip[k]-sz2)/sz1
                }
	}
}
mean(abs(p1-p2))
mean(abs(pmax1-pmax2))


ld<-cor(gcent[,hipip])^2
ldobs<-mean(ld[upper.tri(ld)])

N<-dim(gfmat)[2]
## loop over radnomization = no LD
rgRan<-matrix(NA,nrow=N,ncol=50)
gbRan<-matrix(NA,nrow=N,ncol=50)
ldRan<-rep(NA,50)
gran<-gcent
for(i in 1:50){
	for(k in 1:L){
		y<-sum(gfmat[hipip[k],]);n<-(N*2)-y
		sg<-apply(matrix(sample(c(rep(1,y),rep(0,n)),size=2*N,replace=FALSE),ncol=2),1,sum)
		gran[,hipip[k]]<-sg-mean(gfmat[hipip[k],])
	}
	tempRG<-betaRG %*% t(gran)
	rgRan[,i]<-tempRG+mean(ph[,1],na.rm=TRUE)
	tempGB<-betaGB %*% t(gran)
	gbRan[,i]<-tempGB+mean(ph[,2],na.rm=TRUE)
	ld<-cor(gran[,hipip])^2
	ldRan[i]<-mean(ld[upper.tri(ld)])
}


## loop over anti-radnomization = maxing LD
rgMax<-matrix(NA,nrow=N,ncol=50)
gbMax<-matrix(NA,nrow=N,ncol=50)
ldMax<-rep(NA,50)
gmax<-gcent

for(i in 1:50){
	for(k in 1:L){
		pp<-as.integer(round(pmax1[k]*sz1*N*2));qq<-as.integer(round((1-pmax1[k])*sz1*N*2))
		sgG1<-sort(apply(matrix(sample(c(rep(1,pp),rep(0,qq)),size=as.integer(round(sz1*N*2)),replace=FALSE),ncol=2),1,sum))
		pp<-as.integer(round(pmax2[k]*sz2*N*2));qq<-as.integer(round((1-pmax2[k])*sz2*N*2))
		sgG2<-sort(apply(matrix(sample(c(rep(1,pp),rep(0,qq)),size=as.integer(round(sz2*N*2)),replace=FALSE),ncol=2),1,sum))
		sgG1<-sgG1-mean(gfmat[hipip[k],])
		sgG2<-sgG2-mean(gfmat[hipip[k],])
		gtemp<-gcent[,hipip[k]]
		gbBase<-jitter(gtemp*betaGB[hipip[k]])
		rgBase<-jitter(gtemp*betaRG[hipip[k]])
		if(waxis[k]==1 & hipipB[k] > 0){# RG SNP, grp 1 higher, effect pos
			gmax[which(grp==1),hipip[k]]<-sgG1[order(order(rgBase[which(grp==1)]))]
			gmax[which(grp==2),hipip[k]]<-sgG2[order(order(rgBase[which(grp==2)]))]
		} else if(waxis[k]==1 & hipipB[k] <= 0){ ## RG SNP, grp 1 higher, effect neg
			gmax[which(grp==1),hipip[k]]<-sgG1[order(rev(order(rgBase[which(grp==1)])))]
			gmax[which(grp==2),hipip[k]]<-sgG2[order(rev(order(rgBase[which(grp==2)])))]
		} else if(waxis[k]==2 & hipipB[k] > 0){# GB SNP, grp 2 higher, effect neg
                        gmax[which(grp==1),hipip[k]]<-sgG1[order(rev(order(gbBase[which(grp==1)])))]
                        gmax[which(grp==2),hipip[k]]<-sgG2[order(rev(order(gbBase[which(grp==2)])))]
		} else{ # GB SNP, grp 1 higher, effect pos
                        gmax[which(grp==1),hipip[k]]<-sgG1[order(order(gbBase[which(grp==1)]))]
                        gmax[which(grp==2),hipip[k]]<-sgG2[order(order(gbBase[which(grp==2)]))]
		}

	}
	tempRG<-betaRG %*% t(gmax)
	rgMax[,i]<-tempRG+mean(ph[,1],na.rm=TRUE)
	tempGB<-betaGB %*% t(gmax)
	gbMax[,i]<-tempGB+mean(ph[,2],na.rm=TRUE)
	ld<-cor(gmax[,hipip])^2
	ldMax[i]<-mean(ld[upper.tri(ld)])
}	


## get selection data
load("/uufs/chpc.utah.edu/common/home/gompert-group4/projects/timema_genarch/color_exp_2019/bayes_color_2019.rdat")
# selection coefficients from selection gradient, using alphas, hierarchical level effects
bMM<-est_hcl_mm[1:6,3]
bAC<-est_hcl_ac[1:6,3]
est_hcl_ac[1:6,]
#               2.5%        25%        50%        75%       97.5%
#alpha[1] -2.3032038 -1.4616983 -1.1175156 -0.7782439 -0.06763818
#alpha[2] -0.8420637  0.1809069  0.4940390  0.8370925  2.40637623
#alpha[3] -1.2436702  0.1479737  0.4571522  0.7292981  1.63911571
#alpha[4] -2.2058414 -0.4798853 -0.1369185  0.1542778  1.27441568
#alpha[5] -2.6472714 -1.4042134 -0.9599061 -0.5329878  0.99378836
#alpha[6] -3.4294349 -2.8363119 -2.5581673 -2.2846823 -1.71081218
est_hcl_mm[1:6,]
#               2.5%         25%         50%        75%      97.5%
#alpha[1] -1.1124779 -0.21068750  0.08981357  0.4013856  1.4471706
#alpha[2] -0.6450221 -0.18926332 -0.00817932  0.1884889  0.6599292
#alpha[3] -0.2601696  0.17985564  0.36269085  0.5542720  0.9423904
#alpha[4] -0.3181178  0.08062346  0.25099129  0.4340146  0.8467638
#alpha[5] -0.4705779  0.17011384  0.43190912  0.6952102  1.3068882
#alpha[6] -2.7900761 -2.13395724 -1.87939975 -1.6238560 -1.0600209


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
sd(pgAC)
mfAC<-mean(pgAC)
mfMM<-mean(pgMM)

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


library(scales)
ff<-paste("relativeFitnessRanCh8sort_green_",pip,".pdf",sep="")
pdf(ff,width=6,height=5)

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
mean(relf_ran_AC > 1)
mean(relf_max_MM > 1)
mean(relf_max_AC > 1)
mean(relf_ext_MM > 1)
mean(relf_ext_AC > 1)

## median, 95ETPI and 90ETIP
quantile(relf_ran_MM,probs=c(.5,.025,.975,.05,.95))
quantile(relf_ran_AC,probs=c(.5,.025,.975,.05,.95))
quantile(relf_max_MM,probs=c(.5,.025,.975,.05,.95))
quantile(relf_max_AC,probs=c(.5,.025,.975,.05,.95))
quantile(relf_ext_MM,probs=c(.5,.025,.975,.05,.95))
quantile(relf_ext_AC,probs=c(.5,.025,.975,.05,.95))

library(scales)
cl<-1.6;cm<-1.6;ca<-1.1
ff<-paste("relativeFitnessUncCh8sort_green_",pip,".pdf",sep="")
pdf(ff,width=6,height=11)
par(mfrow=c(3,1))
par(mar=c(5,4.5,2.5,1))
plot(density(relf_ran_MM,adj=1.2),xlim=c(.6,4),ylim=c(0,5),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dmm<-density(relf_ran_MM,adj=1.2)
polygon(c(dmm$x,rev(dmm$x)),c(rep(0,length(dmm$y)),rev(dmm$y)),col=alpha("cadetblue",.5),border=NA)
dac<-density(relf_ran_AC,adj=1.2)
polygon(c(dac$x,rev(dac$x)),c(rep(0,length(dac$y)),rev(dac$y)),col=alpha("orange",.5),border=NA)
abline(v=1,lty=3)
legend(3.5,5,c("MM","A/C"),fill=c(alpha("cadetblue",.5),alpha("orange",.5)),bty='n',cex=1.5)
mtext(paste("MM prob. > 1 = ",round(mean(relf_ran_MM > 1),2),sep=""),1,line=-8,adj=.85,cex=ca)
mtext(paste("A/C prob. > 1 = ",round(mean(relf_ran_AC > 1),2),sep=""),1,line=-6,adj=.85,cex=ca)
title(main="(a) Observed vs minimum LD",cex.main=cm)

plot(density(relf_max_MM,adj=1.2),xlim=c(.6,4),ylim=c(0,4),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl)
dmm<-density(relf_max_MM,adj=1.2)
polygon(c(dmm$x,rev(dmm$x)),c(rep(0,length(dmm$y)),rev(dmm$y)),col=alpha("cadetblue",.5),border=NA)
dac<-density(relf_max_AC,adj=1.2)
polygon(c(dac$x,rev(dac$x)),c(rep(0,length(dac$y)),rev(dac$y)),col=alpha("orange",.5),border=NA)
abline(v=1,lty=3)
mtext(paste("MM prob. > 1 = ",round(mean(relf_max_MM > 1),2),sep=""),1,line=-8,adj=.85,cex=ca)
mtext(paste("A/C prob. > 1 = ",round(mean(relf_max_AC > 1),2),sep=""),1,line=-6,adj=.85,cex=ca)
title(main="(b) Maximum vs observed LD",cex.main=cm)

plot(density(relf_ext_MM,adj=1.2),xlim=c(.6,4),ylim=c(0,2),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl)
dmm<-density(relf_ext_MM,adj=1.2)
polygon(c(dmm$x,rev(dmm$x)),c(rep(0,length(dmm$y)),rev(dmm$y)),col=alpha("cadetblue",.5),border=NA)
dac<-density(relf_ext_AC,adj=1.2)
polygon(c(dac$x,rev(dac$x)),c(rep(0,length(dac$y)),rev(dac$y)),col=alpha("orange",.5),border=NA)
abline(v=1,lty=3)
mtext(paste("MM prob. > 1 = ",round(mean(relf_ext_MM > 1),2),sep=""),1,line=-8,adj=.85,cex=ca)
mtext(paste("A/C prob. > 1 = ",round(mean(relf_ext_AC > 1),2),sep=""),1,line=-6,adj=.85,cex=ca)
title(main="(c) Maximum vs minimum LD",cex.main=cm)
dev.off()
rm(list=c("gcent","Gf","gfmat","gmax"))

ff<-paste("predictCh8_sort_green_",pip,".rdat",sep="")
save(list=ls(),file=ff)

