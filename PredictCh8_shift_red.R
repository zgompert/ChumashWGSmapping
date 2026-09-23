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



## shifts means 

Gf<-fread("ch8_mlpntest_tchum_red.geno",header=FALSE,sep=",")
ph<-read.table("ph_tchum_red.txt",header=FALSE)
gfmat<-as.matrix(Gf[,-c(1:3)])
pos<-as.numeric(unlist(Gf[,1]))

## remove genotypic means to match data processing by gemma
center<-function(x=NA){x<-x-mean(x);return(x)}
gcent<-apply(gfmat,1,center)

## get effects, RG
rg<-fread("output/mav_o_ch8_poly_bslmm_red_ph1_ch0.param.txt",header=FALSE)
betaRG<-as.numeric(rg$V2)
epos<-as.numeric(unlist(rg$V1))
keep<-which(pos %in% epos)

gfmat<-gfmat[keep,]
gcent<-gcent[,keep]
L<-length(keep)

## compare to observed
pscoreRG<-betaRG %*% t(gcent)
pscoreRG<-pscoreRG+mean(ph[,1],na.rm=TRUE) ## add back mean to make comparable
cor.test(ph[,1],as.numeric(pscoreRG))

summary(lm(ph[,1]~as.numeric(pscoreRG)))

## get effects, GB
gb<-fread("output/mav_o_ch8_poly_bslmm_red_ph2_ch0.param.txt",header=FALSE)
betaGB<-as.numeric(gb$V2)
#
## compare to observed
pscoreGB<-betaGB %*% t(gcent)
pscoreGB<-pscoreGB+mean(ph[,2],na.rm=TRUE) ## add back mean to make comparable
cor.test(ph[,2],as.numeric(pscoreGB))

summary(lm(ph[,2]~as.numeric(pscoreGB)))

phi<-apply(gfmat,1,mean)/2

## get the pips
pipRG<-fread("output/pips_o_ch8_poly_bslmm_red_ph1_ch0.param.txt",header=FALSE)
pipGB<-fread("output/pips_o_ch8_poly_bslmm_red_ph2_ch0.param.txt",header=FALSE)

## what about Va
vaRG<-2*phi*(1-phi)*pipRG$V2^2
vaGB<-2*phi*(1-phi)*pipGB$V2^2


hipipRG<-which(pipRG$V2 > pip)
hipipGB<-which(pipGB$V2 > pip)

## check for ovelap
overlap<-hipipRG[hipipRG %in% hipipGB]
if(length(overlap)>0){ ## some overlap
	for(k in 1:length(overlap)){
		if(pipRG$V2[overlap[k]] > pipGB$V2[overlap[k]]){
			a<-which(hipipGB==overlap[k])
			hipipGB<-hipipGB[-a]
		} else{
                        a<-which(hipipRG==overlap[k])
                        hipipRG<-hipipRG[-a]
		}
	}
}



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
nona<-which(is.na(ph[,1])==FALSE)
ko<-kmeans(x=ph[nona,1:2],center=2,iter.max=20,nstart=20)
if(sum(ko$cluster==1) > sum(ko$cluster==2)){ ## flip the assignments
	ko$cluster<-3-ko$cluster
}
plot(ph$V1[nona],ph$V2[nona],pch=19,col=c("brown","green")[ko$cluster])

grp<-rep(NA,dim(ph)[1])
grp[nona]<-ko$cluster
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

## compute phenotypic shifts

## maximize LD ##

## mean genotypes for each group
mnGg1<-apply(gcent[which(grp==1),],2,mean)
mnGg2<-apply(gcent[which(grp==2),],2,mean)

## adjustments to mean genotypes, twice af difference
adj1<-2*(pmax1-p1)
adj2<-2*(pmax2-p2)
mnGg1max<-mnGg1
mnGg2max<-mnGg2
mnGg1max[hipip]<-mnGg1max[hipip]+adj1
mnGg2max[hipip]<-mnGg2max[hipip]+adj2

deltaRGg1x<-(betaRG %*% mnGg1max)-(betaRG %*% mnGg1)
deltaRGg2x<-(betaRG %*% mnGg2max)-(betaRG %*% mnGg2)
deltaGBg1x<-(betaGB %*% mnGg1max)-(betaGB %*% mnGg1)
deltaGBg2x<-(betaGB %*% mnGg2max)-(betaGB %*% mnGg2)

rgMax<-pscoreRG
gbMax<-pscoreGB
rgMax[which(grp==1)]<-rgMax[which(grp==1)]+as.vector(deltaRGg1x)
rgMax[which(grp==2)]<-rgMax[which(grp==2)]+as.vector(deltaRGg2x)
gbMax[which(grp==1)]<-gbMax[which(grp==1)]+as.vector(deltaGBg1x)
gbMax[which(grp==2)]<-gbMax[which(grp==2)]+as.vector(deltaGBg2x)

## minimize LD ##


## adjustments to mean genotypes, twice af difference
adj1<-2*(phi[hipip]-p1)
adj2<-2*(phi[hipip]-p2)
mnGg1min<-mnGg1
mnGg2min<-mnGg2
mnGg1min[hipip]<-mnGg1min[hipip]+adj1
mnGg2min[hipip]<-mnGg2min[hipip]+adj2

deltaRGg1n<-(betaRG %*% mnGg1)-(betaRG %*% mnGg1min)
deltaRGg2n<-(betaRG %*% mnGg2)-(betaRG %*% mnGg2min)
deltaGBg1n<-(betaGB %*% mnGg1)-(betaGB %*% mnGg1min)
deltaGBg2n<-(betaGB %*% mnGg2)-(betaGB %*% mnGg2min)

rgMin<-pscoreRG
gbMin<-pscoreGB
rgMin[which(grp==1)]<-rgMin[which(grp==1)]-as.vector(deltaRGg1n)
rgMin[which(grp==2)]<-rgMin[which(grp==2)]-as.vector(deltaRGg2n)
gbMin[which(grp==1)]<-gbMin[which(grp==1)]-as.vector(deltaGBg1n)
gbMin[which(grp==2)]<-gbMin[which(grp==2)]-as.vector(deltaGBg2n)


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
		sRG<-(rgMin-pmns[1])/psds[1]
		sGB<-(gbMin-pmns[2])/psds[2]

		lpgMM<-mcb_MM[x,6] + mcb_MM[x,1] * sRG + mcb_MM[x,2] * sGB + mcb_MM[x,3] * sRG^2 + 
			mcb_MM[x,4] * sGB^2 + mcb_MM[x,5] * sRG * sGB
		rpgMM<-1/(1+exp(-1 * lpgMM))

		lpgAC<-mcb_AC[x,6] + mcb_AC[x,1] * sRG + mcb_AC[x,2] * sGB + mcb_AC[x,3] * sRG^2 + 
			mcb_AC[x,4] * sGB^2 + mcb_AC[x,5] * sRG * sGB
		rpgAC<-1/(1+exp(-1 * lpgAC))
		ran_mfMM[i,j]<-mean(rpgMM)
		ran_mfAC[i,j]<-mean(rpgAC)

		## antirandom
		sRG<-(rgMax-pmns[1])/psds[1]
		sGB<-(gbMax-pmns[2])/psds[2]
	
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

## prob > 1 .. looks good for 0.05 
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
ff<-paste("relativeFitnessUncCh8shift_red_",pip,".pdf",sep="")
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

ff<-paste("predictCh8_shift_red_",pip,".rdat",sep="")
save(list=ls(),file=ff)


