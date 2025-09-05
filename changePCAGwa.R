## PC clusters for experiment
## change in association
library(data.table)
library(scales)

## graphical parameters
cl<-1.4;ca<-1.1

## get het dat
het<-as.matrix(fread("HetAndP_color.txt",header=TRUE))
## SNPs with het were already dropped in the GWA, drop those first
hetx<-het[het[,1] != 1,]

## read gwa results for RG
resRG<-fread("output/o_RG.assoc.txt",header=TRUE)

## chrom and pos for RG
snpsRG<-as.matrix(fread("SNPsRG.txt",header=FALSE))

## -log10 P from LRG
lpRG<--log10(resRG$p_lrt)
L<-length(lpRG)
##  Chr8 only
sigRG8<-which(lpRG > 2 & snpsRG[,1]==8)
sigRG8lax<-which(lpRG > 1 & snpsRG[,1]==8)

par(mar=c(5,5,1,1))
plot(snpsRG[sigRG8,2],lpRG[sigRG8],ylim=c(0,27),pch=19,xlab="Position (bp)",ylab="RG -Log10(P)",cex.lab=cl)

plot(snpsRG[sigRG8lax,2],lpRG[sigRG8lax],ylim=c(0,27),pch=19,xlim=c(50000000,70000000),xlab="Position (bp)",ylab="RG -Log10(P)",cex.lab=cl)

zRG<-which(hetx[,3] > -.2 & hetx[,3] < .2 & hetx[,2] > .025 & hetx[,2] < .975 & snpsRG[,1]==8 & snpsRG[,2]>50000000 & snpsRG[,2]<70000000)

plot(snpsRG[zRG,2],lpRG[zRG],pch=19,xlab="Position (bp)",ylab="RG -Log10(P)",cex.lab=cl)

zRG<-which(hetx[,3] > -.2 & hetx[,3] < .2 & hetx[,2] > .025 & hetx[,2] < .975 & snpsRG[,1]==8 & snpsRG[,2]>55000000 & snpsRG[,2]<63000000)

plot(snpsRG[zRG,2],lpRG[zRG],pch=19,xlab="Position (bp)",ylab="RG -Log10(P)",cex.lab=cl)

## read gwa results for GB
resGB<-fread("output/o_GB.assoc.txt",header=TRUE)

## chrom and pos for RG
snpsGB<-as.matrix(fread("SNPsGB.txt",header=FALSE))

## -log10 P from LRG
lpGB<--log10(resGB$p_lrt)
L<-length(lpGB)

## Chr8 only
sigGB8<-which(lpGB > 2 & snpsGB[,1]==8)

par(mar=c(5,5,1,1))
plot(snpsGB[sigGB8,2],lpGB[sigGB8],ylim=c(0,33.5),pch=19,xlab="Position (bp)",ylab="GB -Log10(P)",cex.lab=cl)

zGB<-which(hetx[,3] > -.2 & hetx[,3] < .2 & hetx[,2] > .025 & hetx[,2] < .975 & snpsGB[,1]==8 & snpsGB[,2]>55000000 & snpsGB[,2]<63000000)

points(snpsGB[zGB,2],lpGB[zGB],pch=19,col="cadetblue")



### combined plots
plot(snpsRG[zRG,2],lpRG[zRG],pch=19,xlab="Position (bp)",ylab="RG -Log10(P)",cex.lab=cl,col="darkgray")
points(snpsGB[zGB,2],lpGB[zGB],pch=19,col="orange3")

## combine set in region
comblp<-c(lpRG[zRG],lpGB[zGB])
combSNPs<-c(snpsRG[zRG,2],snpsGB[zGB,2])

plot(sort(comblp))
qq<-quantile(comblp,probs=c(.5,.8,.9,.95,.99))
#       50%        80%        90%        95%        99%
# 0.4035478  1.1106558  2.3519926  6.1441750 16.9720125
abline(h=qq)

sort_comblp<-comblp[order(combSNPs)]
sort_combSNPs<-combSNPs[order(combSNPs)]

plot(sort_combSNPs[(sort_comblp > 10)])


## r1
r1lb<-min(sort_combSNPs[(sort_comblp > 10)])
r1ub<-57380207

## r2
r2lb<-58460759
r2ub<-59734122
## r3
r3lb<-62179278
r3ub<-max(sort_combSNPs[(sort_comblp > 10)])

plot(sort_combSNPs,sort_comblp,pch=19,col="gray")
abline(v=c(r1lb,r1ub,r2lb,r2ub,r3lb,r3ub),col=rep(c("orange3","red","violet"),each=2))

G<-fread("hsub_mlpntest_EXP_tchum.geno",header=FALSE,sep=",")
ph<-read.table("ph_tchum2019.txt",header=FALSE)
gmat<-as.matrix(G[,-c(1:3)])

strSNP<-unlist(strsplit(unlist(G[,1]),split=":"))
ch<-as.numeric(strSNP[seq(1,11995976,2)])
pos<-as.numeric(strSNP[seq(2,11995976,2)])

r1<-which(ch==8 & pos >= r1lb & pos <= r1ub)
r2<-which(ch==8 & pos >= r2lb & pos <= r2ub)
r3<-which(ch==8 & pos >= r3lb & pos <= r3ub)


o1<-prcomp(t(gmat[r1,]),center=TRUE,scale=FALSE)
o2<-prcomp(t(gmat[r2,]),center=TRUE,scale=FALSE)
o3<-prcomp(t(gmat[r3,]),center=TRUE,scale=FALSE)


## grab 1st 2 PCs from each region
pcs<-cbind(o1$x[,1:2],o2$x[,1:2],o3$x[,1:2])

## read treatment
trt<-read.table("sort_phfull_tchum2019.txt",header=TRUE,comment.char="%")
## same order as ph
t_ac<-which(trt$Treatment=="AC")
t_mm<-which(trt$Treatment=="MM")

t_ac_surv<-which(trt$Treatment=="AC" & ph[,3]==1)
t_mm_surv<-which(trt$Treatment=="MM" & ph[,3]==1)

##
cors_ac<-abs(cor(pcs[t_ac,]))
cors_mm<-abs(cor(pcs[t_mm,]))
cors_sac<-abs(cor(pcs[t_ac_surv,]))
cors_smm<-abs(cor(pcs[t_mm_surv,]))

obs_ac<-mean(cors_ac[upper.tri(cors_ac)])
obs_mm<-mean(cors_mm[upper.tri(cors_mm)])
obs_acs<-mean(cors_sac[upper.tri(cors_sac)])
obs_mms<-mean(cors_smm[upper.tri(cors_smm)])


## increased, but what is null??
sran<-ph[,3]
null_acs<-rep(NA,1000)
for(i in 1:1000){
        sran[t_ac]<-sample(sran[t_ac],length(t_ac),replace=FALSE)
        r_ac_surv<-which(trt$Treatment=="AC" & sran==1)
        rcors_acs<-abs(cor(pcs[r_ac_surv,]))
        rcorsbar_acs<-mean(rcors_acs[upper.tri(rcors_acs)])

        null_acs[i]<-rcorsbar_acs
}
mean(null_acs >= obs_acs)
#[1] 0.002
obs_acs
#[1] 0.4856303
obs_ac
#[1] 0.3705282
obs_acs-obs_ac
#[1] 0.1151021
## null change
summary(null_acs-obs_ac)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.066379 -0.001296  0.016038  0.017976  0.036287  0.123883 
## null end
summary(null_acs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3041  0.3692  0.3866  0.3885  0.4068  0.4944
(obs_acs-obs_ac)/mean(null_acs-obs_ac)
#[1] 6.40315 ## 6X bigger increase than mean null
mean(null_acs-obs_ac)
#[1] 0.01797586
obs_acs-obs_ac
#[1] 0.1151021


sran<-ph[,3]
null_mms<-rep(NA,1000)
for(i in 1:1000){
        sran[t_mm]<-sample(sran[t_mm],length(t_mm),replace=FALSE)
        r_mm_surv<-which(trt$Treatment=="MM" & sran==1)
        rcors_mms<-abs(cor(pcs[r_mm_surv,]))
        rcorsbar_mms<-mean(rcors_mms[upper.tri(rcors_mms)])

        null_mms[i]<-rcorsbar_mms
}
mean(null_mms >= obs_mms)
#[1] 0.446
obs_mms 
#[1] 0.3798282
obs_mm 
#[1] 0.3523624
obs_mms-obs_mm 
#[1] 0.02746576
## null change 
summary(null_mms-obs_mm) 
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.057630  0.003708  0.024062  0.025508  0.045558  0.129579 
## null end 
summary(null_mms) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2947  0.3561  0.3764  0.3779  0.3979  0.4819 
(obs_mms-obs_mm)/mean(null_mms-obs_mm)
#[1] 1.076731 ## bascially 1X expectations
mean(null_mms-obs_mm) 
#[1] 0.02550848
obs_mms-obs_mm 
#[1] 0.02746576
