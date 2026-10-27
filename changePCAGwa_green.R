## PC clusters for experiment
## change in association
library(data.table)
library(scales)

## graphical parameters
cl<-1.4;ca<-1.1

## get het dat
#het<-as.matrix(fread("HetAndP_color.txt",header=TRUE))
## SNPs with het were already dropped in the GWA, drop those first
#hetx<-het[het[,1] != 1,]

## read gwa results for green RG
res_g_RG<-fread("output/o_green_RG.assoc.txt",header=TRUE)
hsubSNPs_green<-fread("hsubSNPs_green.txt",header=FALSE)
snp_green<-cbind(as.numeric(unlist(hsubSNPs_green[,1])),as.numeric(unlist(hsubSNPs_green[,2])))

## -log10 P from LRG
lp_green_RG<--log10(res_g_RG$p_lrt)
L<-length(lp_green_RG)
sig_green_RG<-which(lp_green_RG > 2) ## get rid of extreme hets and rare
sig_green_RG_ch8<-which(lp_green_RG > 2 & snp_green[,1]==4)
all_green_RG_ch8<-which(snp_green[,1]==4)

par(mar=c(5,5,1,1))
plot(snp_green[sig_green_RG_ch8,2],lp_green_RG[sig_green_RG_ch8],ylim=c(0,30),pch=19,xlab="Position (bp)",ylab="Green RG -Log10(P)",cex.lab=cl)
plot(snp_green[sig_green_RG_ch8,2],lp_green_RG[sig_green_RG_ch8],xlim=c(50000000,70000000),ylim=c(0,30),pch=19,xlab="Position (bp)",ylab="Green RG -Log10(P)",cex.lab=cl)

## read gwa results for green GB
res_g_GB<-fread("output/o_green_GB.assoc.txt",header=TRUE)

## -log10 P from LGB
lp_green_GB<--log10(res_g_GB$p_lrt)
L<-length(lp_green_GB)
sig_green_GB<-which(lp_green_GB > 2) ## get rid of extreme hets and rare
sig_green_GB_ch8<-which(lp_green_GB > 2 & snp_green[,1]==4)
all_green_GB_ch8<-which(snp_green[,1]==4)

par(mar=c(5,5,1,1))
plot(snp_green[sig_green_GB_ch8,2],lp_green_GB[sig_green_GB_ch8],ylim=c(0,30),pch=19,xlab="Position (bp)",ylab="Green GB -Log10(P)",cex.lab=cl)
plot(snp_green[sig_green_GB_ch8,2],lp_green_GB[sig_green_GB_ch8],xlim=c(50000000,70000000),ylim=c(0,30),pch=19,xlab="Position (bp)",ylab="Green GB -Log10(P)",cex.lab=cl)


### combined plots
plot(snp_green[sig_green_RG_ch8,2],lp_green_RG[sig_green_RG_ch8],xlim=c(50000000,70000000),ylim=c(0,30),pch=19,xlab="Position (bp)",ylab="Green RG -Log10(P)",cex.lab=cl,col="firebrick")

points(snp_green[sig_green_GB_ch8,2],lp_green_GB[sig_green_GB_ch8],pch=19,col="cadetblue")

## combine set in region
zGB<-which(lp_green_GB > 2 & snp_green[,1]==4 & snp_green[,2] >58000000 & snp_green[,2] <68000000)
zRG<-which(lp_green_RG > 2 & snp_green[,1]==4 & snp_green[,2] >58000000 & snp_green[,2] <68000000)
comblp<-c(lp_green_RG[zRG],lp_green_GB[zGB])
combSNPs<-c(snp_green[zRG,2],snp_green[zGB,2])

## RG GB indexes
cc<-rep(c(1,2),c(length(zRG),length(zGB)))
length(cc)
#[1] 35862

plot(sort(comblp))
qq<-quantile(comblp,probs=c(.5,.8,.9,.95,.99))
abline(h=qq)

sort_comblp<-comblp[order(combSNPs)]
sort_combSNPs<-combSNPs[order(combSNPs)]
sort_cc<-cc[order(combSNPs)]

plot(sort_combSNPs[(sort_comblp > 10)])


## r1
r1lb<-min(sort_combSNPs[(sort_comblp > 10)])#60060346
r1ub<-60682251

## r2
r2lb<-62703024
r2ub<-63185074

## r3
r3lb<-64589795
r3ub<-64862675

cols<-alpha(c("firebrick","royalblue3"),.5)

pdf("regions.pdf",width=8,height=4)
par(mar=c(5,5,1,1))
plot(sort_combSNPs,sort_comblp,type="n",xlab="Position (bp)",ylab="-log10(P)",cex.lab=1.3)
polygon(c(r1lb,r1ub,r1ub,r1lb),c(-5,-5,45,45),col=alpha("brown",.3),border=NA)
polygon(c(r2lb,r2ub,r2ub,r2lb),c(-5,-5,45,45),col=alpha("brown",.3),border=NA)
polygon(c(r3lb,r3ub,r3ub,r3lb),c(-5,-5,45,45),col=alpha("brown",.3),border=NA)
points(sort_combSNPs,sort_comblp,pch=19,col=alpha("gray30",.4))
mns<-c((r1lb+r1ub)/2,(r2lb+r2ub)/2,(r3lb+r3ub)/2)
text(mns,c(33,33,33),c("r1","r2","r3"),cex=1.3)
dev.off()

G<-fread("hsub_mlpntest_tchum_green.geno",header=FALSE,sep=",")
ph<-read.table("ph_tchum2019.txt",header=FALSE)
gmat<-as.matrix(G[,-c(1:3)])
ch<-snp_green[,1]
pos<-snp_green[,2]

## 4 is ch8

o1<-prcomp(t(gmat[r1,]),center=TRUE,scale=FALSE)
o2<-prcomp(t(gmat[r2,]),center=TRUE,scale=FALSE)
o3<-prcomp(t(gmat[r3,]),center=TRUE,scale=FALSE)


r1<-which(ch==4 & pos >= r1lb & pos <= r1ub)
r2<-which(ch==4 & pos >= r2lb & pos <= r2ub)
r3<-which(ch==4 & pos >= r3lb & pos <= r3ub)



## grab 1st 2 PCs from each region
pcs<-cbind(o1$x[,1:2],o2$x[,1:2],o3$x[,1:2])

## read treatment
trt<-read.table("sort_green_phfull_tchum2019.txt",header=TRUE,comment.char="%",sep=",")
## same order as ph
t_ac<-which(trt$Treatment=="AC")
t_mm<-which(trt$Treatment=="MM")

t_ac_surv<-which(trt$Treatment=="AC" & trt$Survival==1)
t_mm_surv<-which(trt$Treatment=="MM" & trt$Survival==1)

##
cors_ac<-abs(cor(pcs[t_ac,]))
cors_mm<-abs(cor(pcs[t_mm,]))
cors_sac<-abs(cor(pcs[t_ac_surv,]))
cors_smm<-abs(cor(pcs[t_mm_surv,]))

obs_ac<-mean(cors_ac[upper.tri(cors_ac)])
obs_mm<-mean(cors_mm[upper.tri(cors_mm)])
obs_acs<-mean(cors_sac[upper.tri(cors_sac)])
obs_mms<-mean(cors_smm[upper.tri(cors_smm)])

## relative to random mortality
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
#[1] 0.008
obs_acs
#[1] 0.4602588
obs_ac
#[1]  0.371145
obs_acs-obs_ac
#[1] 0.08911382
## null change
summary(null_acs-obs_ac)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.042474 -0.008203  0.005998  0.008821  0.023257  0.148506 
## null end
summary(null_acs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3287  0.3629  0.3771  0.3800  0.3944  0.5197
(obs_acs-obs_ac)/mean(null_acs-obs_ac)
#[1] 10.10243 ## 10X bigger increase than mean null
mean(null_acs-obs_ac)
#[1] 0.008821029
obs_acs-obs_ac
#[1] 0.08911382

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
#[1] 0.449
obs_mms
#[1] 0.3709778
obs_mm
#[1] 0.3510696
obs_mms-obs_mm
#[1] 0.01990821
## null change
summary(null_mms-obs_mm)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#-0.0402548  0.0008169  0.0161588  0.0203591  0.0346106  0.1379995 
## null end
summary(null_mms)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3108  0.3519  0.3672  0.3714  0.3857  0.4891 
(obs_mms-obs_mm)/mean(null_mms-obs_mm)
#[1] 0.9778532 ## bascially 1X expectations
mean(null_mms-obs_mm)
#[1] 0.0203591
obs_mms-obs_mm
#[1] 0.01990821

## increased, but what is null??, allow for
## linear selection

## AC
null_acs_lrg<-rep(NA,1000)
null_acs_lgb<-rep(NA,1000)

o_ac_rg<-glm(ph[t_ac,3] ~ ph[t_ac,1],family="binomial")
o_ac_gb<-glm(ph[t_ac,3] ~ ph[t_ac,2],family="binomial")

p_ac_rg<-ph[t_ac,3]
p_ac_rg[is.na(p_ac_rg)==FALSE]<-o_ac_rg$fitted.values
p_ac_rg[is.na(p_ac_rg)==TRUE]<-mean(o_ac_rg$fitted.values,na.rm=TRUE)

p_ac_gb<-ph[t_ac,3]
p_ac_gb[is.na(p_ac_gb)==FALSE]<-o_ac_gb$fitted.values
p_ac_gb[is.na(p_ac_gb)==TRUE]<-mean(o_ac_gb$fitted.values,na.rm=TRUE)


for(i in 1:1000){
        sri<-sample(t_ac,size=sum(ph[t_ac,3]==1,na.rm=TRUE),replace=FALSE,prob=p_ac_rg)
        sran<-rep(0,length(t_ac)+length(t_mm))
        sran[sri]<-1
        r_ac_surv<-which(sran==1)
        rcors_acs<-abs(cor(pcs[r_ac_surv,]))
        rcorsbar_acs<-mean(rcors_acs[upper.tri(rcors_acs)])
        null_acs_lrg[i]<-rcorsbar_acs

        sri<-sample(t_ac,size=sum(ph[t_ac,3]==1,na.rm=TRUE),replace=FALSE,prob=p_ac_gb)
        sran<-rep(0,length(t_ac)+length(t_mm))
        sran[sri]<-1
        r_ac_surv<-which(sran==1)
        rcors_acs<-abs(cor(pcs[r_ac_surv,]))
        rcorsbar_acs<-mean(rcors_acs[upper.tri(rcors_acs)])
        null_acs_lgb[i]<-rcorsbar_acs
}

## RG
mean(null_acs_lrg >= obs_acs)
#[1] 0.008
obs_acs
#[1] 0.4602588
obs_ac
#[1] 0.371145
obs_acs-obs_ac
#[1] 0.08911382
## null change
summary(null_acs_lrg-obs_ac)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.041445 -0.005538  0.008124  0.011547  0.025933  0.115231 
## null end
summary(null_acs_lrg)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3297  0.3656  0.3793  0.3827  0.3971  0.4864 
(obs_acs-obs_ac)/mean(null_acs_lrg-obs_ac)
#[1] 7.717554 ## almost 8X bigger increase than mean null
mean(null_acs_lrg-obs_ac)
#[1] 0.0115469
obs_acs-obs_ac
#[1] 0.08911382

## GB
mean(null_acs_lgb >= obs_acs)
#[1] 0.009
obs_acs
#[1] 0.4602588
obs_ac
#[1] 0.371145
obs_acs-obs_ac
#[1] 0.08911382
## null change
summary(null_acs_lgb-obs_ac)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.047875 -0.004816  0.007581  0.011911  0.024899  0.126256 
## null end
summary(null_acs_lgb)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3233  0.3663  0.3787  0.3831  0.3960  0.4974 
(obs_acs-obs_ac)/mean(null_acs_lgb-obs_ac)
#[1] 7.481587 ## >7X bigger increase than mean null
mean(null_acs_lgb-obs_ac)
#[1] 0.01191108
obs_acs-obs_ac
#[1] 0.08911382

### MM
null_mms_lrg<-rep(NA,1000)
null_mms_lgb<-rep(NA,1000)

o_mm_rg<-glm(ph[t_mm,3] ~ ph[t_mm,1],family="binomial")
o_mm_gb<-glm(ph[t_mm,3] ~ ph[t_mm,2],family="binomial")

p_mm_rg<-ph[t_mm,3]
p_mm_rg[is.na(p_mm_rg)==FALSE]<-o_mm_rg$fitted.values
p_mm_rg[is.na(p_mm_rg)==TRUE]<-mean(o_mm_rg$fitted.values,na.rm=TRUE)

p_mm_gb<-ph[t_mm,3]
p_mm_gb[is.na(p_mm_gb)==FALSE]<-o_mm_gb$fitted.values
p_mm_gb[is.na(p_mm_gb)==TRUE]<-mean(o_mm_gb$fitted.values,na.rm=TRUE)


for(i in 1:1000){
        sri<-sample(t_mm,size=sum(ph[t_mm,3]==1,na.rm=TRUE),replace=FALSE,prob=p_mm_rg)
        sran<-rep(0,length(t_ac)+length(t_mm))
        sran[sri]<-1
        r_mm_surv<-which(sran==1)
        rcors_mms<-abs(cor(pcs[r_mm_surv,]))
        rcorsbar_mms<-mean(rcors_mms[upper.tri(rcors_mms)])
        null_mms_lrg[i]<-rcorsbar_mms

        sri<-sample(t_mm,size=sum(ph[t_mm,3]==1,na.rm=TRUE),replace=FALSE,prob=p_mm_gb)
        sran<-rep(0,length(t_ac)+length(t_mm))
        sran[sri]<-1
        r_mm_surv<-which(sran==1)
        rcors_mms<-abs(cor(pcs[r_mm_surv,]))
        rcorsbar_mms<-mean(rcors_mms[upper.tri(rcors_mms)])
        null_mms_lgb[i]<-rcorsbar_mms
}

## RG
mean(null_mms_lrg >= obs_mms)
#[1] 0.510
obs_mms 
#[1] 0.3709778
obs_mm 
#[1] 0.3510696
obs_mms-obs_mm 
#[1] 0.01990821
## null change 
summary(null_mms_lrg-obs_mm) 
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.047970  0.003862  0.020813  0.023604  0.039531  0.125431
## null end 
summary(null_mms_lrg) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3031  0.3549  0.3719  0.3747  0.3906  0.4765 
(obs_mms-obs_mm)/mean(null_mms_lrg-obs_mm)
#[1] 0.8434296 ## bascially 1X expectations
mean(null_mms_lrg-obs_mm) 
#[1] 0.02360388
obs_mms-obs_mm 
#[1] 0.01990821

## GB
mean(null_mms_lgb >= obs_mms)
#[1] 0.472
obs_mms 
#[1] 0.3709778
obs_mm 
#[1] 0.3510696
obs_mms-obs_mm 
#[1] 0.01990821
## null change 
summary(null_mms_lgb-obs_mm) 
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.034940  0.002953  0.017864  0.021530  0.035609  0.116433
## null end 
summary(null_mms_lgb) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3161  0.3540  0.3689  0.3726  0.3867  0.4675
(obs_mms-obs_mm)/mean(null_mms_lgb-obs_mm)
#[1] 0.9246647 ## bascially 1X expectations
mean(null_mms_lgb-obs_mm) 
#[1] 0.0215302
obs_mms-obs_mm 
#[1] 0.01990821

save(list=ls(),file="LD_green.rdat")

#### NOT UPDATED #####

########### within regions vs among regions ####################
cors_ac<-abs(cor(pcs[t_ac,]))
cors_mm<-abs(cor(pcs[t_mm,]))
cors_sac<-abs(cor(pcs[t_ac_surv,]))
cors_smm<-abs(cor(pcs[t_mm_surv,]))

ttt<-c("A","A","B","B","C","C")
ttm<-matrix(NA,nrow=6,ncol=6)
for(i in 1:6){for(j in 1:6){if(i!=j){ttm[i,j]<-as.numeric(ttt[i]==ttt[j])}}}
ttm ## 1 = within region, 0 = between region
#     [,1] [,2] [,3] [,4] [,5] [,6]
#[1,]   NA    1    0    0    0    0
#[2,]    1   NA    0    0    0    0
#[3,]    0    0   NA    1    0    0
#[4,]    0    0    1   NA    0    0
#[5,]    0    0    0    0   NA    1
#[6,]    0    0    0    0    1   NA



obs_acW<-mean(cors_ac[upper.tri(cors_ac)][ttm[upper.tri(ttm)]==1])
obs_mmW<-mean(cors_mm[upper.tri(cors_mm)][ttm[upper.tri(ttm)]==1])
obs_acsW<-mean(cors_sac[upper.tri(cors_sac)][ttm[upper.tri(ttm)]==1])
obs_mmsW<-mean(cors_smm[upper.tri(cors_smm)][ttm[upper.tri(ttm)]==1])


### increased, but what is null??
for(i in 1:1000){
        sri<-sample(t_ac,size=sum(ph[t_ac,3]==1,na.rm=TRUE),replace=FALSE,prob=p_ac_rg)
        sran<-rep(0,length(t_ac)+length(t_ac))
        sran[sri]<-1
        r_ac_surv<-which(sran==1)
	rcors_acs<-abs(cor(pcs[r_ac_surv,]))
        rcorsbar_acs<-mean(rcors_acs[upper.tri(rcors_acs)][ttm[upper.tri(ttm)]==1])
        null_acs_lrg[i]<-rcorsbar_acs

	sri<-sample(t_ac,size=sum(ph[t_ac,3]==1,na.rm=TRUE),replace=FALSE,prob=p_ac_gb)
        sran<-rep(0,length(t_ac)+length(t_ac))
        sran[sri]<-1
        r_ac_surv<-which(sran==1)
        rcors_acs<-abs(cor(pcs[r_ac_surv,]))
        rcorsbar_acs<-mean(rcors_acs[upper.tri(rcors_acs)][ttm[upper.tri(ttm)]==1])
        null_acs_lgb[i]<-rcorsbar_acs
}

## RG
mean(null_acs_lrg >= obs_acsW)
#[1] 0.001
obs_acsW
#[1] 0.2590319
obs_acW
#[1] 0.04927402
obs_acsW-obs_acW
#[1] 0.2097579
summary(null_acs_lrg)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.005544 0.057572 0.082646 0.089741 0.116349 0.274767 
(obs_acsW-obs_acW)/mean(null_acs_lrg-obs_acW)
#[1] 5.332236 ## 5X bigger increase than mean null
mean(null_acs_lrg-obs_acW)
#[1] 0.0393377
obs_acsW-obs_acW
#[1]  0.2097579

## GB
mean(null_acs_lgb >= obs_acsW)
#[1] 0.002
obs_acsW
#[1] 0.2590319
obs_acW
#[1] 0.04927402
obs_acsW-obs_acW
#[1] 0.2097579
summary(null_acs_lgb)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.006562 0.056270 0.083955 0.091612 0.118613 0.299865  
(obs_acsW-obs_acW)/mean(null_acs_lgb-obs_acW)
#[1] 4.954414 ## 5X bigger increase than mean null
mean(null_acs_lgb-obs_acW)
#[1] 0.04233758
obs_acsW-obs_acW
#[1]  0.2097579


## MM
for(i in 1:1000){
        sri<-sample(t_mm,size=sum(ph[t_mm,3]==1,na.rm=TRUE),replace=FALSE,prob=p_mm_rg)
        sran<-rep(0,length(t_ac)+length(t_mm))
        sran[sri]<-1
        r_mm_surv<-which(sran==1)
	rcors_mms<-abs(cor(pcs[r_mm_surv,]))
        rcorsbar_mms<-mean(rcors_mms[upper.tri(rcors_mms)][ttm[upper.tri(ttm)]==1])
        null_mms_lrg[i]<-rcorsbar_mms
        
	sri<-sample(t_mm,size=sum(ph[t_mm,3]==1,na.rm=TRUE),replace=FALSE,prob=p_mm_gb)
        sran<-rep(0,length(t_ac)+length(t_mm))
        sran[sri]<-1
        r_mm_surv<-which(sran==1)
	rcors_mms<-abs(cor(pcs[r_mm_surv,]))
        rcorsbar_mms<-mean(rcors_mms[upper.tri(rcors_mms)][ttm[upper.tri(ttm)]==1])
        null_mms_lgb[i]<-rcorsbar_mms
}

## RG
mean(null_mms_lrg >= obs_mmsW)
#[1] 0.632
obs_mmsW 
#[1] 0.07916729
obs_mmW 
#[1] 0.04215004
obs_mmsW-obs_mmW 
#[1] 0.03701724
## null end 
summary(null_mms_lrg) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.006809 0.066128 0.096106 0.101812 0.130770 0.309235   
(obs_mmsW-obs_mmW)/mean(null_mms_lrg-obs_mmW)
#[1] 0.6204501 ## <1X expectations
mean(null_mms_lrg-obs_mmW) 
#[1] 0.05966192
obs_mmsW-obs_mmW 
#[1] 0.03701724

## GB
mean(null_mms_lgb >= obs_mmsW)
#[1] 0.645
obs_mmsW 
#[1] 0.07916729
obs_mmW 
#[1] 0.04215004
obs_mmsW-obs_mmW 
#[1] 0.03701724
## null end 
summary(null_mms_lgb) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00559 0.06578 0.09693 0.10348 0.12975 0.31666   
(obs_mmsW-obs_mmW)/mean(null_mms_lgb-obs_mmW)
#[1] 0.6035303 ## <1X expectations
mean(null_mms_lgb-obs_mmW) 
#[1] 0.06133452
obs_mmsW-obs_mmW 
#[1] 0.03701724


obs_acB<-mean(cors_ac[upper.tri(cors_ac)][ttm[upper.tri(ttm)]==0])
obs_mmB<-mean(cors_mm[upper.tri(cors_mm)][ttm[upper.tri(ttm)]==0])
obs_acsB<-mean(cors_sac[upper.tri(cors_sac)][ttm[upper.tri(ttm)]==0])
obs_mmsB<-mean(cors_smm[upper.tri(cors_smm)][ttm[upper.tri(ttm)]==0])


### increased, but what is null??

##AC
for(i in 1:1000){
        sri<-sample(t_ac,size=sum(ph[t_ac,3]==1,na.rm=TRUE),replace=FALSE,prob=p_ac_rg)
        sran<-rep(0,length(t_ac)+length(t_ac))
        sran[sri]<-1
        r_ac_surv<-which(sran==1)       
       	rcors_acs<-abs(cor(pcs[r_ac_surv,]))
        rcorsbar_acs<-mean(rcors_acs[upper.tri(rcors_acs)][ttm[upper.tri(ttm)]==0])
        null_acs_lrg[i]<-rcorsbar_acs
        
	sri<-sample(t_ac,size=sum(ph[t_ac,3]==1,na.rm=TRUE),replace=FALSE,prob=p_ac_gb)
        sran<-rep(0,length(t_ac)+length(t_ac))
        sran[sri]<-1
        r_ac_surv<-which(sran==1)       
       	rcors_acs<-abs(cor(pcs[r_ac_surv,]))
        rcorsbar_acs<-mean(rcors_acs[upper.tri(rcors_acs)][ttm[upper.tri(ttm)]==0])
        null_acs_lgb[i]<-rcorsbar_acs
}

## RG
mean(null_acs_lrg >= obs_acsB)
#[1] 0.001
obs_acsB
#[1] 0.5422799
obs_acB
#[1] 0.4508417
obs_acsB-obs_acB
#[1]  0.09143818
summary(null_acs_lrg)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.3688  0.4428  0.4625  0.4618  0.4821  0.5447 
(obs_acsB-obs_acB)/mean(null_acs_lrg-obs_acB)
#[1] 8.369867 ## 8X bigger increase than mean null
mean(null_acs_lrg-obs_acB)
#[1]  0.01092469
obs_acsB-obs_acB
#[1] 0.09143818

## GB
mean(null_acs_lgb >= obs_acsB)
#[1] 0.002
obs_acsB
#[1] 0.5422799
obs_acB
#[1] 0.4508417
obs_acsB-obs_acB
#[1]  0.09143818
summary(null_acs_lgb)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.3900  0.4524  0.4695  0.4703  0.4881  0.5494 
(obs_acsB-obs_acB)/mean(null_acs_lgb-obs_acB)
#[1] 4.709238 ## 5X bigger increase than mean null
mean(null_acs_lgb-obs_acB)
#[1] 0.01941677
obs_acsB-obs_acB
#[1] 0.09143818


#MM
for(i in 1:1000){
        sri<-sample(t_mm,size=sum(ph[t_mm,3]==1,na.rm=TRUE),replace=FALSE,prob=p_mm_rg)
        sran<-rep(0,length(t_ac)+length(t_mm))
        sran[sri]<-1
        r_mm_surv<-which(sran==1)        
	rcors_mms<-abs(cor(pcs[r_mm_surv,]))
        rcorsbar_mms<-mean(rcors_mms[upper.tri(rcors_mms)][ttm[upper.tri(ttm)]==0])
        null_mms_lrg[i]<-rcorsbar_mms
        
	sri<-sample(t_mm,size=sum(ph[t_mm,3]==1,na.rm=TRUE),replace=FALSE,prob=p_mm_gb)
        sran<-rep(0,length(t_ac)+length(t_mm))
        sran[sri]<-1
        r_mm_surv<-which(sran==1)        
	rcors_mms<-abs(cor(pcs[r_mm_surv,]))
        rcorsbar_mms<-mean(rcors_mms[upper.tri(rcors_mms)][ttm[upper.tri(ttm)]==0])
        null_mms_lgb[i]<-rcorsbar_mms
}

## RG
mean(null_mms_lrg >= obs_mmsB)
#[1] 0.452
obs_mmsB 
#[1] 0.4549934
obs_mmB 
#[1] 0.4299155
obs_mmsB-obs_mmB 
#[1] 0.02507789
## null end 
summary(null_mms_lrg) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.3527  0.4302  0.4512  0.4519  0.4739  0.5728   
(obs_mmsB-obs_mmB)/mean(null_mms_lrg-obs_mmB)
#[1] 1.141933 ## 1X expectations
mean(null_mms_lrg-obs_mmB) 
#[1] 0.02196091 
obs_mmsB-obs_mmB 
#[1] 0.02507789

## GB
mean(null_mms_lgb >= obs_mmsB)
#[1] 0.43
obs_mmsB 
#[1] 0.4549934
obs_mmB 
#[1] 0.4299155
obs_mmsB-obs_mmB 
#[1] 0.02507789
## null end 
summary(null_mms_lgb) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.3399  0.4241  0.4479  0.4483  0.4738  0.5398    
(obs_mmsB-obs_mmB)/mean(null_mms_lgb-obs_mmB)
#[1] 1.364702 ## 1X expectations
mean(null_mms_lgb-obs_mmB) 
#[1] 0.0183761 
obs_mmsB-obs_mmB 
#[1] 0.02507789

save(list=ls(),file="ChangePC_LS.rdat")

