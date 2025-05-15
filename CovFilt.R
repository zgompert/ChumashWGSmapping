## compute depth per individual and SNP to identify SNPs and individuals to drop
## the idea is to get rid of low coverage individuals
## and SNPs with either very high coverage (3SD > mean) or very high variance in coverage
## across individuals

library(data.table)
d<-as.matrix(fread("depth.txt",header=FALSE))

## mean and SD by SNP
mnc<-apply(d,1,mean)
sdc<-apply(d,1,sd)

## CV 
cvc<-sdc/mnc
mean(cvc);quantile(cvc,probs=c(.5,.9,.95,.99,.999,1))
#[1] 0.922443
#      50%       90%       95%       99%     99.9%      100%
#0.8930029 1.0995234 1.1805665 1.3386413 1.5796200 4.6469001

mean(mnc)+3*sd(mnc)
#[1] 20.77883

mean(mnc);quantile(mnc,probs=c(.5,.9,.95,.99,.999,1))
#[1] 8.162421
#       50%        90%        95%        99%      99.9%       100%
#  8.038235  11.272549  11.653922  13.056863  29.240529 772.807843


## for SNPs, keep if CV < 1.5 and mean < 21
keepSNPs<-as.numeric(mnc < 21 & cvc < 1.5)
mean(keepSNPs)
#[1] 0.9959606
sum(keepSNPs)
#[1] 35061459

## for individuals

mni<-apply(d,2,mean)
plot(sort(mni))
summary(mni)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#  0.00376   6.34337   7.27581   8.16242   8.48578 113.62793

quantile(mni,probs=c(.025,.99))
#     2.5%       99%
# 3.988486 38.582537

keepInds<-as.numeric(mni > 4 & mni < 40)
mean(keepInds)
#[1] 0.9647059
sum(keepInds)
#[1] 984

write.table(file="KeepInds.txt",keepInds,row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(file="KeepSNPs.txt",keepSNPs,row.names=FALSE,col.names=FALSE,quote=FALSE)
