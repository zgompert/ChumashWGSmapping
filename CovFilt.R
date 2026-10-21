## compute depth per individual and SNP to identify SNPs and individuals to drop
## the idea is to get rid of low coverage individuals
## and SNPs with either very high coverage (3SD > mean) or very high variance in coverage
## across individuals

library(data.table)
gf<-list.files(pattern="depth_tchum_green")
rf<-list.files(pattern="depth_tchum_red")

##### G
dG<-vector("list",11)
for(i in 1:11){
	dG[[i]]<-as.matrix(fread(gf[i],header=FALSE))
}

dGm<-do.call(rbind,dG)
## mean and SD by SNP, G and R
mncG<-apply(dGm,1,mean)
sdcG<-apply(dGm,1,sd)
cvcG<-sdcG/mncG
mean(cvcG);quantile(cvcG,probs=c(.5,.9,.95,.99,.999,1))
#[1] 0.922443
#      50%       90%       95%       99%     99.9%      100%
#0.8930029 1.0995234 1.1805665 1.3386413 1.5796200 4.6469001

mean(mncG)+3*sd(mncG)
#[1] 17.91721

mean(mncG);quantile(mncG,probs=c(.5,.9,.95,.99,.999,1))
#[1] 9.404392
#       50%        90%        95%        99%      99.9%       100%
#  8.038235  11.272549  11.653922  13.056863  29.240529 772.807843


## for SNPs, keep if CV < 1.5 and mean < 18
keepSNPs<-as.numeric(mncG < 18 & cvcG < 1.5)
mean(keepSNPs)
#[1] 0.9967751
sum(keepSNPs)
#[1]  53378633
write.table(file="KeepSNPs_G.txt",keepSNPs,row.names=FALSE,col.names=FALSE,quote=FALSE)

## for individuals

mni<-apply(dGm,2,mean)
plot(sort(mni))
summary(mni)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.004692  7.824062  9.235907  9.404392 10.685462 20.585886 


quantile(mni,probs=c(.025,.99))
#     2.5%       99%
# 3.948488 17.768417 

keepInds<-as.numeric(mni > 4)
mean(keepInds)
#[1] 0.9736264
sum(keepInds)
#[1] 443

write.table(file="KeepInds_G.txt",keepInds,row.names=FALSE,col.names=FALSE,quote=FALSE)


####
dR<-vector("list",11)
for(i in 1:11){
	dR[[i]]<-as.matrix(fread(rf[i],header=FALSE))
}

dRm<-do.call(rbind,dR)

mncR<-apply(dRm,1,mean)
sdcR<-apply(dRm,1,sd)

cvcR<-sdcR/mncR
mean(cvcR);quantile(cvcR,probs=c(.5,.9,.95,.99,.999,1))
#[1] 0.4808483
#      50%       90%       95%       99%     99.9%      100% 
#0.4607142 0.5624106 0.5976420 0.7144675 0.9435675 4.5407428 
mean(mncR)+3*sd(mncR)
#[1] 20.62528
mean(mncR);quantile(mncR,probs=c(.5,.9,.95,.99,.999,1))
#[1] 9.512111
#       50%        90%        95%        99%      99.9%       100% 
#  9.734066  10.597802  10.848352  14.180220  31.114286 992.984615 


## for SNPs, keep if CV < 1.5 and mean < 21
keepSNPs<-as.numeric(mncR < 21 & cvcR < 1.5)
mean(keepSNPs)
#[1] 0.9971557
sum(keepSNPs)
#[1] 53128874
write.table(file="KeepSNPs_R.txt",keepSNPs,row.names=FALSE,col.names=FALSE,quote=FALSE)

## for individuals

mni<-apply(dRm,2,mean)
plot(sort(mni))
summary(mni)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.004726  7.932072  9.347886  9.512111 10.787687 20.759556 




quantile(mni,probs=c(.025,.99))
#     2.5%       99%
#  4.046159 17.865337 

keepInds<-as.numeric(mni > 4.05) 
mean(keepInds)
#[1] 0.9736264
sum(keepInds)
#[1] 443

write.table(file="KeepInds_R.txt",keepInds,row.names=FALSE,col.names=FALSE,quote=FALSE)
