## summarize BSLMM GWA, green genome
library(data.table)
library(scales)

## graphical parameters
cl<-1.4;ca<-1.1

## read gwa results for RG
resRG<-fread("output/pips_o_poly_bslmm_GR8_ph1_ch0.param.txt",header=FALSE)


pips<-as.numeric(unlist(resRG[,5]))
L<-length(pips)
sigRG<-which(pips > 0.00001) 

bnds<-c(which(resRG[-L,1] != resRG[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=resRG[,1],mean)


## regions from single res
reg_green<-read.csv("ChumashColorRegions-green_genome.csv")
r1lb<-reg_green[1,6]
r1ub<-reg_green[1,7]
r2lb<-reg_green[1,8]
r2ub<-reg_green[1,9]
r3lb<-reg_green[1,10]
r3ub<-reg_green[1,11]

sum(pips[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r1ub])
#[1] 0.863769 
sum(pips[resRG[,1]==4 & resRG[,2] > r2lb & resRG[,2] < r2ub])
#[1] 1.050717
sum(pips[resRG[,1]==4 & resRG[,2] > r3lb & resRG[,2] < r3ub])
#[1] 0.001006667
sum(pips[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r3ub])
#[1]  2.0138
pipsRG<-pips

nr<-dim(reg_green)[1]
pipRG_region<-rep(NA,nr)
for(k in 1:nr){
	sp<-rep(NA,reg_green$num_windows[k])
	for(a in 1:reg_green$num_windows[k]){
		b<-(a-1)*2
		lb<-reg_green[k,6+b]
		ub<-reg_green[k,7+b]
		sp[a]<-sum(pips[resRG[,1]==4 & resRG[,2] > lb & resRG[,2] < ub])
	}
	pipRG_region[k]<-sum(sp)
}

## read gwa results for GB
resGB<-fread("output/pips_o_poly_bslmm_GR8_ph2_ch0.param.txt",header=FALSE)

pips<-as.numeric(unlist(resGB[,5]))
L<-length(pips)
sigGB<-which(pips > 0.00001) 
bnds<-c(which(resGB[-L,1] != resGB[-1,1])+.5,L+.5)
mids<-tapply(X=1:L,INDEX=resGB[,1],mean)


## regions from single res
reg_green<-read.csv("ChumashColorRegions-green_genome.csv")
r1lb<-reg_green[1,6]
r1ub<-reg_green[1,7]
r2lb<-reg_green[1,8]
r2ub<-reg_green[1,9]
r3lb<-reg_green[1,10]
r3ub<-reg_green[1,11]

sum(pips[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r1ub])
#[1] 1.084871
sum(pips[resGB[,1]==4 & resGB[,2] > r2lb & resGB[,2] < r2ub])
#[1] 1.071337
sum(pips[resGB[,1]==4 & resGB[,2] > r3lb & resGB[,2] < r3ub])
#[1] 0.2200295

sum(pips[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r3ub])
#3.244835
pipsGB<-pips

nr<-dim(reg_green)[1]
pipGB_region<-rep(NA,nr)
for(k in 1:nr){
        sp<-rep(NA,reg_green$num_windows[k])
        for(a in 1:reg_green$num_windows[k]){
                b<-(a-1)*2
                lb<-reg_green[k,6+b]
                ub<-reg_green[k,7+b]
                sp[a]<-sum(pips[resGB[,1]==4 & resGB[,2] > lb & resGB[,2] < ub])
        }
        pipGB_region[k]<-sum(sp)
}



## number of QTL total
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r1ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r1ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r1ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r1ub])
#[1] 1.94864 
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r2lb & resGB[,2] < r2ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r2lb & resRG[,2] < r2ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r2lb & resGB[,2] < r2ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r2lb & resRG[,2] < r2ub])
#[1] 2.122052
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r3lb & resGB[,2] < r3ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r3lb & resRG[,2] < r3ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r3lb & resGB[,2] < r3ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r3lb & resRG[,2] < r3ub])
#[1] 0.2210361
sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r3ub]) + sum(pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r3ub]) - sum(pipsGB[resGB[,1]==4 & resGB[,2] > r1lb & resGB[,2] < r3ub] * pipsRG[resRG[,1]==4 & resRG[,2] > r1lb & resRG[,2] < r3ub])
#[1] 5.258634


summary(pipRG_region)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1.334   1.876   1.914   1.856   1.915   1.916  
summary(pipGB_region)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.308   2.121   2.367   2.451   3.151   3.205 

