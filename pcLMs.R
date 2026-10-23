## which PCs exaplain color
library(BMS)

G<-fread("hsub_mlpntest_tchum_green.geno",header=FALSE,sep=",")
gmat<-as.matrix(G[,-c(1:3)])

hsubSNPs_green<-fread("hsubSNPs_green.txt",header=FALSE)
snp_green<-cbind(as.numeric(unlist(hsubSNPs_green[,1])),as.numeric(unlist(hsubSNPs_green[,2])))
ch<-snp_green[,1]
pos<-snp_green[,2]

ch8<-which(ch==4)
gmat_ch8<-gmat[ch8,]
pos_ch8<-pos[ch8]

trt<-read.table("sort_green_phfull_tchum2019.txt",header=TRUE,comment.char="%",sep=",")

regions<-read.csv("ChumashColorRegions-green_genome.csv")
rr<-1
rmat<-matrix(regions[rr,6:(6+regions$num_windows[rr]*2-1)],byrow=TRUE,ncol=2)
rmat
#[1,] 60040001 61640001
#[2,] 62480001 63860001
#[3,] 64490001 64890001


r1<-which(pos_ch8 > rmat[1,1] & pos_ch8 < rmat[1,2])
r2<-which(pos_ch8 > rmat[2,1] & pos_ch8 < rmat[2,2])
r3<-which(pos_ch8 > rmat[3,1] & pos_ch8 < rmat[3,2])

o1<-prcomp(t(gmat_ch8[r1,]),center=TRUE,scale=FALSE)
o2<-prcomp(t(gmat_ch8[r2,]),center=TRUE,scale=FALSE)
o3<-prcomp(t(gmat_ch8[r3,]),center=TRUE,scale=FALSE)

## drop missing
nas<-which(is.na(trt$Survival)==TRUE)

pcs<-cbind(o1$x[-nas,1:2],o2$x[-nas,1:2],o3$x[-nas,1:2])
nms<-c("y","R1PC1","R1PC2","R2PC1","R2PC2","R3PC1","R3PC2")
dfRG<-data.frame(y=trt$RG[-nas],pcs)
dfGB<-data.frame(y=trt$GB[-nas],pcs)
colnames(dfRG)<-nms
colnames(dfGB)<-nms

mnpcc<-rep(NA,6)
pcc<-cor(pcs)
mnpcc[1]<-mean(abs(pcc[1,3:6]))
mnpcc[2]<-mean(abs(pcc[2,3:6]))
mnpcc[3]<-mean(abs(pcc[3,c(1,2,5,6)]))
mnpcc[4]<-mean(abs(pcc[4,c(1,2,5,6)]))
mnpcc[5]<-mean(abs(pcc[5,1:4]))
mnpcc[6]<-mean(abs(pcc[6,1:4]))


## univariate correlations
corsRG<-matrix(NA,nrow=6,ncol=4)
corsGB<-matrix(NA,nrow=6,ncol=4)
for(i in 1:6){
	o<-cor.test(trt$RG[-nas],pcs[,i])
	corsRG[i,]<-c(o$estimate,o$conf.int[1:2],o$p.value)
	o<-cor.test(trt$GB[-nas],pcs[,i])
	corsGB[i,]<-c(o$estimate,o$conf.int[1:2],o$p.value)
}
corsRG
#            [,1]       [,2]        [,3]         [,4]
#[1,] -0.42465059 -0.49755790 -0.34578561 3.972947e-21
#[2,]  0.24930860  0.16057059  0.33404868 8.379020e-08
#[3,] -0.49217631 -0.55917456 -0.41879118 7.877398e-29
#[4,]  0.03930318 -0.05332904  0.13126475 4.055507e-01
#[5,] -0.45483963 -0.52519638 -0.37830696 2.325958e-24
#[6,] -0.11664981 -0.20685775 -0.02447526 1.328293e-02

corsGB
#[1,]  0.1948231  0.1042623  0.2821798 3.164863e-05
#[2,] -0.3415276 -0.4206849 -0.2572093 9.339767e-14
#[3,]  0.2864356  0.1992734  0.3691011 6.036553e-10
#[4,] -0.3336714 -0.4133602 -0.2489104 3.651764e-13
#[5,]  0.2703182  0.1824385  0.3539132 5.636493e-09
#[6,]  0.3233355  0.2380108  0.4037077 2.071155e-12


## bma

mfsRG<-bms(dfRG,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
mfsGB<-bms(dfGB,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)

summary(mfsRG)
#Mean no. regressors               Draws             Burnins                Time 
#           "2.1585"                "64"                 "0"   "0.02771115 secs" 
# No. models visited      Modelspace 2^K           % visited         % Topmodels 
#               "64"                "64"               "100"               "100" 
#           Corr PMP            No. Obs.         Model Prior             g-Prior 
#               "NA"               "450"       "uniform / 3"               "UIP" 
#    Shrinkage-Stats 
#        "Av=0.9978" 

ooRG<-coef(mfsRG,order.by.pip=FALSE)
#	PIP     Post Mean      Post SD Cond.Pos.Sign Idx
#R1PC1 0.10557055  3.680684e-05 1.726101e-04     0.9614042   1
#R1PC2 0.23207646  9.269051e-05 1.967439e-04     1.0000000   2
#R2PC1 0.99589680 -1.383697e-03 2.534039e-04     0.0000000   3
#R2PC2 0.05463332 -1.237957e-06 4.750994e-05     0.3328252   4
#R3PC1 0.07269108 -2.003976e-05 1.033023e-04     0.0000000   5
#R3PC2 0.69759240 -4.579997e-04 3.531630e-04     0.0000000   6

summary(mfsGB)
#Mean no. regressors               Draws             Burnins                Time 
#           "3.2046"                "64"                 "0"   "0.01127553 secs" 
# No. models visited      Modelspace 2^K           % visited         % Topmodels 
#               "64"                "64"               "100"               "100" 
#           Corr PMP            No. Obs.         Model Prior             g-Prior 
#               "NA"               "450"       "uniform / 3"               "UIP" 
#    Shrinkage-Stats 
#        "Av=0.9978" 


ooGB<-coef(mfsGB,order.by.pip=FALSE)
#	PIP     Post Mean      Post SD Cond.Pos.Sign Idx
#R1PC1 0.99967053 -1.422472e-03 2.478468e-04    0.00000180   1
#R1PC2 0.05752162  9.554961e-06 6.998281e-05    0.99678585   2
#R2PC1 0.99999660  2.325828e-03 3.099910e-04    1.00000000   3
#R2PC2 0.99999999 -1.236374e-03 1.445817e-04    0.00000000   4
#R3PC1 0.04643686 -2.786291e-06 5.169123e-05    0.00003779   5
#R3PC2 0.10099500  3.291111e-05 1.260603e-04    1.00000000   6

pdf("figPCPips.pdf",width=7,height=4)
par(mar=c(4.5,4.5,1,1))
barplot(rbind(ooRG[,1],ooGB[,1]),beside=TRUE,col=c("firebrick","royalblue"),names.arg=c("R1.PC1","R1.PC2","R2.PC1","R2.PC2","R3.PC1","R3.PC2"),ylab="PIP",cex.lab=1.3,ylim=c(0,1),xlab="Covariate")
box()
dev.off()

