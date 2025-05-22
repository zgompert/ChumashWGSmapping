## gets rid of extreme het loci from the genotype file
library(data.table)

dat<-fread("mlpntest_EXP_tchum.geno",header=FALSE,sep=",")

dat<-dat[,-433]



het<-as.matrix(fread("HetAndP_color.txt",header=TRUE))

keep<-which(het[,3] > -.2 & het[,3] < .2)
kdat<-dat[keep,]
write.table(file="hsub_mlpntest_EXP_tchum.geno",kdat,row.names=FALSE,col.names=FALSE,quote=FALSE,sep=",")
