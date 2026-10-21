## gets rid of extreme het loci from the genotype file
library(data.table)

dat<-fread("mlpntest_tchum_red.geno",header=FALSE,sep=",")


het<-as.matrix(fread("HetAndP_red.txt",header=TRUE))

keep<-which(het[,3] > -.2 & het[,3] < .2)
kdat<-dat[keep,]
write.table(file="hsub_mlpntest_tchum_red.geno",kdat,row.names=FALSE,col.names=FALSE,quote=FALSE,sep=",")

dat<-fread("mlpntest_tchum_green.geno",header=FALSE,sep=",")


het<-as.matrix(fread("HetAndP_green.txt",header=TRUE))

keep<-which(het[,3] > -.2 & het[,3] < .2)
kdat<-dat[keep,]
write.table(file="hsub_mlpntest_tchum_green.geno",kdat,row.names=FALSE,col.names=FALSE,quote=FALSE,sep=",")
