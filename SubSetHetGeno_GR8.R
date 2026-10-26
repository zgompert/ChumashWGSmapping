## gets rid of extreme het loci from the genotype file
library(data.table)

dat<-fread("mlpntest_GR8_tchum_gr8_green.geno",header=FALSE,sep=",")


het<-as.matrix(fread("HetAndP_GR8.txt",header=TRUE))

keep<-which(het[,3] > -.2 & het[,3] < .2)
kdat<-dat[keep,]
write.table(file="hsub_mlpntest_GR8_tchum_gr8_green.geno",kdat,row.names=FALSE,col.names=FALSE,quote=FALSE,sep=",")

