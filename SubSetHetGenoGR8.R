## gets rid of extreme het loci from the genotype file
## also drop large number of individuals at the end with no trait data
library(data.table)

dat<-fread("mlpntest_GR8.06_tchum.geno",header=FALSE,sep=",")

dat<-dat[,-559]

ph<-read.table("ph_tchumGR8.txt",header=FALSE)



het<-as.matrix(fread("HetAndP_color_GR8.txt",header=TRUE))

keep<-which(het[,3] > -.2 & het[,3] < .2)
kdat<-dat[keep,1:(457+3)]
write.table(file="hsub_mlpntest_GR8_tchum.geno",kdat,row.names=FALSE,col.names=FALSE,quote=FALSE,sep=",")
write.table(file="sub_ph_tchumGR8.txt",ph[1:457,],row.names=FALSE,col.names=FALSE,quote=FALSE)
