## get p and het for individuals with color dat
library(data.table)

dat<-fread("mlpntest_EXP_tchum.geno",header=FALSE,sep=",")

g<-as.matrix(dat[,-c(1:3,433)]) ## last is trailing space

## keep inds with color data
ph<-read.table("ph_tchum2019.txt")
keep<-which(is.na(ph[,1])==FALSE)

het<-apply(g[,keep]==1,1,mean)
## all hets for some SNPs
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.04245 0.09670 0.20591 0.28774 1.00000 
sum(het==1)
#[1] 101


p<-apply(g[,keep],1,mean)/2
ehet<-2*p*(1-p)
hdev<-(het-ehet)

o<-cbind(het,p,hdev)

write.table(o,file="HetAndP_color.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)
