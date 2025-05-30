## get p and het for individuals with color dat
library(data.table)

dat<-fread("mlpntest_GR8.06_tchum.geno",header=FALSE,sep=",")

g<-as.matrix(dat[,-c(1:3,559)]) ## last is trailing space

## keep inds with color data
ph<-read.table("ph_tchumGR8.txt",header=FALSE)
keep<-which(is.na(ph[,1])==FALSE) #keep 455

het<-apply(g[,keep]==1,1,mean)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.04245 0.09670 0.20591 0.28774 1.00000 


p<-apply(g[,keep],1,mean)/2
ehet<-2*p*(1-p)
hdev<-(het-ehet)

o<-cbind(het,p,hdev)

write.table(o,file="HetAndP_color_GR8.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)
