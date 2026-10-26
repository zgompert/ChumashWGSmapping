## get p and het for individuals with color dat
library(data.table)

####### green
dat<-fread("mlpntest_GR8_tchum_gr8_green.geno",header=FALSE,sep=",")

g<-as.matrix(dat[,-c(1:3)])

## keep inds with color data
ph<-read.table("ph_tchumGR8.txt")
keep<-which(is.na(ph[,1])==FALSE)

het<-apply(g[,keep]==1,1,mean)
summary(het)
sum(het==1)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.04933 0.12108 0.17329 0.27578 1.00000 
#[1] 6


p<-apply(g[,keep],1,mean)/2
ehet<-2*p*(1-p)
hdev<-(het-ehet)

o<-cbind(het,p,hdev)

write.table(o,file="HetAndP_GR8.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)


