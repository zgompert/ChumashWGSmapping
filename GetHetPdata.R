## get p and het for individuals with color dat
library(data.table)

####### green
dat<-fread("mlpntest_tchum_green.geno",header=FALSE,sep=",")

g<-as.matrix(dat[,-c(1:3)])

## keep inds with color data
ph<-read.table("ph_tchum_green.txt")
keep<-which(is.na(ph[,1])==FALSE)

het<-apply(g[,keep]==1,1,mean)
## all hets for some SNPs
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.04889 0.11778 0.18418 0.28444 0.99556
sum(het==1)
#[1] 0


p<-apply(g[,keep],1,mean)/2
ehet<-2*p*(1-p)
hdev<-(het-ehet)

o<-cbind(het,p,hdev)

write.table(o,file="HetAndP_green.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)


####### red ###########
dat<-fread("mlpntest_tchum_red.geno",header=FALSE,sep=",")

g<-as.matrix(dat[,-c(1:3)])

## keep inds with color data
ph<-read.table("ph_tchum_red.txt")
keep<-which(is.na(ph[,1])==FALSE)

het<-apply(g[,keep]==1,1,mean)
## all hets for some SNPs
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.04889 0.11556 0.18299 0.28222 1.00000
sum(het==1)
#[1] 9 


p<-apply(g[,keep],1,mean)/2
ehet<-2*p*(1-p)
hdev<-(het-ehet)

o<-cbind(het,p,hdev)

write.table(o,file="HetAndP_red.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)
