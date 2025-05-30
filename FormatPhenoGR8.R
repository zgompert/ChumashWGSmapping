## order color for GR8

## ids in order of genetic data
ids<-read.table("indIdsGR8.txt",header=FALSE)

## set up data table
N<-dim(ids)[1] #555
ph<-matrix(NA,nrow=N,ncol=2)

## trt data
trt<-read.table("2015_Tchumash_color.csv",header=TRUE,sep=",")
for(i in 1:N){
	a<-which(trt$id==ids[i,1])
	if(length(a)==1){
		ph[i,]<-c(trt$latRG[a],trt$latGB[a])
	}
}
apply(ph,2,mean,na.rm=TRUE)
#[1] 0.004728719 0.247120000

write.table(file="ph_tchumGR8.txt",ph,row.names=FALSE,col.names=FALSE,quote=FALSE)
