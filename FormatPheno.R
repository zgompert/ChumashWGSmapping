## order color and surival data

## ids in order of genetic data
ids<-read.table("indIds.txt",header=FALSE)

## set up data table
N<-dim(ids)[1] #429
ph<-matrix(NA,nrow=N,ncol=3)

## trt data
trt<-read.table("2019_Tchumash_transplant_table.csv",header=TRUE,sep=",")
for(i in 1:N){
	a<-which(trt$id==ids[i,1])
	if(length(a)==1){
		ph[i,]<-c(trt$RG[a],trt$GB[a],trt$Survival[a])
	}else if(length(a)==2){
		ph[i,]<-c(trt$RG[a[1]],trt$GB[a[1]],NA)
	}
}
apply(ph,2,mean,na.rm=TRUE)
#[1] 0.01816416 0.17633720 0.22090261
apply(is.na(ph)==FALSE,2,sum)
#[1] 424 424 421

write.table(file="ph_tchum2019.txt",ph,row.names=FALSE,col.names=FALSE,quote=FALSE)
