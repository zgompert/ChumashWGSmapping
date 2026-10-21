## order color and surival data

## ids in order of genetic data
ids<-read.table("indIds_green.txt",header=FALSE)

## set up data table
N<-dim(ids)[1] #455
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
#[1] 0.01955499 0.17498849 0.21476510
apply(is.na(ph)==FALSE,2,sum)
#[1] 450 450 447

write.table(file="ph_tchum_green.txt",ph,row.names=FALSE,col.names=FALSE,quote=FALSE)

## ids in order of genetic data
ids<-read.table("indIds_red.txt",header=FALSE)

## set up data table
N<-dim(ids)[1] #455
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
#[1] 0.01955499 0.17498849 0.21476510
apply(is.na(ph)==FALSE,2,sum)
#[1] 450 450 447

write.table(file="ph_tchum_red.txt",ph,row.names=FALSE,col.names=FALSE,quote=FALSE)
