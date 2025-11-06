## chumash edinburgh
library(data.table)

### old chumash regions
## r1
r1lb<-55676246
r1ub<-57380207

## r2
r2lb<-58460759
r2ub<-59734122
## r3
r3lb<-62179278
r3ub<-62791657


## first compare new genomes to old chumash

## 24_0159 H1
dat<-fread("out_synteny_Ed24_0159_h1.psl",header=FALSE)
dfdat<-as.data.frame(dat)
## identify large scaffolds for each
xx<-table(dfdat[,14]) ## dovetail
dove<-names(xx)[xx>500] ## 10 
xx<-table(dfdat[,10]) ## 24_0159
chum<-names(xx)[xx>500] ## 11
keep<-(dfdat[,14] %in% dove) & (dfdat[,10] %in% chum)
subDfdat<-dfdat[keep,]


tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
dove_sc<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,40,4)])
tchum_sc<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,22,2)],pattern="Chr",replacement=""))

## normalize with respect to dove 

ntab<-tab
for(i in 1:10){
	ntab[,i]<-ntab[,i]/sum(ntab[,i])
}


pdf("SynTchum_ed24_0159h1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. chumash (24_0159 H1)",ylab="T. chumash (Dovetail)",cex.lab=1.4)
axis(2,at=seq(0,10,length.out=10)/10,dove_sc,las=2)
axis(1,at=seq(0,11,length.out=11)/11,tchum_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes, chromosome numbers here based on knulli
chtab<-matrix(c(1,43,1, ## 1 vs 4 vs 9 (dropped) arbitrary
	2,1392,6,
	4,43,2,
	5,56,5,
	6,1469,7,
	7,1510,10,
	8,113,4,
	10,1213,8,
	11,48,9,
	12,1403,11,
	13,1308,3),nrow=11,ncol=3,byrow=TRUE)


chumSize<-c(926504516,142035733,926504516,149761397,151976569,125220160,174643345,135538230,129098588,79636963,208324082)

pdf("AlnPlot_Tchum_ed24_0159h1.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:11){
	tdove<-grep(x=dfdat[,14],pattern=chtab[i,2]) 
	tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
	cc<-tdove[tdove %in% tch]
	subd<-dfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. chumash (24_0159 H1)",ylab="T. chumash (Dovetail)")
	title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],chumSize[i]-subd[j,16:17],col="cadetblue")
		}
	}

}
dev.off()

## get regions
i<-7## ch8
tdove<-grep(x=dfdat[,14],pattern=chtab[i,2])
tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
cc<-tdove[tdove %in% tch]
subd<-dfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])
N<-dim(subd)[1]
xxx<-subd[,12:13]
yyy<-subd[,16:17]
yyy[subd[,9]=="+-",]<-chumSize[i]-subd[subd[,9]=="+-",16:17]
## r1
r1sz<-r1ub-r1lb
r1<-which((yyy[,2] > r1lb & yyy[,2] < r1ub) | (yyy[,1] > r1lb & yyy[,1] < r1ub))
midx<-median(as.matrix(xxx[r1,]))
dx<-abs(midx-apply(xxx[r1,],1,mean))
keep<-which(dx < r1sz*2)
r1_24_0159h1<-quantile(as.vector(as.matrix(xxx[r1[keep],])),probs=c(0.025,.975))
## r2
r2sz<-r2ub-r2lb
r2<-which((yyy[,2] > r2lb & yyy[,2] < r2ub) | (yyy[,1] > r2lb & yyy[,1] < r2ub))
midx<-median(as.matrix(xxx[r2,]))
dx<-abs(midx-apply(xxx[r2,],1,mean))
keep<-which(dx < r2sz*2)
r2_24_0159h1<-quantile(as.vector(as.matrix(xxx[r2[keep],])),probs=c(0.025,.975))
## r3
r3sz<-r3ub-r3lb
r3<-which((yyy[,2] > r3lb & yyy[,2] < r3ub) | (yyy[,1] > r3lb & yyy[,1] < r3ub))
midx<-median(as.matrix(xxx[r2,]))
dx<-abs(midx-apply(xxx[r3,],1,mean))
keep<-which(dx < r3sz*2)
r3_24_0159h1<-quantile(as.vector(as.matrix(xxx[r3[keep],])),probs=c(0.025,.975))

## 24_0159 H2
dat<-fread("out_synteny_Ed24_0159_h2.psl",header=FALSE)
dfdat<-as.data.frame(dat)
## identify large scaffolds for each
xx<-table(dfdat[,14]) ## dovetail
dove<-names(xx)[xx>500] ## 10 
xx<-table(dfdat[,10]) ## 24_0159
chum<-names(xx)[xx>500] ## 11
keep<-(dfdat[,14] %in% dove) & (dfdat[,10] %in% chum)
subDfdat<-dfdat[keep,]


tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
dove_sc<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,40,4)])
tchum_sc<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,22,2)],pattern="Chr",replacement=""))

## normalize with respect to dove 

ntab<-tab
for(i in 1:10){
	ntab[,i]<-ntab[,i]/sum(ntab[,i])
}


pdf("SynTchum_ed24_0159h2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. chumash (24_0159 H2)",ylab="T. chumash (Dovetail)",cex.lab=1.4)
axis(2,at=seq(0,10,length.out=10)/10,dove_sc,las=2)
axis(1,at=seq(0,11,length.out=11)/11,tchum_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes, chromosome numbers here based on knulli
chtab<-matrix(c(1,43,1, ## 1 vs 4 vs 9 (dropped) arbitrary
	2,1392,6,
	4,43,2,
	5,56,5,
	6,1469,7,
	7,1510,10,
	8,113,4,
	10,1213,8,
	11,48,9,
	12,1403,11,
	13,1308,3),nrow=11,ncol=3,byrow=TRUE)


chumSize<-c(926504516,142035733,926504516,149761397,151976569,125220160,174643345,135538230,129098588,79636963,208324082)

pdf("AlnPlot_Tchum_ed24_0159h2.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:11){
	tdove<-grep(x=dfdat[,14],pattern=chtab[i,2]) 
	tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
	cc<-tdove[tdove %in% tch]
	subd<-dfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. chumash (24_0159 H2)",ylab="T. chumash (Dovetail)")
	title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],chumSize[i]-subd[j,16:17],col="cadetblue")
		}
	}

}
dev.off()

## get regions
i<-7## ch8
tdove<-grep(x=dfdat[,14],pattern=chtab[i,2])
tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
cc<-tdove[tdove %in% tch]
subd<-dfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])
N<-dim(subd)[1]
xxx<-subd[,12:13]
yyy<-subd[,16:17]
yyy[subd[,9]=="+-",]<-chumSize[i]-subd[subd[,9]=="+-",16:17]
## r1
r1sz<-r1ub-r1lb
r1<-which((yyy[,2] > r1lb & yyy[,2] < r1ub) | (yyy[,1] > r1lb & yyy[,1] < r1ub))
midx<-median(as.matrix(xxx[r1,]))
dx<-abs(midx-apply(xxx[r1,],1,mean))
keep<-which(dx < r1sz*2)
r1_24_0159h2<-quantile(as.vector(as.matrix(xxx[r1[keep],])),probs=c(0.025,.975))
## r2
r2sz<-r2ub-r2lb
r2<-which((yyy[,2] > r2lb & yyy[,2] < r2ub) | (yyy[,1] > r2lb & yyy[,1] < r2ub))
midx<-median(as.matrix(xxx[r2,]))
dx<-abs(midx-apply(xxx[r2,],1,mean))
keep<-which(dx < r2sz*2)
r2_24_0159h2<-quantile(as.vector(as.matrix(xxx[r2[keep],])),probs=c(0.025,.975))
## r3
r3sz<-r3ub-r3lb
r3<-which((yyy[,2] > r3lb & yyy[,2] < r3ub) | (yyy[,1] > r3lb & yyy[,1] < r3ub))
midx<-median(as.matrix(xxx[r2,]))
dx<-abs(midx-apply(xxx[r3,],1,mean))
keep<-which(dx < r3sz*2)
r3_24_0159h2<-quantile(as.vector(as.matrix(xxx[r3[keep],])),probs=c(0.025,.975))

## 24_0161 H1
dat<-fread("out_synteny_Ed24_0161_h1.psl",header=FALSE)
dfdat<-as.data.frame(dat)
## identify large scaffolds for each
xx<-table(dfdat[,14]) ## dovetail
dove<-names(xx)[xx>500] ## 10 
xx<-table(dfdat[,10]) ## 24_0161
chum<-names(xx)[xx>500] ## 11
keep<-(dfdat[,14] %in% dove) & (dfdat[,10] %in% chum)
subDfdat<-dfdat[keep,]


tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
dove_sc<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,40,4)])
tchum_sc<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,22,2)],pattern="Chr",replacement=""))

## normalize with respect to dove 

ntab<-tab
for(i in 1:10){
	ntab[,i]<-ntab[,i]/sum(ntab[,i])
}


pdf("SynTchum_ed24_0161h1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. chumash (24_0161 H1)",ylab="T. chumash (Dovetail)",cex.lab=1.4)
axis(2,at=seq(0,10,length.out=10)/10,dove_sc,las=2)
axis(1,at=seq(0,11,length.out=11)/11,tchum_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes, chromosome numbers here based on knulli
chtab<-matrix(c(1,43,1, ## 1 vs 4 vs 9 (dropped) arbitrary
	2,1392,7,
	4,43,2,
	5,56,5,
	6,1469,6,
	7,1510,9,
	8,113,4,
	10,1213,8,
	11,48,10,
	12,1403,11,
	13,1308,3),nrow=11,ncol=3,byrow=TRUE)


chumSize<-c(926504516,142035733,926504516,149761397,151976569,125220160,174643345,135538230,129098588,79636963,208324082)

pdf("AlnPlot_Tchum_ed24_0161h1.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:11){
	tdove<-grep(x=dfdat[,14],pattern=chtab[i,2]) 
	tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
	cc<-tdove[tdove %in% tch]
	subd<-dfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. chumash (24_0161 H1)",ylab="T. chumash (Dovetail)")
	title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],chumSize[i]-subd[j,16:17],col="cadetblue")
		}
	}
	if(chtab[i,1]==8){
		abline(h=c(55676246,62791657),lty=2)
	}

}
dev.off()

## get regions
i<-7## ch8
tdove<-grep(x=dfdat[,14],pattern=chtab[i,2])
tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
cc<-tdove[tdove %in% tch]
subd<-dfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])
N<-dim(subd)[1]
xxx<-subd[,12:13]
yyy<-subd[,16:17]
yyy[subd[,9]=="+-",]<-chumSize[i]-subd[subd[,9]=="+-",16:17]
## r1
r1sz<-r1ub-r1lb
r1<-which((yyy[,2] > r1lb & yyy[,2] < r1ub) | (yyy[,1] > r1lb & yyy[,1] < r1ub))
midx<-median(as.matrix(xxx[r1,]))
dx<-abs(midx-apply(xxx[r1,],1,mean))
keep<-which(dx < r1sz*2)
r1_24_0161h1<-quantile(as.vector(as.matrix(xxx[r1[keep],])),probs=c(0.025,.975))
## r2
r2sz<-r2ub-r2lb
r2<-which((yyy[,2] > r2lb & yyy[,2] < r2ub) | (yyy[,1] > r2lb & yyy[,1] < r2ub))
midx<-median(as.matrix(xxx[r2,]))
dx<-abs(midx-apply(xxx[r2,],1,mean))
keep<-which(dx < r2sz*2)
r2_24_0161h1<-quantile(as.vector(as.matrix(xxx[r2[keep],])),probs=c(0.025,.975))
## r3
r3sz<-r3ub-r3lb
r3<-which((yyy[,2] > r3lb & yyy[,2] < r3ub) | (yyy[,1] > r3lb & yyy[,1] < r3ub))
midx<-median(as.matrix(xxx[r2,]))
dx<-abs(midx-apply(xxx[r3,],1,mean))
keep<-which(dx < r3sz*2)
r3_24_0161h1<-quantile(as.vector(as.matrix(xxx[r3[keep],])),probs=c(0.025,.975))

i<-7
pdf("AlnPlot_Tchum_ed24_0161h1-zoom.pdf",width=5,height=5)
par(mfrow=c(1,1))
par(mar=c(4.5,5.5,2.5,1.5))
tdove<-grep(x=dfdat[,14],pattern=chtab[i,2])
tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
cc<-tdove[tdove %in% tch]
subd<-dfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])

plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(35000000,48000000),ylim=c(50000000,65000000),cex.lab=1.4,xlab="T. chumash (24_0161 H1)",ylab="T. chumash (Dovetail)")
title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
N<-dim(subd)[1]
for(j in 1:N){
	if(subd[j,9]=="++"){
        	lines(subd[j,12:13],subd[j,16:17])
        }
        else{
        	lines(subd[j,12:13],chumSize[i]-subd[j,16:17],col="cadetblue")
        }
}
abline(h=c(55676246,62791657),lty=2)
dev.off()

## 24_0161 H2
dat<-fread("out_synteny_Ed24_0161_h2.psl",header=FALSE)
dfdat<-as.data.frame(dat)
## identify large scaffolds for each
xx<-table(dfdat[,14]) ## dovetail
dove<-names(xx)[xx>500] ## 10 
xx<-table(dfdat[,10]) ## 24_0161
chum<-names(xx)[xx>500] ## 11
keep<-(dfdat[,14] %in% dove) & (dfdat[,10] %in% chum)
subDfdat<-dfdat[keep,]


tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
dove_sc<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,40,4)])
tchum_sc<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,22,2)],pattern="Chr",replacement=""))

## normalize with respect to dove 

ntab<-tab
for(i in 1:10){
	ntab[,i]<-ntab[,i]/sum(ntab[,i])
}


pdf("SynTchum_ed24_0161h2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. chumash (24_0161 H2)",ylab="T. chumash (Dovetail)",cex.lab=1.4)
axis(2,at=seq(0,10,length.out=10)/10,dove_sc,las=2)
axis(1,at=seq(0,11,length.out=11)/11,tchum_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes, chromosome numbers here based on knulli
chtab<-matrix(c(1,43,1, ## 1 vs 4 vs 9 (dropped) arbitrary
	2,1392,7,
	4,43,2,
	5,56,5,
	6,1469,6,
	7,1510,9,
	8,113,4,
	10,1213,8,
	11,48,10,
	12,1403,11,
	13,1308,3),nrow=11,ncol=3,byrow=TRUE)


chumSize<-c(926504516,142035733,926504516,149761397,151976569,125220160,174643345,135538230,129098588,79636963,208324082)

pdf("AlnPlot_Tchum_ed24_0161h2.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:11){
	tdove<-grep(x=dfdat[,14],pattern=chtab[i,2]) 
	tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
	cc<-tdove[tdove %in% tch]
	subd<-dfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. chumash (24_0161 H2)",ylab="T. chumash (Dovetail)")
	title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],chumSize[i]-subd[j,16:17],col="cadetblue")
		}
	}
	if(chtab[i,1]==8){
		abline(h=c(55676246,62791657),lty=2)
	}

}
dev.off()

## get regions
i<-7## ch8
tdove<-grep(x=dfdat[,14],pattern=chtab[i,2])
tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
cc<-tdove[tdove %in% tch]
subd<-dfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])
N<-dim(subd)[1]
xxx<-subd[,12:13]
yyy<-subd[,16:17]
yyy[subd[,9]=="+-",]<-chumSize[i]-subd[subd[,9]=="+-",16:17]
## r1
r1sz<-r1ub-r1lb
r1<-which((yyy[,2] > r1lb & yyy[,2] < r1ub) | (yyy[,1] > r1lb & yyy[,1] < r1ub))
midx<-median(as.matrix(xxx[r1,]))
dx<-abs(midx-apply(xxx[r1,],1,mean))
keep<-which(dx < r1sz*2)
r1_24_0161h2<-quantile(as.vector(as.matrix(xxx[r1[keep],])),probs=c(0.025,.975))
## r2
r2sz<-r2ub-r2lb
r2<-which((yyy[,2] > r2lb & yyy[,2] < r2ub) | (yyy[,1] > r2lb & yyy[,1] < r2ub))
midx<-median(as.matrix(xxx[r2,]))
dx<-abs(midx-apply(xxx[r2,],1,mean))
keep<-which(dx < r2sz*2)
r2_24_0161h2<-quantile(as.vector(as.matrix(xxx[r2[keep],])),probs=c(0.025,.975))
## r3
r3sz<-r3ub-r3lb
r3<-which((yyy[,2] > r3lb & yyy[,2] < r3ub) | (yyy[,1] > r3lb & yyy[,1] < r3ub))
midx<-median(as.matrix(xxx[r2,]))
dx<-abs(midx-apply(xxx[r3,],1,mean))
keep<-which(dx < r3sz*2)
r3_24_0161h2<-quantile(as.vector(as.matrix(xxx[r3[keep],])),probs=c(0.025,.975))


## 24_0163 H1
dat<-fread("out_synteny_Ed24_0163_h1.psl",header=FALSE)
dfdat<-as.data.frame(dat)
## identify large scaffolds for each
xx<-table(dfdat[,14]) ## dovetail
dove<-names(xx)[xx>500] ## 10 
xx<-table(dfdat[,10]) ## 24_0163
chum<-names(xx)[xx>500] ## 11
keep<-(dfdat[,14] %in% dove) & (dfdat[,10] %in% chum)
subDfdat<-dfdat[keep,]


tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
dove_sc<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,40,4)])
tchum_sc<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,22,2)],pattern="Chr",replacement=""))

## normalize with respect to dove 

ntab<-tab
for(i in 1:10){
	ntab[,i]<-ntab[,i]/sum(ntab[,i])
}


pdf("SynTchum_ed24_0163h1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. chumash (24_0163 H1)",ylab="T. chumash (Dovetail)",cex.lab=1.4)
axis(2,at=seq(0,10,length.out=10)/10,dove_sc,las=2)
axis(1,at=seq(0,11,length.out=11)/11,tchum_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes, chromosome numbers here based on knulli
chtab<-matrix(c(1,43,1, ## 1 vs 4 vs 9 (dropped) arbitrary
	2,1392,7,
	4,43,2,
	5,56,5,
	6,1469,6,
	7,1510,10,
	8,113,4,
	10,1213,8,
	11,48,9,
	12,1403,11,
	13,1308,3),nrow=11,ncol=3,byrow=TRUE)


chumSize<-c(926504516,142035733,926504516,149761397,151976569,125220160,174643345,135538230,129098588,79636963,208324082)

pdf("AlnPlot_Tchum_ed24_0163h1.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:11){
	tdove<-grep(x=dfdat[,14],pattern=chtab[i,2]) 
	tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
	cc<-tdove[tdove %in% tch]
	subd<-dfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. chumash (24_0163 H1)",ylab="T. chumash (Dovetail)")
	title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],chumSize[i]-subd[j,16:17],col="cadetblue")
		}
	}
	if(chtab[i,1]==8){
		abline(h=c(55676246,62791657),lty=2)
	}

}
dev.off()

## get regions
i<-7## ch8
tdove<-grep(x=dfdat[,14],pattern=chtab[i,2])
tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
cc<-tdove[tdove %in% tch]
subd<-dfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])
N<-dim(subd)[1]
xxx<-subd[,12:13]
yyy<-subd[,16:17]
yyy[subd[,9]=="+-",]<-chumSize[i]-subd[subd[,9]=="+-",16:17]
## r1
r1sz<-r1ub-r1lb
r1<-which((yyy[,2] > r1lb & yyy[,2] < r1ub) | (yyy[,1] > r1lb & yyy[,1] < r1ub))
midx<-median(as.matrix(xxx[r1,]))
dx<-abs(midx-apply(xxx[r1,],1,mean))
keep<-which(dx < r1sz*2)
r1_24_0163h1<-quantile(as.vector(as.matrix(xxx[r1[keep],])),probs=c(0.025,.975))
## r2
r2sz<-r2ub-r2lb
r2<-which((yyy[,2] > r2lb & yyy[,2] < r2ub) | (yyy[,1] > r2lb & yyy[,1] < r2ub))
midx<-median(as.matrix(xxx[r2,]))
dx<-abs(midx-apply(xxx[r2,],1,mean))
keep<-which(dx < r2sz*2)
r2_24_0163h1<-quantile(as.vector(as.matrix(xxx[r2[keep],])),probs=c(0.025,.975))
## r3
r3sz<-r3ub-r3lb
r3<-which((yyy[,2] > r3lb & yyy[,2] < r3ub) | (yyy[,1] > r3lb & yyy[,1] < r3ub))
midx<-median(as.matrix(xxx[r2,]))
dx<-abs(midx-apply(xxx[r3,],1,mean))
keep<-which(dx < r3sz*2)
r3_24_0163h1<-quantile(as.vector(as.matrix(xxx[r3[keep],])),probs=c(0.025,.975))

i<-7
pdf("AlnPlot_Tchum_ed24_0163h1-zoom.pdf",width=5,height=5)
par(mfrow=c(1,1))
par(mar=c(4.5,5.5,2.5,1.5))
tdove<-grep(x=dfdat[,14],pattern=chtab[i,2])
tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
cc<-tdove[tdove %in% tch]
subd<-dfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])

plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(54000000,69000000),ylim=c(50000000,65000000),cex.lab=1.4,xlab="T. chumash (24_016e H1)",ylab="T. chumash (Dovetail)")
title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
N<-dim(subd)[1]
for(j in 1:N){
	if(subd[j,9]=="++"){
        	lines(subd[j,12:13],subd[j,16:17])
        }
        else{
        	lines(subd[j,12:13],chumSize[i]-subd[j,16:17],col="cadetblue")
        }
}
abline(h=c(55676246,62791657),lty=2)
dev.off()

## 24_0163 H2
dat<-fread("out_synteny_Ed24_0163_h2.psl",header=FALSE)
dfdat<-as.data.frame(dat)
## identify large scaffolds for each
xx<-table(dfdat[,14]) ## dovetail
dove<-names(xx)[xx>500] ## 10 
xx<-table(dfdat[,10]) ## 24_0163
chum<-names(xx)[xx>500] ## 11
keep<-(dfdat[,14] %in% dove) & (dfdat[,10] %in% chum)
subDfdat<-dfdat[keep,]


tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
dove_sc<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,40,4)])
tchum_sc<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,22,2)],pattern="Chr",replacement=""))

## normalize with respect to dove 

ntab<-tab
for(i in 1:10){
	ntab[,i]<-ntab[,i]/sum(ntab[,i])
}


pdf("SynTchum_ed24_0163h2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. chumash (24_0163 H2)",ylab="T. chumash (Dovetail)",cex.lab=1.4)
axis(2,at=seq(0,10,length.out=10)/10,dove_sc,las=2)
axis(1,at=seq(0,11,length.out=11)/11,tchum_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes, chromosome numbers here based on knulli
chtab<-matrix(c(1,43,1, ## 1 vs 4 vs 9 (dropped) arbitrary
	2,1392,7,
	4,43,2,
	5,56,5,
	6,1469,6,
	7,1510,10,
	8,113,4,
	10,1213,8,
	11,48,9,
	12,1403,11,
	13,1308,3),nrow=11,ncol=3,byrow=TRUE)


chumSize<-c(926504516,142035733,926504516,149761397,151976569,125220160,174643345,135538230,129098588,79636963,208324082)

pdf("AlnPlot_Tchum_ed24_0163h2.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:11){
	tdove<-grep(x=dfdat[,14],pattern=chtab[i,2]) 
	tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
	cc<-tdove[tdove %in% tch]
	subd<-dfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. chumash (24_0163 H2)",ylab="T. chumash (Dovetail)")
	title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],chumSize[i]-subd[j,16:17],col="cadetblue")
		}
	}
	if(chtab[i,1]==8){
		abline(h=c(55676246,62791657),lty=2)
	}

}
dev.off()

## get regions
i<-7## ch8
tdove<-grep(x=dfdat[,14],pattern=chtab[i,2])
tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
cc<-tdove[tdove %in% tch]
subd<-dfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])
N<-dim(subd)[1]
xxx<-subd[,12:13]
yyy<-subd[,16:17]
yyy[subd[,9]=="+-",]<-chumSize[i]-subd[subd[,9]=="+-",16:17]
## r1
r1sz<-r1ub-r1lb
r1<-which((yyy[,2] > r1lb & yyy[,2] < r1ub) | (yyy[,1] > r1lb & yyy[,1] < r1ub))
midx<-median(as.matrix(xxx[r1,]))
dx<-abs(midx-apply(xxx[r1,],1,mean))
keep<-which(dx < r1sz*2)
r1_24_0163h2<-quantile(as.vector(as.matrix(xxx[r1[keep],])),probs=c(0.025,.975))
## r2
r2sz<-r2ub-r2lb
r2<-which((yyy[,2] > r2lb & yyy[,2] < r2ub) | (yyy[,1] > r2lb & yyy[,1] < r2ub))
midx<-median(as.matrix(xxx[r2,]))
dx<-abs(midx-apply(xxx[r2,],1,mean))
keep<-which(dx < r2sz*2)
r2_24_0163h2<-quantile(as.vector(as.matrix(xxx[r2[keep],])),probs=c(0.025,.975))
## r3
r3sz<-r3ub-r3lb
r3<-which((yyy[,2] > r3lb & yyy[,2] < r3ub) | (yyy[,1] > r3lb & yyy[,1] < r3ub))
midx<-median(as.matrix(xxx[r2,]))
dx<-abs(midx-apply(xxx[r3,],1,mean))
keep<-which(dx < r3sz*2)
r3_24_0163h2<-quantile(as.vector(as.matrix(xxx[r3[keep],])),probs=c(0.025,.975))

save(list=ls(),file="alnsTchum.rdat")

## plots
prepDat<-function(dat=adat,chtab=chtab,id1="Tchum1",id2="Tchum2",
		  bndtab1=NA,bndtab2=NA){
	dfdat<-as.data.frame(dat)
	## identify large scaffolds for each
	xx<-table(dfdat[,14]) ## 24_0163 H1
	dove<-names(xx)[xx>500] ## 10 
	xx<-table(dfdat[,10]) ## 24_0163 H2
	chum<-names(xx)[xx>500] ## 11
	keep<-(dfdat[,14] %in% dove) & (dfdat[,10] %in% chum)
	subDfdat<-dfdat[keep,]


	tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
	tchum_sc<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,22,2)],pattern="Chr",replacement=""))
	dove_sc<-as.numeric(gsub(x=unlist(strsplit(x=colnames(tab),split="_"))[seq(2,22,2)],pattern="Chr",replacement=""))

	tdove<-grep(x=dfdat[,14],pattern=chtab[2])
	tch<-grep(x=dfdat[,10],pattern=chtab[3])
	cc<-tdove[tdove %in% tch]
	subd<-dfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab=id1,ylab=id2)

	title(main=paste("Chrom.",chtab[1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
	        	lines(subd[j,12:13],subd[j,16:17])
       		}
        	else{
        		lines(subd[j,12:13],yub-subd[j,16:17])
        	}
	}
	xlb<-min(bndtab1[,1]);xub<-max(bndtab1[,2])
	ylb<-min(bndtab2[,1]);yub<-max(bndtab2[,2])

	abline(h=c(ylb,yub),lty=2) 
	abline(v=c(xlb,xub),lty=2) 
	
	#xub<-max(subd[,13]);yub<-max(subd[,17])
	xlb<-min(bndtab1[,1]-5e6)
	xub<-min(bndtab1[,2]+5e6)
	ylb<-min(bndtab2[,1]-5e6)
	yub<-min(bndtab2[,2]+5e6)
	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(xlb,xub),ylim=c(ylb,yub),cex.lab=1.4,xlab=id1,ylab=id2)
	xub<-max(subd[,13]);yub<-max(subd[,17])
	title(main=paste("Zoom chrom.",chtab[1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
	        	lines(subd[j,12:13],subd[j,16:17])
       		}
        	else{
        		lines(subd[j,12:13],yub-subd[j,16:17])
        	}
	}

	abline(v=as.vector(t(bndtab1)),col=rep(c("firebrick","cadetblue","forestgreen"),each=2),lty=2) 
	abline(h=as.vector(t(bndtab2)),col=rep(c("firebrick","cadetblue","forestgreen"),each=2),lty=2) 
}

## 24_0163 H1 (called dove below to keep code simple) vs 24_0163 H2
adat<-fread("out_synteny_E240163H1_E240163H2.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240163H1_E240163H2.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0163 H2)",id2="T. chumash (0163 H1)",
                  bndtab1=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2),bndtab2=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1))
dev.off()

## 24_0161 H1 vs 24_0161 H2 = waiting
adat<-fread("out_synteny_E240161H1_E240161H2.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240151H1_E240161H2.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0161 H2)",id2="T. chumash (0161 H1)",
                  bndtab1=rbind(r1_24_0161h2,r2_24_0161h2,r3_24_0161h2),bndtab2=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1))
dev.off()

## 24_0159 H1 vs 24_0159 H2
adat<-fread("out_synteny_E240159H1_E240159H2.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240159H1_E240159H2.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0159 H2)",id2="T. chumash (0159 H1)",
                  bndtab1=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2),bndtab2=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1))
dev.off()

## 24_0159 H1 vs 24_0161 H1
adat<-fread("out_synteny_E240159H1_E240161H1.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240159H1_E240161H1.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0161 H1)",id2="T. chumash (0159 H1)",
                  bndtab1=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1),bndtab2=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1))
dev.off()

## 24_0159 H1 vs 24_0161 H2
adat<-fread("out_synteny_E240159H1_E240161H2.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240159H1_E240161H2.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0161 H2)",id2="T. chumash (0159 H1)",
                  bndtab1=rbind(r1_24_0161h2,r2_24_0161h2,r3_24_0161h2),bndtab2=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1))
dev.off()

## 24_0159 H1 vs 24_0163 H1
adat<-fread("out_synteny_E240159H1_E240163H1.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240159H1_E240163H1.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0163 H1)",id2="T. chumash (0159 H1)",
                  bndtab1=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1),bndtab2=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1))
dev.off()

## 24_0159 H1 vs 24_0163 H2
adat<-fread("out_synteny_E240159H1_E240163H2.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240159H1_E240163H2.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0163 H2)",id2="T. chumash (0159 H1)",
                  bndtab1=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2),bndtab2=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1))
dev.off()

## 24_0159 H2 vs 24_0161 H1
adat<-fread("out_synteny_E240159H2_E240161H1.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240159H2_E240161H1.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0161 H1)",id2="T. chumash (0159 H2)",
                  bndtab1=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1),bndtab2=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2))
dev.off()

## 24_0159 H2 vs 24_0161 H2
adat<-fread("out_synteny_E240159H2_E240161H2.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240159H2_E240161H2.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0161 H2)",id2="T. chumash (0159 H2)",
                  bndtab1=rbind(r1_24_0161h2,r2_24_0161h2,r3_24_0161h2),bndtab2=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2))
dev.off()

## 24_0159 H2 vs 24_0163 H1
adat<-fread("out_synteny_E240159H2_E240163H1.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240159H2_E240163H1.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0163 H1)",id2="T. chumash (0159 H2)",
                  bndtab1=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1),bndtab2=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2))
dev.off()

## 24_0159 H2 vs 24_0163 H2 == colinear
adat<-fread("out_synteny_E240159H2_E240163H2.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240159H2_E240163H2.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0163 H2)",id2="T. chumash (0159 H2)",
                  bndtab1=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2),bndtab2=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2))
dev.off()

## 24_0161 H1 vs 24_0163 H1
adat<-fread("out_synteny_E240161H1_E240163H1.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240161H1_E240163H1.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0163 H1)",id2="T. chumash (0161 H1)",
                  bndtab1=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1),bndtab2=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1))
dev.off()

## 24_0161 H1 vs 24_0163 H2
adat<-fread("out_synteny_E240161H1_E240163H2.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240161H1_E240163H2.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0163 H2)",id2="T. chumash (0161 H1)",
                  bndtab1=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2),bndtab2=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1))
dev.off()

## 24_0161 H2 vs 24_0163 H1
adat<-fread("out_synteny_E240161H2_E240163H1.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240161H2_E240163H1.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0163 H1)",id2="T. chumash (0161 H2)",
                  bndtab1=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1),bndtab2=rbind(r1_24_0161h2,r2_24_0161h2,r3_24_0161h2))
dev.off()

## 24_0161 H2 vs 24_0163 H2
adat<-fread("out_synteny_E240161H2_E240163H2.psl",header=FALSE)
chtab<-c(8,4,4)
pdf("tchum_E240161H2_E240163H2.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2.5,1))
prepDat(dat=adat,chtab=chtab,id1="T. chumash (0163 H2)",id2="T. chumash (0161 H2)",
                  bndtab1=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2),bndtab2=rbind(r1_24_0161h2,r2_24_0161h2,r3_24_0161h2))
dev.off()
######################## initial comparative alignments #############
## 24_0163 H1 (called dove below to keep code simple) vs 24_0163 H2
dat<-fread("out_synteny_E240163H1_E240163H2.psl",header=FALSE)
dfdat<-as.data.frame(dat)
## identify large scaffolds for each
xx<-table(dfdat[,14]) ## 24_0163 H1
dove<-names(xx)[xx>500] ## 10 
xx<-table(dfdat[,10]) ## 24_0163 H2
chum<-names(xx)[xx>500] ## 11
keep<-(dfdat[,14] %in% dove) & (dfdat[,10] %in% chum)
subDfdat<-dfdat[keep,]


tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tchum_sc<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,22,2)],pattern="Chr",replacement=""))
dove_sc<-as.numeric(gsub(x=unlist(strsplit(x=colnames(tab),split="_"))[seq(2,22,2)],pattern="Chr",replacement=""))

## normalize with respect to dove 

ntab<-tab
for(i in 1:11){
	ntab[,i]<-ntab[,i]/sum(ntab[,i])
}


pdf("SynTchum_ed24_0163h1_24_0163h2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. chumash (24_0163 H2)",ylab="T. chumash (24_0163 H1)",cex.lab=1.4)
axis(2,at=seq(0,10,length.out=11)/11,dove_sc,las=2)
axis(1,at=seq(0,11,length.out=11)/11,tchum_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes, chromosome numbers here based on knulli
chtab<-matrix(c(1,43,1, ## 1 vs 4 vs 9 (dropped) arbitrary
	2,1392,7,
	4,43,2,
	5,56,5,
	6,1469,6,
	7,1510,10,
	8,113,4,
	10,1213,8,
	11,48,9,
	12,1403,11,
	13,1308,3),nrow=11,ncol=3,byrow=TRUE)



pdf("AlnPlot_Tchum_ed24_0163h1_24_0163h2.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:11){
	tdove<-grep(x=dfdat[,14],pattern=chtab[i,3]) 
	tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
	cc<-tdove[tdove %in% tch]
	subd<-dfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. chumash (24_0163 H2)",ylab="T. chumash (24_0163 H1)")
	title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17],col="cadetblue")
		}
	}
	if(chtab[i,1]==8){
		abline(h=c(56676246,60791657),lty=2) ## bounds super approximte
	}

}
dev.off()

i<-7
pdf("AlnPlot_Tchum_ed24_0163h1_0163h2-zoom.pdf",width=5,height=5)
par(mfrow=c(1,1))
par(mar=c(4.5,5.5,2.5,1.5))
tdove<-grep(x=dfdat[,14],pattern=chtab[i,3])
tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
cc<-tdove[tdove %in% tch]
subd<-dfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])

plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(54000000,69000000),ylim=c(54000000,69000000),cex.lab=1.4,xlab="T. chumash (24_0163 H2)",ylab="T. chumash (24_0163 H1)")
title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
N<-dim(subd)[1]
for(j in 1:N){
	if(subd[j,9]=="++"){
        	lines(subd[j,12:13],subd[j,16:17])
        }
        else{
        	lines(subd[j,12:13],yub-subd[j,16:17],col="cadetblue")
        }
}
abline(h=c(56676246,60791657),lty=2) ## bounds super approximte
dev.off()


## 24_0159 H1 (called dove below to keep code simple) vs 24_0159 H2
dat<-fread("out_synteny_E240159H1_E240159H2.psl",header=FALSE)
dfdat<-as.data.frame(dat)
## identify large scaffolds for each
xx<-table(dfdat[,14]) ## 24_0159 H1
dove<-names(xx)[xx>500] ## 10 
xx<-table(dfdat[,10]) ## 24_0159 H2
chum<-names(xx)[xx>500] ## 11
keep<-(dfdat[,14] %in% dove) & (dfdat[,10] %in% chum)
subDfdat<-dfdat[keep,]


tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tchum_sc<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,22,2)],pattern="Chr",replacement=""))
dove_sc<-as.numeric(gsub(x=unlist(strsplit(x=colnames(tab),split="_"))[seq(2,22,2)],pattern="Chr",replacement=""))

## normalize with respect to dove 

ntab<-tab
for(i in 1:11){
	ntab[,i]<-ntab[,i]/sum(ntab[,i])
}


pdf("SynTchum_ed24_0159h1_24_0159h2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. chumash (24_0159 H2)",ylab="T. chumash (24_0159 H1)",cex.lab=1.4)
axis(2,at=seq(0,10,length.out=11)/11,dove_sc,las=2)
axis(1,at=seq(0,11,length.out=11)/11,tchum_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes, chromosome numbers here based on knulli
chtab<-matrix(c(1,43,1, ## 1 vs 4 vs 9 (dropped) arbitrary
	2,1392,7,
	4,43,2,
	5,56,5,
	6,1469,6,
	7,1510,10,
	8,113,4,
	10,1213,8,
	11,48,9,
	12,1403,11,
	13,1308,3),nrow=11,ncol=3,byrow=TRUE)



pdf("AlnPlot_Tchum_ed24_0159h1_24_0159h2.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:11){
	tdove<-grep(x=dfdat[,14],pattern=chtab[i,3]) 
	tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
	cc<-tdove[tdove %in% tch]
	subd<-dfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. chumash (24_0159 H2)",ylab="T. chumash (24_0159 H1)")
	title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17],col="cadetblue")
		}
	}
	if(chtab[i,1]==8){
		abline(h=c(56676246,60791657),lty=2) ## bounds super approximte
	}

}
dev.off()

i<-7
pdf("AlnPlot_Tchum_ed24_0159h1_0159h2-zoom.pdf",width=5,height=5)
par(mfrow=c(1,1))
par(mar=c(4.5,5.5,2.5,1.5))
tdove<-grep(x=dfdat[,14],pattern=chtab[i,3])
tch<-grep(x=dfdat[,10],pattern=chtab[i,3])
cc<-tdove[tdove %in% tch]
subd<-dfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])

#plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(8000000,153000000),ylim=c(7000000,172000000),cex.lab=1.4,xlab="T. chumash (24_0159 H2)",ylab="T. chumash (24_0159 H1)")
plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(35000000,50000000),ylim=c(57000000,72000000),cex.lab=1.4,xlab="T. chumash (24_0159 H2)",ylab="T. chumash (24_0159 H1)")
title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
N<-dim(subd)[1]
for(j in 1:N){
	if(subd[j,9]=="++"){
        	lines(subd[j,12:13],subd[j,16:17])
        }
        else{
        	lines(subd[j,12:13],yub-subd[j,16:17],col="cadetblue")
        }
}
abline(h=c(62000000,68000000),lty=2) ## bounds super approximte
dev.off()
