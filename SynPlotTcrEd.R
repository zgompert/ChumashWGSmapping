library(data.table)

#240016H1 = VP green from C

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, gsh2
chtab<-matrix(c(1,8483,12,13,13,
	2,14640,6,5,6,
	3,42935,2,3,2,
	4,42912,1,1,1,
	5,18722,7,12,12,
	6,9928,8,4,5,
	7,10660,10,10,8,
	8,7748,11,11,4,
	9,16151,5,8,9,
	10,14160,4,7,7,
	11,12033,9,9,10,
	12,12380,13,6,11,
	13,14101,3,2,3),nrow=13,ncol=5,byrow=TRUE)

## read synteny dat
#########################################################################
## gsh1 x 240016H1 
##############
dat_gsh1_vpgush1<-fread("out_synteny_TcrGSH1_TcrE240016H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsh1_vpgush1)

## target = GSH1
## query = VPGUSH1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
vpgush1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gshxx<-table(dfdat[,10])
gsh1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% vpgush1Ch) & (dfdat[,10] %in% gsh1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gsh1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
tc_vpgush1<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for tcr


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSH1_VPGUSH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (VP GUS1)",ylab="T. cristinae (H GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsh1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_vpgush1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, vpgush1
chtab<-matrix(c(1,8483,12,13,9,
	2,14640,6,5,6,
	3,42935,2,3,3,
	4,42912,1,1,2,
	5,18722,7,12,5,
	6,9928,8,4,8,
	7,10660,10,10,13,
	8,7748,11,11,1,
	9,16151,5,8,10,
	10,14160,4,7,12,
	11,12033,9,9,11,
	12,12380,13,6,7,
	13,14101,3,2,4),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GSH1_VPGUSH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gsh1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep="")) 
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_vpgush1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_gsh1[tcr_gsh1 %in% tcr_vpgush1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="VP GUS1",ylab="H GS1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsh1==chtab[i,4])]-subd[j,16:17])
		}
	}
}
dev.off()

#### comparing two new haplotypes
## read synteny dat
#########################################################################
## 240016H1 x 240016H2
##############
dat_vp_gush1_gush2<-fread("out_synteny_TcrE240016H1_TcrE240016H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_vp_gush1_gush2)

## target = VPGUSH1
## query = VPGUSH2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gush2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gshxx<-table(dfdat[,10])
gush1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gush2Ch) & (dfdat[,10] %in% gush1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_vp_gush1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_vp_gush2<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement=""))

#sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for tcr


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_VP_GUSH1_GUSH2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (VP GUS2)",ylab="T. cristinae (VP GUS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_vp_gush1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_vp_gush2,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, vpgush1=vpgush2
chtab<-matrix(c(1,8483,12,13,9,
	2,14640,6,5,6,
	3,42935,2,3,3,
	4,42912,1,1,2,
	5,18722,7,12,5,
	6,9928,8,4,8,
	7,10660,10,10,13,
	8,7748,11,11,1,
	9,16151,5,8,10,
	10,14160,4,7,12,
	11,12033,9,9,11,
	12,12380,13,6,7,
	13,14101,3,2,4),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_VP_GUSH1_GUSH2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_vp_gush2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
	tcr_vp_gush1<-grep(x=subDfdat[,14],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_vp_gush1[tcr_vp_gush1 %in% tcr_vp_gush2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="VP GUS2",ylab="VP GUS 1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

### we have a new melanic genome, looking at that now
## read synteny dat
#########################################################################
## gsh1 x 240038H1 
##############
dat_gsh1_vpmh1<-fread("out_synteny_TcrGSH1_TcrE240038H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsh1_vpmh1)

## target = GSH1
## query = VPMH1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
vpmh1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gshxx<-table(dfdat[,10])
gsh1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% vpmh1Ch) & (dfdat[,10] %in% gsh1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gsh1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
tc_vpmh1<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for tcr


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSH1_VPMH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (VP M1)",ylab="T. cristinae (H GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsh1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_vpmh1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, vpmh1
chtab<-matrix(c(1,8483,12,13,8,
	2,14640,6,5,5,
	3,42935,2,3,1,
	4,42912,1,1,2,
	5,18722,7,12,13,
	6,9928,8,4,6,
	7,10660,10,10,12,
	8,7748,11,11,4,
	9,16151,5,8,9,
	10,14160,4,7,10,
	11,12033,9,9,11,
	12,12380,13,6,7,
	13,14101,3,2,3),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GSH1_VPMH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        tcr_gsh1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_vpmh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
        cc<-tcr_gsh1[tcr_gsh1 %in% tcr_vpmh1]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="VP M1",ylab="H GS1")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],sizes[which(tc_gsh1==chtab[i,4])]-subd[j,16:17])
                }
        }
}
dev.off()

## zoomed version 2 to 4e7 on ch8
pdf("AlnPlotTcris_GSH1_VPMH1_zoom.pdf",width=4,height=4)
par(mar=c(4.5,4.5,1,1))
i<-8        
tcr_gsh1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep=""))
if(chtab[i,5] > 9){
	val<-paste(chtab[i,5],sep="")
}else{
	val<-paste("0",chtab[i,5],sep="")
}
tcr_vpmh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
cc<-tcr_gsh1[tcr_gsh1 %in% tcr_vpmh1]
subd<-subDfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])

plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(2.4e7,3e7),ylim=c(2.4e7,3e7),cex.lab=1.4,xlab="VP M1",ylab="H GS1")
N<-dim(subd)[1]
for(j in 2:N){
	if(subd[j,9]=="++"){
        	lines(subd[j,12:13],subd[j,16:17])
        }
        else{
        	lines(subd[j,12:13],sizes[which(tc_gsh1==chtab[i,4])]-subd[j,16:17])
        }
}
abline(a=0,b=1,lty=3,col="gray")
dev.off()

#########################################################################
## gush2 x 240038H1 
##############
dat_gush2_vpmh1<-fread("out_synteny_TcrGUSH2_TcrE240038H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gush2_vpmh1)

## target = GUSH2
## query = VPMH1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
vpmh1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gush2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% vpmh1Ch) & (dfdat[,10] %in% gush2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gush2<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,53,4)])
tc_vpmh1<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GUSH2_VPMH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (VP M1)",ylab="T. cristinae (H GUS2)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gush2,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_vpmh1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gush2, vpmh1
chtab<-matrix(c(1,8483,12,15,8,
	2,14640,6,1,5,
	3,42935,2,3,1,
	4,42912,1,35,2,
	5,18722,7,10,13,
	6,9928,8,44,6,
	7,10660,10,7,12,
	8,7748,11,23,4,
	9,16151,5,21,9,
	10,14160,4,16,10,
	11,12033,9,12,11,
	12,12380,13,36,7,
	13,14101,3,8,3),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GUSH2_VPMH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        tcr_gush2<-grep(x=subDfdat[,14],pattern=paste("45T_",chtab[i,4],"_",sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_vpmh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
        cc<-tcr_gush2[tcr_gush2 %in% tcr_vpmh1]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="VP M1",ylab="H GUS2")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()

pdf("AlnPlotTcris_GUSH2_VPMH1_zoom.pdf",width=4,height=4)
par(mar=c(4.5,4.5,1,1))
i<-8        
tcr_gush2<-grep(x=subDfdat[,14],pattern=paste("45T_",chtab[i,4],"_",sep=""))
if(chtab[i,5] > 9){
	val<-paste(chtab[i,5],sep="")
}else{
	val<-paste("0",chtab[i,5],sep="")
}
tcr_vpmh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
cc<-tcr_gush2[tcr_gush2 %in% tcr_vpmh1]
subd<-subDfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])

plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(2.4e7,5e7),ylim=c(2.4e7,5e7),cex.lab=1.4,xlab="VP M1",ylab="H GUS2")
N<-dim(subd)[1]
for(j in 2:N){
	if(subd[j,9]=="++"){
        	lines(subd[j,12:13],subd[j,16:17])
        }
        else{
        	lines(subd[j,12:13],yub-subd[j,16:17])
        }
}
abline(a=0,b=1,lty=3,col="gray")
abline(v=c(2.45e7,4.05e7),col="red",lty=2)
dev.off()

pdf("AlnPlotTcris_GUSH2_VPMH1_zoom2.pdf",width=4,height=4)
par(mar=c(4.5,4.5,1,1))
i<-8        
tcr_gush2<-grep(x=subDfdat[,14],pattern=paste("45T_",chtab[i,4],"_",sep=""))
if(chtab[i,5] > 9){
	val<-paste(chtab[i,5],sep="")
}else{
	val<-paste("0",chtab[i,5],sep="")
}
tcr_vpmh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
cc<-tcr_gush2[tcr_gush2 %in% tcr_vpmh1]
subd<-subDfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])

plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(2.4e7,3e7),ylim=c(3.3e7,3.9e7),cex.lab=1.4,xlab="VP M1",ylab="H GUS2")
N<-dim(subd)[1]
for(j in 2:N){
	if(subd[j,9]=="++"){
        	lines(subd[j,12:13],subd[j,16:17])
        }
        else{
        	lines(subd[j,12:13],yub-subd[j,16:17])
        }
}
abline(a=8.57e6,b=1,lty=3,col="gray")
dev.off()

#### comparing the two melanic haplotypes
## read synteny dat
#########################################################################
## 240038H1 x 240038H2
##############
dat_vp_mh1_mh2<-fread("out_synteny_TcrE240038H1_TcrE240038H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_vp_mh1_mh2)

## target = VPMH1
## query = VPMH2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
mh2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
mh1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% mh2Ch) & (dfdat[,10] %in% mh1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_vp_mh1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_vp_mh2<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement=""))

#sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for tcr


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_VP_MH1_MH2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (VP M2)",ylab="T. cristinae (VP M1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_vp_mh1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_vp_mh2,las=2)
box()
dev.off()

## chrom number, gs,gsr1, gush2, vpmh1=vpmh2
chtab<-matrix(c(1,8483,12,15,8,
	2,14640,6,1,5,
	3,42935,2,3,1,
	4,42912,1,35,2,
	5,18722,7,10,13,
	6,9928,8,44,6,
	7,10660,10,7,12,
	8,7748,11,23,4,
	9,16151,5,21,9,
	10,14160,4,16,10,
	11,12033,9,12,11,
	12,12380,13,36,7,
	13,14101,3,8,3),nrow=13,ncol=5,byrow=TRUE)




pdf("AlnPlotTcris_VP_MH1_MH2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_vp_mh2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
	tcr_vp_mh1<-grep(x=subDfdat[,14],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_vp_mh1[tcr_vp_mh1 %in% tcr_vp_mh2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="VP M2",ylab="VP M1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

### FHA COMPARISONS ######

library(data.table)

#240175H1 = FHA Stripe

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, gsh2
chtab<-matrix(c(1,8483,12,13,13,
	2,14640,6,5,6,
	3,42935,2,3,2,
	4,42912,1,1,1,
	5,18722,7,12,12,
	6,9928,8,4,5,
	7,10660,10,10,8,
	8,7748,11,11,4,
	9,16151,5,8,9,
	10,14160,4,7,7,
	11,12033,9,9,10,
	12,12380,13,6,11,
	13,14101,3,2,3),nrow=13,ncol=5,byrow=TRUE)

## read synteny data
#########################################################################
## gsh1 x 240175H1 
##############
dat_gsh1_fhagsh1<-fread("out_synteny_TcrGSH1_TcrE240175H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsh1_fhagsh1)

## target = GSH1
## query =  FHA GSH1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
fhagsh1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gsh1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% fhagsh1Ch) & (dfdat[,10] %in% gsh1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gsh1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
tc_fhagsh1<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for tcr


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSH1_FHA0175_GSH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (FHA 0175 GS1)",ylab="T. cristinae (H GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsh1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_fhagsh1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, fha1075gsh1
chtab<-matrix(c(1,8483,12,13,7,
	2,14640,6,5,5,
	3,42935,2,3,2,
	4,42912,1,1,1,
	5,18722,7,12,13,
	6,9928,8,4,6,
	7,10660,10,10,12,
	8,7748,11,11,4,
	9,16151,5,8,9,
	10,14160,4,7,10,
	11,12033,9,9,11,
	12,12380,13,6,8,
	13,14101,3,2,3),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GSH1_FHA0175_GSH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gsh1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep="")) 
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_fhagsh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_gsh1[tcr_gsh1 %in% tcr_fhagsh1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="FHA 0175 GS1",ylab="H GS1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsh1==chtab[i,4])]-subd[j,16:17])
		}
	}
}
dev.off()

####################
## gush2 x 240175H1 
##############
dat_gush2_fhagsh1<-fread("out_synteny_TcrGUSH2_TcrE240175H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gush2_fhagsh1)

## target = GUSH2
## query = FHAGSH1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
fhagsh1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gush2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% fhagsh1Ch) & (dfdat[,10] %in% gush2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gush2<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,53,4)])
tc_fhagsh1<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GUSH2_FHA0175_GSH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (FHA0175 GS1)",ylab="T. cristinae (H GUS2)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gush2,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_fhagsh1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gush2, fhagsh1
chtab<-matrix(c(1,8483,12,15,7,
	2,14640,6,1,5,
	3,42935,2,3,2,
	4,42912,1,35,1,
	5,18722,7,10,13,
	6,9928,8,44,6,
	7,10660,10,7,12,
	8,7748,11,23,4,
	9,16151,5,21,9,
	10,14160,4,16,10,
	11,12033,9,12,11,
	12,12380,13,36,8,
	13,14101,3,8,3),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GUSH2_FHA0175_GSH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        tcr_gush2<-grep(x=subDfdat[,14],pattern=paste("45T_",chtab[i,4],"_",sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_fhagsh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
        cc<-tcr_gush2[tcr_gush2 %in% tcr_fhagsh1]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="FHA0175 GS1",ylab="H GUS2")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()

###############################################

pdf("AlnPlotTcris_GUSH2_VPMH1_zoom.pdf",width=4,height=4)
par(mar=c(4.5,4.5,1,1))
i<-8        
tcr_gush2<-grep(x=subDfdat[,14],pattern=paste("45T_",chtab[i,4],"_",sep=""))
if(chtab[i,5] > 9){
	val<-paste(chtab[i,5],sep="")
}else{
	val<-paste("0",chtab[i,5],sep="")
}
tcr_vpmh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
cc<-tcr_gush2[tcr_gush2 %in% tcr_vpmh1]
subd<-subDfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])

plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(2.4e7,5e7),ylim=c(2.4e7,5e7),cex.lab=1.4,xlab="VP M1",ylab="H GUS2")
N<-dim(subd)[1]
for(j in 2:N){
	if(subd[j,9]=="++"){
        	lines(subd[j,12:13],subd[j,16:17])
        }
        else{
        	lines(subd[j,12:13],yub-subd[j,16:17])
        }
}
abline(a=0,b=1,lty=3,col="gray")
abline(v=c(2.45e7,4.05e7),col="red",lty=2)
dev.off()

pdf("AlnPlotTcris_GUSH2_VPMH1_zoom2.pdf",width=4,height=4)
par(mar=c(4.5,4.5,1,1))
i<-8        
tcr_gush2<-grep(x=subDfdat[,14],pattern=paste("45T_",chtab[i,4],"_",sep=""))
if(chtab[i,5] > 9){
	val<-paste(chtab[i,5],sep="")
}else{
	val<-paste("0",chtab[i,5],sep="")
}
tcr_vpmh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
cc<-tcr_gush2[tcr_gush2 %in% tcr_vpmh1]
subd<-subDfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])

plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(2.4e7,3e7),ylim=c(3.3e7,3.9e7),cex.lab=1.4,xlab="VP M1",ylab="H GUS2")
N<-dim(subd)[1]
for(j in 2:N){
	if(subd[j,9]=="++"){
        	lines(subd[j,12:13],subd[j,16:17])
        }
        else{
        	lines(subd[j,12:13],yub-subd[j,16:17])
        }
}
abline(a=8.57e6,b=1,lty=3,col="gray")
dev.off()

#### comparing two new FHA stripe haplotypes
## read synteny dat
#########################################################################
## 240175H1 x 240175H2
##############
dat_fha_gsh1_gsh2<-fread("out_synteny_TcrE240175H1_TcrE240175H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_fha_gsh1_gsh2)

## target = FHA 175 GSH1
## query = FHA 175 GSH2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsh2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gsh1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsh2Ch) & (dfdat[,10] %in% gsh1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_fha_gsh1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_fha_gsh2<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement=""))



## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_FHA0175_GSH1_GSH2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (FHA0175 GS2)",ylab="T. cristinae (FHA0175 GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_fha_gsh1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_fha_gsh2,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, fhagsh1=fhagsh2
chtab<-matrix(c(1,8483,12,13,7,
        2,14640,6,5,5,
        3,42935,2,3,2,
        4,42912,1,1,1,
        5,18722,7,12,13,
        6,9928,8,4,6,
        7,10660,10,10,12,
        8,7748,11,11,4,
        9,16151,5,8,9,
        10,14160,4,7,10,
        11,12033,9,9,11,
        12,12380,13,6,8,
        13,14101,3,2,3),nrow=13,ncol=5,byrow=TRUE)

pdf("AlnPlotTcris_FHA_GSH1_GSH2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_fha_gsh2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
	tcr_fha_gsh1<-grep(x=subDfdat[,14],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_fha_gsh1[tcr_fha_gsh1 %in% tcr_fha_gsh2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	
	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="FHA0175 GS2",ylab="FHA0175 GS 1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

################ out_synteny_E240175H1_TcrE240176H1.psl ###############
## 240175H1 x 240176H1
##############
dat_fha_gsh1_gsh1<-fread("out_synteny_E240175H1_TcrE240176H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_fha_gsh1_gsh1)

## target = FHA 175 GSH1
## query = FHA 176 GSH1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gshBCh<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gshACh<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gshBCh) & (dfdat[,10] %in% gshACh)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_fha_gshA<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_fha_gshB<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement="")) 



## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_FHA0175_GSH1_FHA0176_GSH2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (FHA0176 GS1)",ylab="T. cristinae (FHA0175 GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_fha_gshA,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_fha_gshB,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, fhagsh1 176, fhagsh1 175
chtab<-matrix(c(1,8483,12,8,7,
        2,14640,6,5,5,
        3,42935,2,2,2,
        4,42912,1,1,1,
        5,18722,7,13,13,
        6,9928,8,6,6,
        7,10660,10,12,12,
        8,7748,11,4,4,
        9,16151,5,10,9,
        10,14160,4,9,10,
        11,12033,9,11,11,
        12,12380,13,7,8,
        13,14101,3,3,3),nrow=13,ncol=5,byrow=TRUE)

pdf("AlnPlotTcris_FHA0175_GSH1_0176GSH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_fha_gshB<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
	tcr_fha_gshA<-grep(x=subDfdat[,14],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_fha_gshA[tcr_fha_gshA %in% tcr_fha_gshB]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="FHA0176 GS1",ylab="FHA0175 GS 1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

#out_synteny_TcrE240176H1_TcrE240176H2.ps

################ out_synteny_E240176H1_TcrE240176H2.psl ###############
## 240176H1 x 240176H2
##############
dat_fha_gsh1_gsh2<-fread("out_synteny_TcrE240176H1_TcrE240176H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_fha_gsh1_gsh2)

## target = FHA 176 GSH1
## query = FHA 176 GSH2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsh2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gsh1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsh2Ch) & (dfdat[,10] %in% gsh1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_fha_gsh1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_fha_gsh2<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement="")) 



## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_FHA0176_GSH1_GSH2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (FHA0176 GS2)",ylab="T. cristinae (FHA0176 GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_fha_gshA,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_fha_gshB,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, fhagsh1 176 h1 and h2, fhagsh1 175
chtab<-matrix(c(1,8483,12,8,7,
        2,14640,6,5,5,
        3,42935,2,2,2,
        4,42912,1,1,1,
        5,18722,7,13,13,
        6,9928,8,6,6,
        7,10660,10,12,12,
        8,7748,11,4,4,
        9,16151,5,10,9,
        10,14160,4,9,10,
        11,12033,9,11,11,
        12,12380,13,7,8,
        13,14101,3,3,3),nrow=13,ncol=5,byrow=TRUE)

pdf("AlnPlotTcris_FHA0176GSH1_0176GSH2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	if(chtab[i,4] > 9){
		val<-paste(chtab[i,4],sep="")
	}else{
		val<-paste("0",chtab[i,4],sep="")
	}
	tcr_fha_gsh2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
	tcr_fha_gsh1<-grep(x=subDfdat[,14],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_fha_gsh1[tcr_fha_gsh1 %in% tcr_fha_gsh2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="FHA0176 GS2",ylab="FHA0176 GS 2")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

## 176 is a het for the inverted translocation ##


## read synteny dat
#########################################################################
## gsh1 x 240087H1 
##############
dat_gsh1_vpgush1<-fread("out_synteny_TcrGSH1_TcrE240087H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsh1_vpgush1)

## target = GSH1
## query = VPGUSH1 = 0087

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
vpgush1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gshxx<-table(dfdat[,10])
gsh1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% vpgush1Ch) & (dfdat[,10] %in% gsh1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gsh1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
tc_vpgush1<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for tcr


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSH1_VP0087_GSH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (VP0087 GS1)",ylab="T. cristinae (H GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsh1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_vpgush1,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, vp0087 GSH1
chtab<-matrix(c(1,8483,12,13,5,
	2,14640,6,5,4,
	3,42935,2,3,2,
	4,42912,1,1,1,
	5,18722,7,12,12,
	6,9928,8,4,6,
	7,10660,10,10,11,
	8,7748,11,11,3,
	9,16151,5,8,9,
	10,14160,4,7,7,
	11,12033,9,9,10,
	12,12380,13,6,8,
	13,14101,3,2,13),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GSH1_VP0087_GSH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gsh1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep="")) 
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_vpgush1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_gsh1[tcr_gsh1 %in% tcr_vpgush1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="VP0087 GS1",ylab="H GS1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsh1==chtab[i,4])]-subd[j,16:17])
		}
	}
}
dev.off()


## read synteny dat
#########################################################################
## gsh1 x 240087H2 
##############
dat_gsh1_vpgush2<-fread("out_synteny_TcrGSH1_TcrE240087H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsh1_vpgush2)

## target = GSH1
## query = VPGUSH1 = 0087

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
vpgush2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gshxx<-table(dfdat[,10])
gsh1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% vpgush2Ch) & (dfdat[,10] %in% gsh1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gsh1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
tc_vpgush1<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement=""))

sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for tcr


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSH1_VP0087_GSH2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (VP0087 GS2)",ylab="T. cristinae (H GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsh1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_vpgush1,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, vp0087 GSH2
chtab<-matrix(c(1,8483,12,13,5,
	2,14640,6,5,4,
	3,42935,2,3,2,
	4,42912,1,1,1,
	5,18722,7,12,12,
	6,9928,8,4,6,
	7,10660,10,10,11,
	8,7748,11,11,3,
	9,16151,5,8,9,
	10,14160,4,7,7,
	11,12033,9,9,10,
	12,12380,13,6,8,
	13,14101,3,2,13),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GSH1_VP0087_GSH2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gsh1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep="")) 
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_vpgsh2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
	cc<-tcr_gsh1[tcr_gsh1 %in% tcr_vpgsh2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="VP0087 GS1",ylab="H GS1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsh1==chtab[i,4])]-subd[j,16:17])
		}
	}
}
dev.off()


####################
## gush2 x 240087H1 
##############
dat_gush2_vpgsh1<-fread("out_synteny_TcrGUSH2_TcrE240087H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gush2_vpgsh1)

## target = GUSH2
## query = VPGSH1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
vpgsh1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gush2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% vpgsh1Ch) & (dfdat[,10] %in% gush2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gush2<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,53,4)])
tc_vpgsh1<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GUSH2_FHA0087_GSH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (FHA0087 GS1)",ylab="T. cristinae (H GUS2)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gush2,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_vpgsh1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gush2, vpgsh1
chtab<-matrix(c(1,8483,12,15,5,
	2,14640,6,1,4,
	3,42935,2,3,2,
	4,42912,1,35,1,
	5,18722,7,10,12,
	6,9928,8,44,6,
	7,10660,10,7,12,
	8,7748,11,23,3,
	9,16151,5,21,9,
	10,14160,4,16,7,
	11,12033,9,12,10,
	12,12380,13,36,8,
	13,14101,3,8,13),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GUSH2_FHA0087_GSH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        tcr_gush2<-grep(x=subDfdat[,14],pattern=paste("45T_",chtab[i,4],"_",sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_vpgsh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
        cc<-tcr_gush2[tcr_gush2 %in% tcr_vpgsh1]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="VP0087 GS1",ylab="H GUS2")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()


pdf("AlnPlotTcris_GUSH2_VP0087_GSH1_zoom.pdf",width=4,height=4)
par(mar=c(4.5,4.5,1,1))
i<-8        
tcr_gush2<-grep(x=subDfdat[,14],pattern=paste("45T_",chtab[i,4],"_",sep=""))
if(chtab[i,5] > 9){
	val<-paste(chtab[i,5],sep="")
}else{
	val<-paste("0",chtab[i,5],sep="")
}
tcr_vpgsh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
cc<-tcr_gush2[tcr_gush2 %in% tcr_vpgsh1]
subd<-subDfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])

plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(2.4e7,5e7),ylim=c(2.4e7,5e7),cex.lab=1.4,xlab="VP GS1",ylab="H GUS2")
N<-dim(subd)[1]
for(j in 2:N){
	if(subd[j,9]=="++"){
        	lines(subd[j,12:13],subd[j,16:17])
        }
        else{
        	lines(subd[j,12:13],yub-subd[j,16:17])
        }
}
abline(a=0,b=1,lty=3,col="gray")
abline(v=c(2.45e7,4.05e7),col="red",lty=2)
dev.off()

pdf("AlnPlotTcris_GUSH2_VP0087_GSH1_zoom2.pdf",width=4,height=4)
par(mar=c(4.5,4.5,1,1))
i<-8        
tcr_gush2<-grep(x=subDfdat[,14],pattern=paste("45T_",chtab[i,4],"_",sep=""))
if(chtab[i,5] > 9){
	val<-paste(chtab[i,5],sep="")
}else{
	val<-paste("0",chtab[i,5],sep="")
}
tcr_vpgsh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
cc<-tcr_gush2[tcr_gush2 %in% tcr_vpgsh1]
subd<-subDfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])

plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(3.2e7,4e7),ylim=c(3.2e7,4e7),cex.lab=1.4,xlab="VP GS1",ylab="H GUS2")
N<-dim(subd)[1]
for(j in 2:N){
	if(subd[j,9]=="++"){
        	lines(subd[j,12:13],subd[j,16:17])
        }
        else{
        	lines(subd[j,12:13],yub-subd[j,16:17])
        }
}
abline(a=0,b=1,lty=3,col="gray")
dev.off()


#### comparing the two VP 0087 GS haplotypes
## read synteny dat
#########################################################################
## 240087H1 x 240087H2
##############
dat_vp_gsh1_gsh2<-fread("out_synteny_TcrE240087H1_TcrE240087H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_vp_gsh1_gsh2)

## target = GSH1 0087
## query = GSH2 0087

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsh2Ch<-names(xx)[xx>200]
xx<-table(dfdat[,10])
gsh1Ch<-names(xx)[xx>200]
keep<-(dfdat[,14] %in% gsh2Ch) & (dfdat[,10] %in% gsh1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_vp_gsh1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_vp_gsh2<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement=""))

#sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for tcr


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_VP0087_GSH1_GSH2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (VP0087 GS2)",ylab="T. cristinae (VP0087 GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_vp_gsh1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_vp_gsh2,las=2)
box()
dev.off()

## chrom number, gs,gsr1, gush2, vpgsh1=2
chtab<-matrix(c(1,8483,12,15,5,
	2,14640,6,1,4,
	3,42935,2,3,2,
	4,42912,1,35,1,
	5,18722,7,10,12,
	6,9928,8,44,6,
	7,10660,10,7,12,
	8,7748,11,23,3,
	9,16151,5,21,9,
	10,14160,4,16,7,
	11,12033,9,12,10,
	12,12380,13,36,8,
	13,14101,3,8,13),nrow=13,ncol=5,byrow=TRUE)



pdf("AlnPlotTcris_VP0087_GSH1_GSH2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_vp_gsh2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
	tcr_vp_gsh1<-grep(x=subDfdat[,14],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_vp_gsh1[tcr_vp_gsh1 %in% tcr_vp_gsh2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="VP0087 GS2",ylab="VP0087 GS1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

## zoom 
pdf("AlnPlotTcris_VP0087_GSH1_GSH2_zoom.pdf",width=11,height=11)
par(mar=c(4.5,4.5,1,1))
i<-8
if(chtab[i,5] > 9){
	val<-paste(chtab[i,5],sep="")
}else{
	val<-paste("0",chtab[i,5],sep="")
}
tcr_vp_gsh2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
tcr_vp_gsh1<-grep(x=subDfdat[,14],pattern=paste("hap1_Chr",val,sep=""))
cc<-tcr_vp_gsh1[tcr_vp_gsh1 %in% tcr_vp_gsh2]
subd<-subDfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(2.4e7,5e7),ylim=c(2.4e7,5e7),cex.lab=1.4,xlab="VP0087 GS2",ylab="VP0087 GS1")
N<-dim(subd)[1]
for(j in 2:N){
	if(subd[j,9]=="++"){
		lines(subd[j,12:13],subd[j,16:17])
	}
	else{
		lines(subd[j,12:13],yub-subd[j,16:17])
	}
}

dev.off()

################ out_synteny_TcrE240038H1_TcrE240087H1.psl ###############
## 240038H1 x 240087H1
##############
dat_vp_mh1_gsh1<-fread("out_synteny_TcrE240038H1_TcrE240087H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_vp_mh1_gsh1)

## target = VP 038 MH1
## query = VP 087 GSH1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsh1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
mh1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% mh1Ch) & (dfdat[,10] %in% gsh1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_vp_gsh1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_vp_mh1<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement="")) 



## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_VP0038_MH1_VP0087_GSH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (VP0087 GS1)",ylab="T. cristinae (VP0038 M1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_vp_mh1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_vp_gsh1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, vpm, vpgsh1
chtab<-matrix(c(1,8483,12,8,5,
	2,14640,6,5,4,
	3,42935,2,1,2,
	4,42912,1,2,1,
	5,18722,7,13,12,
	6,9928,8,6,6,
	7,10660,10,12,11,
	8,7748,11,4,3,
	9,16151,5,9,9,
	10,14160,4,10,7,
	11,12033,9,11,10,
	12,12380,13,7,8,
	13,14101,3,3,13),nrow=13,ncol=5,byrow=TRUE)

pdf("AlnPlotTcris_VP0038_MH1_VP0087_GSH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}

	tcr_vp_gsh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
	if(chtab[i,4] > 9){
		val<-paste(chtab[i,4],sep="")
	}else{
		val<-paste("0",chtab[i,4],sep="")
	}
	tcr_vp_mh1<-grep(x=subDfdat[,14],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_vp_mh1[tcr_vp_mh1 %in% tcr_vp_gsh1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="VP0038 M1",ylab="VP0087 GS 1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

##### new refugio alignments ###

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, gsh2
chtab<-matrix(c(1,8483,12,13,13,
	2,14640,6,5,6,
	3,42935,2,3,2,
	4,42912,1,1,1,
	5,18722,7,12,12,
	6,9928,8,4,5,
	7,10660,10,10,8,
	8,7748,11,11,4,
	9,16151,5,8,9,
	10,14160,4,7,7,
	11,12033,9,9,10,
	12,12380,13,6,11,
	13,14101,3,2,3),nrow=13,ncol=5,byrow=TRUE)

## read synteny dat
#########################################################################
## out_synteny_TcrGSR1_TcrE240072H1.psl = stripe x stripe
##############
dat_gsr1_gs072r1<-fread("out_synteny_TcrGSR1_TcrE240072H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsr1_gs072r1)

## target = GSR1
## query = R12 0072 GSR1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gs072r1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gsr1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gs072r1Ch) & (dfdat[,10] %in% gsr1Ch)
subDfdat<-dfdat[keep,]## retains all, 12 or 13 CH

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gsr1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
tc_gs072r1<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for tcr


## normalize  
ntab<-tab
for(i in 1:12){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSR1_R0072_GSR1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (R0072 GS1)",ylab="T. cristinae (R GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsr1,las=2)
axis(1,at=seq(0,12,length.out=12)/12,tc_gs072r1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, r0072 gsr1
chtab<-matrix(c(1,8483,12,13,11,
	2,14640,6,5,7,
	3,42935,2,3,1,
	4,42912,1,1,1,
	5,18722,7,12,6,
	6,9928,8,4,5,
	7,10660,10,10,12,
	8,7748,11,11,3,
	9,16151,5,8,8,
	10,14160,4,7,9,
	11,12033,9,9,10,
	12,12380,13,6,4,
	13,14101,3,2,2),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GSR1_R0072_GSR1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gsr1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,3],"_",sep="")) 
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_gs072r1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_gsr1[tcr_gsr1 %in% tcr_gs072r1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="R0072 GS1",ylab="R GS1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsr1==chtab[i,3])]-subd[j,16:17])
		}
	}
}
dev.off()

## zoom to look at a funny bit of 8

pdf("AlnPlotTcris_GSR1_R0072_GSR1_zoom.pdf",width=11,height=11)
par(mfrow=c(1,1))
par(mar=c(4.5,4.5,1,1))
i<-8
tcr_gsr1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,3],"_",sep="")) 
if(chtab[i,5] > 9){
	val<-paste(chtab[i,5],sep="")
}else{
	val<-paste("0",chtab[i,5],sep="")
}
tcr_gs072r1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
cc<-tcr_gsr1[tcr_gsr1 %in% tcr_gs072r1]
subd<-subDfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])	
	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(1.1e7,2.1e7),ylim=c(1.1e7,2.1e7),cex.lab=1.4,xlab="R0072 GS1",ylab="R GS1")
N<-dim(subd)[1]
for(j in 2:N){
	if(subd[j,9]=="++"){
		lines(subd[j,12:13],subd[j,16:17])
	}
	else{
		lines(subd[j,12:13],sizes[which(tc_gsr1==chtab[i,3])]-subd[j,16:17])
	}
}
dev.off()

#### comparing two new haplotypes from R striped
## read synteny dat
#########################################################################
## 240072H1 x 240072H2
##############
dat_r_gsr1_gsr2<-fread("out_synteny_TcrE240072H1_TcrE240072H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_r_gsr1_gsr2)

## target = r12 0072 GSR1
## query = r12 0072 GSR2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsr2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gsr1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsr2Ch) & (dfdat[,10] %in% gsr1Ch)
subDfdat<-dfdat[keep,]## retains all, 12 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_r_gsr1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_r_gsr2<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:12){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_R_GSR1_GSR2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (R0072 GS2)",ylab="T. cristinae (R0072 GS1)",cex.lab=1.4)
axis(2,at=seq(0,12,length.out=12)/12,tc_r_gsr1,las=2)
axis(1,at=seq(0,12,length.out=12)/12,tc_r_gsr2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, r0072 gsr1=2
chtab<-matrix(c(1,8483,12,13,11,
	2,14640,6,5,7,
	3,42935,2,3,1,
	4,42912,1,1,1,
	5,18722,7,12,6,
	6,9928,8,4,5,
	7,10660,10,10,12,
	8,7748,11,11,3,
	9,16151,5,8,8,
	10,14160,4,7,9,
	11,12033,9,9,10,
	12,12380,13,6,4,
	13,14101,3,2,2),nrow=13,ncol=5,byrow=TRUE)

pdf("AlnPlotTcris_R_GSR1_GSR2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_r_gsr2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
	tcr_r_gsr1<-grep(x=subDfdat[,14],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_r_gsr1[tcr_r_gsr1 %in% tcr_r_gsr2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="R0072 GS2",ylab="R0072 GS 1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

## 2 might be the melanic, it has the same inverted translocation I think from main mountain and maybe even lacks the big external inversion that we thought sepearted mountains, in other words maybe the between mountain inversion does not include melanics

## look at green and then come back to this

#########################################################################
## out_synteny_TcrGSR1_TcrE240073H1.psl = stripe x green
##############
dat_gsr1_gus073r1<-fread("out_synteny_TcrGSR1_TcrE240073H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsr1_gus073r1)

## target = GSR1
## query = R12 0073 GUSR1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gus073r1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gsr1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gus073r1Ch) & (dfdat[,10] %in% gsr1Ch)
subDfdat<-dfdat[keep,]## retains all, 12 or 13 CH 

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gsr1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
tc_gus073r1<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for tcr


## normalize  
ntab<-tab
for(i in 1:12){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSR1_R0073_GUSR1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (R0073 GUS1)",ylab="T. cristinae (R GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsr1,las=2)
axis(1,at=seq(0,12,length.out=12)/12,tc_gus073r1,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, r0073 gusr1
chtab<-matrix(c(1,8483,12,13,9,
	2,14640,6,5,6,
	3,42935,2,3,1,
	4,42912,1,1,1,
	5,18722,7,12,5,
	6,9928,8,4,4,
	7,10660,10,10,11,
	8,7748,11,11,3,
	9,16151,5,8,7,
	10,14160,4,7,8,
	11,12033,9,9,10,
	12,12380,13,6,12,
	13,14101,3,2,2),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GSR1_R0073_GUSR1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gsr1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,3],"_",sep="")) 
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_gus073r1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_gsr1[tcr_gsr1 %in% tcr_gus073r1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="R0073 GUS1",ylab="R GS1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsr1==chtab[i,3])]-subd[j,16:17])
		}
	}
}
dev.off()

#### comparing two new haplotypes from R green
## read synteny dat
#########################################################################
## 240073H1 x 240073H2
##############
dat_r_gusr1_gusr2<-fread("out_synteny_TcrE240073H1_TcrE240073H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_r_gusr1_gusr2)

## target = r12 0073 GUSR1
## query = r12 0073 GUSR2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gusr2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gusr1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gusr2Ch) & (dfdat[,10] %in% gusr1Ch)
subDfdat<-dfdat[keep,]## retains all, 12 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_r_gusr1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_r_gusr2<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:12){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_R0073_GUSR1_GUSR2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (R0073 GUS2)",ylab="T. cristinae (R0073 GUS1)",cex.lab=1.4)
axis(2,at=seq(0,12,length.out=12)/12,tc_r_gusr1,las=2)
axis(1,at=seq(0,12,length.out=12)/12,tc_r_gusr2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, r0073 gusr1 = 2
chtab<-matrix(c(1,8483,12,13,9,
	2,14640,6,5,6,
	3,42935,2,3,1,
	4,42912,1,1,1,
	5,18722,7,12,5,
	6,9928,8,4,4,
	7,10660,10,10,11,
	8,7748,11,11,3,
	9,16151,5,8,7,
	10,14160,4,7,8,
	11,12033,9,9,10,
	12,12380,13,6,12,
	13,14101,3,2,2),nrow=13,ncol=5,byrow=TRUE)

pdf("AlnPlotTcris_R0073_GUSR1_GUSR2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_r_gusr2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
	tcr_r_gusr1<-grep(x=subDfdat[,14],pattern=paste("hap1_Chr",val,sep=""))
	cc<-tcr_r_gusr1[tcr_r_gusr1 %in% tcr_r_gusr2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="R0073 GUS2",ylab="R0073 GUS 1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

### putative melanic allele from Refugio vs main mountain melanic
## AHH, have not done this

#### SEE
#out_synteny_TcrE240072H2_TcrE240038H1.psl
#out_synteny_TcrE240072H2_TcrE240039H1.psl
#out_synteny_TcrE240072H1_TcrE240038H1.psl
#out_synteny_TcrE240072H1_TcrE240039H1.psl

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, r0072 gsr1
chtab<-matrix(c(1,8483,12,13,11,
        2,14640,6,5,7,
        3,42935,2,3,1,
        4,42912,1,1,1,
        5,18722,7,12,6,
        6,9928,8,4,5,
        7,10660,10,10,12,
        8,7748,11,11,3,
        9,16151,5,8,8,
        10,14160,4,7,9,
        11,12033,9,9,10,
        12,12380,13,6,4,
        13,14101,3,2,2),nrow=13,ncol=5,byrow=TRUE)

#### comparing refguio putative mealnic to vp melanics 
## read synteny dat
#########################################################################
## 240072H2 x 240038H1
##############
dat_gsr2_mh1<-fread("out_synteny_TcrE240072H2_TcrE240038H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsr2_mh1)

## target = r12 0072 GSR2
## query = vp 0038 MH1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsr2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
mh1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsr2Ch) & (dfdat[,10] %in% mh1Ch)
subDfdat<-dfdat[keep,]## retains 12 or 13 chs

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_r_gsr2<-as.numeric(gsub(pattern="hap2_Chr",x=colnames(tab),replacement=""))
tc_vp_mh1<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_R0072_GSR2_VP0038MH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (VP0038 M1)",ylab="T. cristinae (R0072 GS2)",cex.lab=1.4)
axis(2,at=seq(0,12,length.out=12)/12,tc_r_gsr2,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_vp_mh1,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, vp0038 mh1, r0072 gsr1
chtab<-matrix(c(1,8483,12,10,11,
        2,14640,6,8,7,
        3,42935,2,1,1, ## order on which is 3 vs 4 for 38 currently arbitrary
        4,42912,1,2,1,
        5,18722,7,6,6,
        6,9928,8,5,5,
        7,10660,10,12,12,
        8,7748,11,4,3,
        9,16151,5,7,8,
        10,14160,4,9,9,
        11,12033,9,11,10,
        12,12380,13,13,4,
        13,14101,3,3,2),nrow=13,ncol=5,byrow=TRUE)

pdf("AlnPlotTcris_R0073_GSR2_VP0038_MH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	if(chtab[i,4] > 9){
		val<-paste(chtab[i,4],sep="")
	}else{
		val<-paste("0",chtab[i,4],sep="")
	}
	tcr_vp_mh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_r_gsr2<-grep(x=subDfdat[,14],pattern=paste("hap2_Chr",val,sep=""))
	cc<-tcr_r_gsr2[tcr_r_gsr2 %in% tcr_vp_mh1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="VP0038 M1",ylab="R0072 GS 2")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

## 240072H2 x 240039H1
##############
dat_gsr2_mh1<-fread("out_synteny_TcrE240072H2_TcrE240039H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsr2_mh1)

## target = r12 0072 GSR2
## query = vp 0039 MH1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsr2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
mh1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsr2Ch) & (dfdat[,10] %in% mh1Ch)
subDfdat<-dfdat[keep,]## retains 12 or 13 chs

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_r_gsr2<-as.numeric(gsub(pattern="hap2_Chr",x=colnames(tab),replacement=""))
tc_vp_mh1<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_R0072_GSR2_VP0039MH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (VP0039 M1)",ylab="T. cristinae (R0072 GS2)",cex.lab=1.4)
axis(2,at=seq(0,12,length.out=12)/12,tc_r_gsr2,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_vp_mh1,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, vp0039 mh1, r0072 gsr1
chtab<-matrix(c(1,8483,12,5,11,
        2,14640,6,7,7,
        3,42935,2,1,1, ## order on which is 3 vs 4 for 39 currently arbitrary
        4,42912,1,3,1,
        5,18722,7,13,6,
        6,9928,8,4,5,
        7,10660,10,10,12,
        8,7748,11,12,3,
        9,16151,5,6,8,
        10,14160,4,8,9,
        11,12033,9,9,10,
        12,12380,13,11,4,
        13,14101,3,2,2),nrow=13,ncol=5,byrow=TRUE)

pdf("AlnPlotTcris_R0073_GSR2_VP0039_MH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	if(chtab[i,4] > 9){
		val<-paste(chtab[i,4],sep="")
	}else{
		val<-paste("0",chtab[i,4],sep="")
	}
	tcr_vp_mh1<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
	if(chtab[i,5] > 9){
		val<-paste(chtab[i,5],sep="")
	}else{
		val<-paste("0",chtab[i,5],sep="")
	}
	tcr_r_gsr2<-grep(x=subDfdat[,14],pattern=paste("hap2_Chr",val,sep=""))
	cc<-tcr_r_gsr2[tcr_r_gsr2 %in% tcr_vp_mh1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="VP0039 M1",ylab="R0072 GS 2")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()
