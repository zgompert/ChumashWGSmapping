## synteny plots for thinking about color, focused on the phased 
## cristinae genomes, the old melanic genome and the chumash genome

library(data.table)
library(scales)

## read synteny dat
dat_gsh1_mel<-fread("out_synteny_TcrGSH1_TcrM.psl",header=FALSE)
dat_gush2_mel<-fread("out_synteny_TcrGUSH2_TcrM.psl",header=FALSE)
dat_gsh1_chum<-fread("out_synteny_TcrGSH1_Tchum.psl",header=FALSE)
dat_gush2_chum<-fread("out_synteny_TcrGUSH2_Tchum.psl",header=FALSE)

################# Tcris GSH1 vs melanic #####
dfdat<-as.data.frame(dat_gsh1_mel)

## target = GSH1
## query = melanic

## get Chromosome 8 and scaffolds 2963, 702.1, 128, and 1845 in melanic
## ch8 in GSH1 = Scaffold_11__2_contigs__length_91821751
## mel-stripe, chr 8
ms<-c("ScRvNkF_2963","ScRvNkF_702.1","ScRvNkF_128","ScRvNkF_1845")
keep<-dfdat[,14]=="Scaffold_11__2_contigs__length_91821751" & dfdat[,10] %in% ms
subDfdat<-dfdat[keep,]
ch8sz<-91821751

pdf("AlnPlotColorGSH1xMelnic_Ch8.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:4){
	cc<-which(subDfdat[,10]==ms[i])
	subd<-subDfdat[cc,]
	if(i==1){## flip 2963
		subd[,12]<-17133644-subd[,12]
		subd[,13]<-17133644-subd[,13]
	}
	if(i==2){## flip 702.1
		subd[,12]<-14171514-subd[,12]
		subd[,13]<-14171514-subd[,13]
	}
	yub<-max(subd[,13]);xub<-max(subd[,17])	

	plot(y=as.numeric(subd[1,12:13]),x=as.numeric(subd[1,16:17]),type='n',xlim=c(0,ch8sz),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (GSH1)",ylab="T. cristinae (M)")
	title(main=ms[i],cex.main=1.4)
	#polygon(c(1096929,12418814,12418814,1096929),c(0,0,yub*1.1,yub*1.1),border=NA,col=alpha("gray",.5))
	#polygon(c(45229604,84397122,84397122,45229604),c(0,0,yub*1.1,yub*1.1),border=NA,col=alpha("gray",.5))
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(y=subd[j,12:13],x=subd[j,16:17])
		}
		else{
			#lines(y=subd[j,12:13],x=subd[j,16:17],col="cadetblue")
			lines(y=subd[j,12:13],x=ch8sz-subd[j,16:17],col="cadetblue")
		}
	}
	if(i==2){
		abline(h=c(14171514,14171514-4139489),col="red")
	}
	if(i==3){
		abline(h=c(1,6000000),col="red")
	}
}
dev.off()

################# Tcris GUSH2 vs melanic #####
dfdat<-as.data.frame(dat_gush2_mel)

## target = GUSH2
## query = melanic

## get Chromosome 8 and scaffolds 2963, 702.1, 128, and 1845 in melanic
## ch8 in GUSH2 = ScrX45T_23_HRSCAF_264
## mel-stripe, chr 8
ms<-c("ScRvNkF_2963","ScRvNkF_702.1","ScRvNkF_128","ScRvNkF_1845")
keep<-dfdat[,14]=="ScrX45T_23_HRSCAF_264" & dfdat[,10] %in% ms
subDfdat<-dfdat[keep,]
ch8sz<-max(as.numeric(dfdat[dfdat[,14]=="ScrX45T_23_HRSCAF_264",17])) ## approximate

pdf("AlnPlotColorGUSH2xMelnic_Ch8.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:4){
	cc<-which(subDfdat[,10]==ms[i])
	subd<-subDfdat[cc,]
	if(i==1){## flip 2963
		subd[,12]<-17133644-subd[,12]
		subd[,13]<-17133644-subd[,13]
	}
	if(i==2){## flip 702.1
		subd[,12]<-14171514-subd[,12]
		subd[,13]<-14171514-subd[,13]
	}
	yub<-max(subd[,13]);xub<-max(subd[,17])	

	plot(y=as.numeric(subd[1,12:13]),x=as.numeric(subd[1,16:17]),type='n',xlim=c(0,ch8sz),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (GUSH2)",ylab="T. cristinae (M)")
	title(main=ms[i],cex.main=1.4)
	#polygon(c(1096929,12418814,12418814,1096929),c(0,0,yub*1.1,yub*1.1),border=NA,col=alpha("gray",.5))
	#polygon(c(45229604,84397122,84397122,45229604),c(0,0,yub*1.1,yub*1.1),border=NA,col=alpha("gray",.5))
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(y=subd[j,12:13],x=subd[j,16:17])
		}
		else{
			#lines(y=subd[j,12:13],x=subd[j,16:17],col="cadetblue")
			lines(y=subd[j,12:13],x=ch8sz-subd[j,16:17],col="cadetblue")
		}
	}
	if(i==2){
		abline(h=c(14171514,14171514-4139489),col="red")
	}
	if(i==3){
		abline(h=c(1,6000000),col="red")
	}
}
dev.off()

################# combo plot #####
dfdat<-as.data.frame(dat_gsh1_mel)

## target = GSH1
## query = melanic

## get Chromosome 8 and scaffolds 2963, 702.1, 128, and 1845 in melanic
## ch8 in GSH1 = Scaffold_11__2_contigs__length_91821751
## mel-stripe, chr 8
ms<-c("ScRvNkF_2963","ScRvNkF_702.1","ScRvNkF_128","ScRvNkF_1845")
keep<-dfdat[,14]=="Scaffold_11__2_contigs__length_91821751" & dfdat[,10] %in% ms
subDfdat1<-dfdat[keep,]
ch8sz1<-91821751

dfdat<-as.data.frame(dat_gush2_mel)

## target = GUSH2
## query = melanic

## get Chromosome 8 and scaffolds 2963, 702.1, 128, and 1845 in melanic
## ch8 in GUSH2 = ScrX45T_23_HRSCAF_264
## mel-stripe, chr 8
ms<-c("ScRvNkF_2963","ScRvNkF_702.1","ScRvNkF_128","ScRvNkF_1845")
keep<-dfdat[,14]=="ScrX45T_23_HRSCAF_264" & dfdat[,10] %in% ms
subDfdat2<-dfdat[keep,]
ch8sz2<-max(as.numeric(dfdat[dfdat[,14]=="ScrX45T_23_HRSCAF_264",17])) ## approximate

cs<-c("brown","cadetblue","orange","violet")
pdf("AlnPlotColorTcrComboxMelnic_Ch8.pdf",width=10,height=16)
par(mfrow=c(8,1))
par(mar=c(4.5,5.5,2.5,1.5))
ch8sz<-ch8sz1;subDfdat<-subDfdat1
#plot(y=1,x=1,type='n',xlim=c(0,ch8sz),ylim=c(0,1.8e7),cex.lab=1.4,xlab="T. cristinae (GSH1)",ylab="T. cristinae (M)")
for(i in rev(1:4)){
	cc<-which(subDfdat[,10]==ms[i])
	subd<-subDfdat[cc,]
	if(i==1){## flip 2963
		subd[,12]<-17133644-subd[,12]
		subd[,13]<-17133644-subd[,13]
	}
	if(i==2){## flip 702.1
		subd[,12]<-14171514-subd[,12]
		subd[,13]<-14171514-subd[,13]
	}
	yub<-max(subd[,13]);xub<-max(subd[,17])	
	plot(y=1,x=1,type='n',xlim=c(0,ch8sz),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (GSH1)",ylab="T. cristinae (M)")
	N<-dim(subd)[1]
	if(i==2){
		polygon(c(-1,15e7,15e7,-1),c(14171514,14171514,14171514-4139489,14171514-4139489),border=NA,col=alpha("gray",.5))
	}
	if(i==3){
		polygon(c(-1,15e7,15e7,-1),c(1,1,6000000,6000000),border=NA,col=alpha("gray",.5))
	}
	#polygon(c(45229604,84397122,84397122,45229604),c(0,0,yub*1.1,yub*1.1),border=NA,col=alpha("gray",.5))
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(y=subd[j,12:13],x=subd[j,16:17],col=cs[i])
		}
		else{
			lines(y=subd[j,12:13],x=ch8sz-subd[j,16:17],col=cs[i])
		}
	}
}
ch8sz<-ch8sz2;subDfdat<-subDfdat2
#plot(y=1,x=1,type='n',xlim=c(0,ch8sz),ylim=c(0,1.8e7),cex.lab=1.4,xlab="T. cristinae (GUSH2)",ylab="T. cristinae (M)")
for(i in rev(1:4)){
	cc<-which(subDfdat[,10]==ms[i])
	subd<-subDfdat[cc,]
	if(i==1){## flip 2963
		subd[,12]<-17133644-subd[,12]
		subd[,13]<-17133644-subd[,13]
	}
	if(i==2){## flip 702.1
		subd[,12]<-14171514-subd[,12]
		subd[,13]<-14171514-subd[,13]
	}
	yub<-max(subd[,13]);xub<-max(subd[,17])	
	plot(y=1,x=1,type='n',xlim=c(0,ch8sz),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (GUSH2)",ylab="T. cristinae (M)")
	N<-dim(subd)[1]
	if(i==2){
		polygon(c(-1,15e7,15e7,-1),c(14171514,14171514,14171514-4139489,14171514-4139489),border=NA,col=alpha("gray",.5))
	}
	if(i==3){
		polygon(c(-1,15e7,15e7,-1),c(1,1,6000000,6000000),border=NA,col=alpha("gray",.5))
	}
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(y=subd[j,12:13],x=subd[j,16:17],col=cs[i])
		}
		else{
			lines(y=subd[j,12:13],x=ch8sz-subd[j,16:17],col=cs[i])
		}
	}
}
dev.off()

################# combo plot with chumash #####
dfdat<-as.data.frame(dat_gsh1_chum)

## target = GSH1
## query = chum

##  chr 8
keep<-dfdat[,14]=="Scaffold_11__2_contigs__length_91821751" & dfdat[,10] == "ScN4ago_113_HRSCAF_504" 
subDfdat1<-dfdat[keep,]
ch8sz1<-91821751

dfdat<-as.data.frame(dat_gush2_chum)

## target = GUSH2
## query = chum

## chr 8
keep<-dfdat[,14]=="ScrX45T_23_HRSCAF_264" & dfdat[,10] == "ScN4ago_113_HRSCAF_504" 
subDfdat2<-dfdat[keep,]
ch8sz2<-max(as.numeric(dfdat[dfdat[,14]=="ScrX45T_23_HRSCAF_264",17])) ## approximate

pdf("AlnPlotColorTcrxChum_Ch8.pdf",width=10,height=5)
par(mfrow=c(1,2))
par(mar=c(4.5,5.5,2.5,1.5))
ch8sz<-ch8sz1;subDfdat<-subDfdat1
#plot(y=1,x=1,type='n',xlim=c(0,ch8sz),ylim=c(0,1.8e7),cex.lab=1.4,xlab="T. cristinae (GSH1)",ylab="T. cristinae (M)")
subd<-subDfdat
yub<-max(subd[,13]);xub<-max(subd[,17])	
plot(y=1,x=1,type='n',xlim=c(0,ch8sz),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (GSH1)",ylab="T. chumash")
## approximate gwa signal in chumash
polygon(c(-1,2e8,2e8,-1),c(5.5e7,5.5e7,7e7,7e7),border=NA,col=alpha("gray",.5))
polygon(c(2.5e7,4.5e7,4.5e7,2.5e7),c(-1,-1,2e8,2e8),border=NA,col=alpha("tan",.5))
N<-dim(subd)[1]
for(j in 1:N){
	if(subd[j,9]=="++"){
		lines(y=subd[j,12:13],x=subd[j,16:17])
	}
	else{
		lines(y=subd[j,12:13],x=ch8sz-subd[j,16:17])
	}
}

ch8sz<-ch8sz2;subDfdat<-subDfdat2
subd<-subDfdat
yub<-max(subd[,13]);xub<-max(subd[,17])	
plot(y=1,x=1,type='n',xlim=c(0,ch8sz),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (GUSH2)",ylab="T. chumash")
N<-dim(subd)[1]
polygon(c(-1,2e8,2e8,-1),c(5.5e7,5.5e7,7e7,7e7),border=NA,col=alpha("gray",.5))
polygon(c(2.5e7,4.5e7,4.5e7,2.5e7),c(-1,-1,2e8,2e8),border=NA,col=alpha("tan",.5))
for(j in 1:N){
	if(subd[j,9]=="++"){
		lines(y=subd[j,12:13],x=subd[j,16:17])
	}
	else{
		lines(y=subd[j,12:13],x=ch8sz-subd[j,16:17])
	}
}
dev.off()
