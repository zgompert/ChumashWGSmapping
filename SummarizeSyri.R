## summarize SV calls from SyRI
library(data.table)

sv<-fread("syri_24_0159_hap1_24_0163_hap1syri.out",header=FALSE)

table(sv[,11])
#V11
#    CPG     CPL     DEL     DUP   DUPAL     HDR     INS     INV   INVAL   INVDP 
#     63      40      22      73      98    4796      11      11      32      56 
#INVDPAL   INVTR INVTRAL   NOTAL     SYN   SYNAL     TDM   TRANS TRANSAL 
#     68      52      70     513     100    4998       9      50      67 

## inversions
inv<-which(sv[,11]=="INV")

cbind(sv[inv,2],as.numeric(unlist(sv[inv,3]))-as.numeric(unlist(sv[inv,2]))
#          V2     V2
#      <char>  <num>
# 1:  8537248    306
# 2: 11415890    921
# 3: 45398808    855
# 4: 46415456   1164
# 5: 62470076  10871
# 6: 62880396 724534
# 7: 64961996    404
# 8: 67014532    563
# 9: 71993690   1231
#10: 79585655   8561
#11: 95225149 541537

## translocations
trans<-which(sv[,11]=="TRANS")

cbind(sv[trans,2],as.numeric(unlist(sv[trans,3]))-as.numeric(unlist(sv[trans,2]))) ##50 of them


## summarize in focal region
lb<-60000000;ub<-65000000
lpos<-as.numeric(unlist(sv[,2]));upos<-as.numeric(unlist(sv[,3]))
creg<-which((lpos > lb & sv[,3] < upos) | (lpos > lb & lpos < ub) | (upos > lb & upos < ub) | (lpos < lb & upos > ub))

## annotations of interest
annos<-c("INV","TRNAS","INVTR","DUP","INVDP","CPG","CPL","INS","TDM","DEL","HDR")
cs<-c("black","cadetblue","darkblue","goldenrod","gold","darkslategray","darkslategray","darkorchid","darkred","deeppink2","darkolivegreen3")
keep<-which(unlist(sv[creg,11]) %in% annos)
kcreg<-creg[keep]

pdf("fig_svannot.pdf",width=8,height=3)
par(mar=c(4.5,4.5,.5,.5))
y<-rep(.15,length(annos));y[length(annos)]<-.4
plot(c(lb,ub),c(0,.5),type="n",xlab="Position (bp)",ylab="Variant annotation",axes=FALSE)
for(k in 1:length(annos)){
	a<-which(sv[kcreg,11]==annos[k])
	if(length(a) > 0){
		for(x in a){
			lines(c(lpos[kcreg][x],upos[kcreg][x]),c(y[k],y[k]),col=cs[k],lwd=8)
		}
	}
}
axis(1)
dev.off()

## now try for all genomes
library(scales)
ff<-c("syri_24_0163_hap1_24_0159_hap1syri.out","syri_24_0163_hap1_24_0159_hap2syri.out","syri_24_0163_hap1_24_0160_hap1syri.out","syri_24_0163_hap1_24_0160_hap2syri.out","syri_24_0163_hap1_24_0161_hap1syri.out","syri_24_0163_hap1_24_0161_hap2syri.out","syri_24_0163_hap1_24_0162_hap1syri.out","syri_24_0163_hap1_24_0162_hap2syri.out","syri_24_0163_hap1_24_0163_hap2syri.out")

annos<-c("HDR","INV","TRNAS","INVTR","DUP","INVDP","CPG","CPL","INS","TDM","DEL")
cs<-c("azure2","black","cadetblue","darkblue","goldenrod","gold","darkslategray","darkslategray","darkorchid","darkred","deeppink2")

y1<-seq(0.025,0.825,length.out=9)
y2<-y1+.025
y3<-(y1+y2)/2

pdf("fig_svannot_full.pdf",width=8,height=8)
par(mar=c(4.5,6,.5,.5))
lb<-60000000;ub<-65000000
plot(c(lb,ub),c(0,1),type="n",xlab="Position (bp)",ylab="",axes=FALSE)
for(fx in 1:length(ff)){
	sv<-fread(ff[fx],header=FALSE)
	lpos<-as.numeric(unlist(sv[,2]));upos<-as.numeric(unlist(sv[,3]))
	creg<-which((lpos > lb & sv[,3] < upos) | (lpos > lb & lpos < ub) | (upos > lb & upos < ub) | (lpos < lb & upos > ub))
	keep<-which(unlist(sv[creg,11]) %in% annos)
	kcreg<-creg[keep]
	y<-rep(y3[fx],length(annos));y[length(annos)]<-y3[fx]

	for(k in 1:length(annos)){
		a<-which(sv[kcreg,11]==annos[k])
		if(length(a) > 0){
			for(x in a){
				lines(c(lpos[kcreg][x],upos[kcreg][x]),c(y[k],y[k]),col=cs[k],lwd=8)
			}
		}
	}
}
axis(1)
axis(2,at=(y1+y2)/2,c("159_h1","159_h2","160_h1","160_h2","161_h1","161_h2","162_h1","162_h2","163_h2"),las=2,tick=FALSE,lwd=0)
legend(6.05e7,.98,legend=annos,fill=cs,ncol=4,bty='n')

dev.off()

pdf("fig_svannot_wide.pdf",width=8,height=8)
par(mar=c(4.5,6,.5,.5))
lb<-55000000;ub<-70000000
plot(c(lb,ub),c(0,1),type="n",xlab="Position (bp)",ylab="",axes=FALSE)
for(fx in 1:length(ff)){
	sv<-fread(ff[fx],header=FALSE)
	lpos<-as.numeric(unlist(sv[,2]));upos<-as.numeric(unlist(sv[,3]))
	creg<-which((lpos > lb & sv[,3] < upos) | (lpos > lb & lpos < ub) | (upos > lb & upos < ub) | (lpos < lb & upos > ub))
	keep<-which(unlist(sv[creg,11]) %in% annos)
	kcreg<-creg[keep]
	y<-rep(y3[fx],length(annos));y[length(annos)]<-y3[fx]

	for(k in 1:length(annos)){
		a<-which(sv[kcreg,11]==annos[k])
		if(length(a) > 0){
			for(x in a){
				lines(c(lpos[kcreg][x],upos[kcreg][x]),c(y[k],y[k]),col=cs[k],lwd=8)
			}
		}
	}
}
axis(1)
axis(2,at=(y1+y2)/2,c("159_h1","159_h2","160_h1","160_h2","161_h1","161_h2","162_h1","162_h2","163_h2"),las=2,tick=FALSE,lwd=0)
legend(6.05e7,.98,legend=annos,fill=cs,ncol=4,bty='n')

dev.off()

pdf("fig_svannot_verywide.pdf",width=8,height=8)
par(mar=c(4.5,6,.5,.5))
lb<-30000000;ub<-95000000
plot(c(lb,ub),c(0,1),type="n",xlab="Position (bp)",ylab="",axes=FALSE)
for(fx in 1:length(ff)){
	sv<-fread(ff[fx],header=FALSE)
	lpos<-as.numeric(unlist(sv[,2]));upos<-as.numeric(unlist(sv[,3]))
	creg<-which((lpos > lb & sv[,3] < upos) | (lpos > lb & lpos < ub) | (upos > lb & upos < ub) | (lpos < lb & upos > ub))
	keep<-which(unlist(sv[creg,11]) %in% annos)
	kcreg<-creg[keep]
	y<-rep(y3[fx],length(annos));y[length(annos)]<-y3[fx]

	for(k in 1:length(annos)){
		a<-which(sv[kcreg,11]==annos[k])
		if(length(a) > 0){
			for(x in a){
				lines(c(lpos[kcreg][x],upos[kcreg][x]),c(y[k],y[k]),col=cs[k],lwd=5)
			}
		}
	}
}
axis(1)
axis(2,at=(y1+y2)/2,c("159_h1","159_h2","160_h1","160_h2","161_h1","161_h2","162_h1","162_h2","163_h2"),las=2,tick=FALSE,lwd=0)
legend(6.05e7,.98,legend=annos,fill=cs,ncol=4,bty='n')

dev.off()


### hdr profiles
pdf("fig_hdr_green.pdf",width=10,heigh=5)
par(mar=c(5,5,1,1))
sv<-fread(ff[fx],header=FALSE)
library(RColorBrewer)
lch<-brewer.pal(n=9,"PuBu")
uu<-max(as.numeric(unlist(sv[,3])),na.rm=TRUE)
plot(c(0,1),c(0,1),xlim=c(0,uu),ylim=c(0,1),type="n",xlab="Position",ylab="% HDR",cex.lab=1.3)
ux<-seq(1,uu,by=100000)
lx<-length(ux)
xx<-(ux[-c(1:10)]+ux[1:(lx-10)])/2
for(fx in 1:length(ff)){
	sv<-fread(ff[fx],header=FALSE)
	uu<-max(as.numeric(unlist(sv[,3])),na.rm=TRUE)
	ux<-seq(1,uu,by=100000)
	lx<-length(ux)
	lpos<-as.numeric(unlist(sv[,2]));upos<-as.numeric(unlist(sv[,3]))
	hdr_win<-rep(0,(lx-10))
	for(i in 1:(lx-10)){
		hdr<-which(lpos > ux[i] & upos < ux[i+10] & sv[,11]=="HDR")
		if(length(hdr) > 0){
			hdr_win[i]<-0+sum(upos[hdr]-lpos[hdr])
		}
	}
	hdr_win<-hdr_win/1e6
	lines(xx,hdr_win,col=lch[fx])
}
dev.off()

#notal
annos<-c("NOTAL","INV","TRNAS","INVTR","DUP","INVDP","CPG","CPL","INS","TDM","DEL")
cs<-c("azure2","black","cadetblue","darkblue","goldenrod","gold","darkslategray","darkslategray","darkorchid","darkred","deeppink2")

y1<-seq(0.025,0.825,length.out=9)
y2<-y1+.025
y3<-(y1+y2)/2

pdf("fig_svannot_notal.pdf",width=8,height=8)
par(mar=c(4.5,6,.5,.5))
lb<-60000000;ub<-65000000
plot(c(lb,ub),c(0,1),type="n",xlab="Position (bp)",ylab="",axes=FALSE)
for(fx in 1:length(ff)){
	sv<-fread(ff[fx],header=FALSE)
	lpos<-as.numeric(unlist(sv[,2]));upos<-as.numeric(unlist(sv[,3]))
	creg<-which((lpos > lb & sv[,3] < upos) | (lpos > lb & lpos < ub) | (upos > lb & upos < ub) | (lpos < lb & upos > ub))
	keep<-which(unlist(sv[creg,11]) %in% annos)
	kcreg<-creg[keep]
	y<-rep(y3[fx],length(annos));y[length(annos)]<-y3[fx]

	for(k in 1:length(annos)){
		a<-which(sv[kcreg,11]==annos[k])
		if(length(a) > 0){
			for(x in a){
				lines(c(lpos[kcreg][x],upos[kcreg][x]),c(y[k],y[k]),col=cs[k],lwd=8)
			}
		}
	}
}
axis(1)
axis(2,at=(y1+y2)/2,c("159_h1","159_h2","160_h1","160_h2","161_h1","161_h2","162_h1","162_h2","163_h2"),las=2,tick=FALSE,lwd=0)
legend(6.05e7,.98,legend=annos,fill=cs,ncol=4,bty='n')

dev.off()

### test summaries/table
library(xtable)

ff<-c("syri_24_0163_hap1_24_0159_hap1syri.out","syri_24_0163_hap1_24_0159_hap2syri.out","syri_24_0163_hap1_24_0160_hap1syri.out","syri_24_0163_hap1_24_0160_hap2syri.out","syri_24_0163_hap1_24_0161_hap1syri.out","syri_24_0163_hap1_24_0161_hap2syri.out","syri_24_0163_hap1_24_0162_hap1syri.out","syri_24_0163_hap1_24_0162_hap2syri.out","syri_24_0163_hap1_24_0163_hap2syri.out")

nms<-c("24_0159H1","24_0159H2","24_0160H1","24_0160H2","24_0161H1","24_0161H2","24_0162H1","24_0162H2","24_0163H2")
annos<-c("INV","TRNAS","INVTR","DUP","INVDP","CPG","CPL","INS","TDM","DEL")

lb<-60000000;ub<-65000000
ofull<-vector("list",length(ff))
obig<-ofull
minsz<-1000
for(fx in 1:length(ff)){
	sv<-fread(ff[fx],header=FALSE)
	lpos<-as.numeric(unlist(sv[,2]));upos<-as.numeric(unlist(sv[,3]))
	creg<-which((lpos > lb & sv[,3] < upos) | (lpos > lb & lpos < ub) | (upos > lb & upos < ub) | (lpos < lb & upos > ub))
	keep<-which(unlist(sv[creg,11]) %in% annos)
	kcreg<-creg[keep]
	ofull[[fx]]<-cbind(rep(nms[fx],length(keep)),sv[kcreg,c(11,2,3)],as.numeric(unlist(sv[kcreg,3]))-as.numeric(unlist(sv[kcreg,2])))
	## fix size for ins
	ins<-which(ofull[[fx]][,2]=="INS")
	if(length(ins) > 0){
		for(a in ins){
			isz<-as.numeric(unlist(sv[creg[keep[a]],8]))-
				as.numeric(unlist(sv[creg[keep[a]],7]))
			ofull[[fx]][a,5]<-isz
		}
	}
	big<-which(ofull[[fx]][,5] >= minsz)
	obig[[fx]]<-ofull[[fx]][big,]
}

combined_obig<-do.call(rbind,obig)
combined_ofull<-do.call(rbind,ofull)
table(combined_obig[,2])
#  DUP   INV INVDP INVTR
#    5    14     4     8 
table(combined_ofull[,2])
#  CPL   DEL   DUP   INS   INV INVDP INVTR
#    2     3     6     1    19    30    15
dim(combined_ofull) ## 53 out of 76 unique
#[1] 76  5
dim(combined_obig) ## 28 out of 32 distinct
#[1] 32  5

xtable(combined_obig)
write.table(combined_ofull,file="syri_sv_green_full.csv",sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE)
