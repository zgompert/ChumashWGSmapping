## summarize SV calls from SyRI
library(data.table)

## red genome

## now try for all genomes
library(scales)
ff<-c("syri_24_0159_hap1_24_0159_hap2syri.out","syri_24_0159_hap1_24_0160_hap1syri.out","syri_24_0159_hap1_24_0160_hap2syri.out","syri_24_0159_hap1_24_0161_hap1syri.out","syri_24_0159_hap1_24_0161_hap2syri.out","syri_24_0159_hap1_24_0162_hap1syri.out","syri_24_0159_hap1_24_0162_hap2syri.out","syri_24_0159_hap1_24_0163_hap1syri.out","syri_24_0159_hap1_24_0163_hap2syri.out")

annos<-c("HDR","INV","TRNAS","INVTR","DUP","INVDP","CPG","CPL","INS","TDM","DEL")
cs<-c("azure2","black","cadetblue","darkblue","goldenrod","gold","darkslategray","darkslategray","darkorchid","darkred","deeppink2")

y1<-seq(0.025,0.825,length.out=9)
y2<-y1+.025
y3<-(y1+y2)/2

pdf("fig_svannot_full_red.pdf",width=8,height=8)
par(mar=c(4.5,6,.5,.5))
lb<-60000000;ub<-64000000
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

ff<-c("syri_24_0159_hap1_24_0159_hap2syri.out","syri_24_0159_hap1_24_0160_hap1syri.out","syri_24_0159_hap1_24_0160_hap2syri.out","syri_24_0159_hap1_24_0161_hap1syri.out","syri_24_0159_hap1_24_0161_hap2syri.out","syri_24_0159_hap1_24_0162_hap1syri.out","syri_24_0159_hap1_24_0162_hap2syri.out","syri_24_0159_hap1_24_0163_hap1syri.out","syri_24_0159_hap1_24_0163_hap2syri.out")

nms<-c("24_0159H2","24_0160H1","24_0160H2","24_0161H1","24_0161H2","24_0162H1","24_0162H2","24_0163H1","24_0163H2")
annos<-c("INV","TRNAS","INVTR","DUP","INVDP","CPG","CPL","INS","TDM","DEL")

lb<-60000000;ub<-64000000
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
#  CPL   DUP   INS   INV INVDP INVTR 
#    2    11     1    12     4     2 
table(combined_ofull[,2])
#  CPG   CPL   DEL   DUP   INS   INV INVDP INVTR 
#    3     4     1    14     2    13     8     6 
dim(combined_ofull) ## 47 out of 51 unique
#[1] 76  5
dim(combined_obig) ## 29 out of 32 distinct
#[1] 32  5

xtable(combined_obig)
write.table(combined_ofull,file="syri_sv_red_full.csv",sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE)

