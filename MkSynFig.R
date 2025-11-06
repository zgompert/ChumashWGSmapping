## chumash edinburgh
library(data.table)
library(scales)

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


load(file="alnsTchum.rdat")

## plots
panelPlot<-function(dat=adat,chtab=chtab,bndtab1=NA,bndtab2=NA,flipY=FALSE,flipX=FALSE,scale=FALSE){
	dfdat<-as.data.frame(dat)
	## identify large scaffolds for each
	xx<-table(dfdat[,14]) 
	dove<-names(xx)[xx>500]  
	xx<-table(dfdat[,10]) 
	chum<-names(xx)[xx>500] 
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
	if(flipX==TRUE){
		bndtab2<-xub-as.matrix(bndtab2)
		bndtab2<-bndtab2[,c(2,1)]
		subd[,12]<-xub-as.numeric(subd[,12])
		subd[,13]<-xub-as.numeric(subd[,13])
	}
	if(flipY==TRUE){
		bndtab1<-yub-as.matrix(bndtab1)
		bndtab1<-bndtab1[,c(2,1)]
		subd[,16]<-yub-as.numeric(subd[,16])
		subd[,17]<-yub-as.numeric(subd[,17])
	}
	xlb<-min(bndtab2)
	xub<-max(bndtab2)
	mid<-(xlb+xub)/2
	xlb<-mid-3e6
	xub<-mid+3e6
	ylb<-min(bndtab1)
	yub<-max(bndtab1)
	mid<-(ylb+yub)/2
	ylb<-mid-3e6
	yub<-mid+3e6
	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',axes=FALSE,xlim=c(xlb,xub),ylim=c(ylb,yub),cex.lab=1.4,xlab="",ylab="")
	if(scale==TRUE){
		lines(c(xlb,xlb+1e6),c(yub-3e5,yub-3e5),lwd=1.5,col="firebrick")
		text(xlb+.5e6,yub-1e6,"1 Mb")
	}
	polygon(y=c(bndtab1[1,],rev(bndtab1[1,])),x=c(rep(bndtab2[1,],each=2)),col=alpha("aquamarine3",.2),border=NA)
	polygon(y=c(bndtab1[2,],rev(bndtab1[2,])),x=c(rep(bndtab2[2,],each=2)),col=alpha("gold3",.2),border=NA)
	polygon(y=c(bndtab1[3,],rev(bndtab1[3,])),x=c(rep(bndtab2[3,],each=2)),col=alpha("coral3",.2),border=NA)
	xub<-max(unlist(subd[,12:13]));yub<-max(unlist(subd[,16:17]))
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
	        	lines(subd[j,12:13],subd[j,16:17])
       		}
        	else{
        		lines(subd[j,12:13],yub-subd[j,16:17])
        	}
	}
	box()
}

chtab<-c(8,4,4)
ct<-1.2
pdf("exfig_dots.pdf",width=11,height=11)
par(mfrow=c(5,5))
par(mar=c(.7,.7,.7,.7))

## row 1
plot(0:1,0:1,type='n',xlab="",ylab="",axes=FALSE)
text(0.5,0.5,"24_0159 H1",cex=ct)
adat<-fread("out_synteny_E240159H1_E240159H2.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=FALSE,flipY=TRUE,scale=TRUE,
                  bndtab1=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1),bndtab2=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2))
adat<-fread("out_synteny_E240159H1_E240161H1.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=FALSE,flipY=TRUE,
                  bndtab1=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1),bndtab2=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1))
adat<-fread("out_synteny_E240159H1_E240163H1.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,flipY=TRUE,
                  bndtab1=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1),bndtab2=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1))
adat<-fread("out_synteny_E240159H1_E240163H2.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,flipY=TRUE,
                  bndtab1=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1),bndtab2=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2))
## row 2
adat<-fread("out_synteny_E240159H1_E240159H2.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=FALSE,flipY=TRUE,
                  bndtab1=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1),bndtab2=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2))
plot(0:1,0:1,type='n',xlab="",ylab="",axes=FALSE)
text(0.5,0.5,"24_0159 H1",cex=ct)
adat<-fread("out_synteny_E240159H2_E240161H1.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=FALSE,flipY=FALSE,
                  bndtab1=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2),bndtab2=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1))
adat<-fread("out_synteny_E240159H2_E240163H1.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,flipY=FALSE,
                  bndtab1=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2),bndtab2=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1))
adat<-fread("out_synteny_E240159H2_E240163H2.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,flipY=FALSE,
                  bndtab1=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2),bndtab2=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2))

## row 3 
adat<-fread("out_synteny_E240159H1_E240161H1.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipY=TRUE,
                  bndtab1=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1),bndtab2=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1))
adat<-fread("out_synteny_E240159H2_E240161H1.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=FALSE,flipY=FALSE,
                  bndtab1=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2),bndtab2=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1))
plot(0:1,0:1,type='n',xlab="",ylab="",axes=FALSE)
text(0.5,0.5,"24_0161 H1",cex=ct)
adat<-fread("out_synteny_E240161H1_E240163H1.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,
                  bndtab1=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1),bndtab2=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1))
adat<-fread("out_synteny_E240161H1_E240163H2.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,
                  bndtab1=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1),bndtab2=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2))

## row 4
adat<-fread("out_synteny_E240159H1_E240163H1.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,flipY=TRUE,
                  bndtab1=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1),bndtab2=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1))
adat<-fread("out_synteny_E240159H2_E240163H1.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,flipY=FALSE,
                  bndtab1=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2),bndtab2=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1))
adat<-fread("out_synteny_E240161H1_E240163H1.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,
                  bndtab1=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1),bndtab2=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1))
plot(0:1,0:1,type='n',xlab="",ylab="",axes=FALSE)
text(0.5,0.5,"24_0163 H1",cex=ct)
adat<-fread("out_synteny_E240163H1_E240163H2.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,flipY=TRUE,
                  bndtab1=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1),bndtab2=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2))

## row 5
adat<-fread("out_synteny_E240159H1_E240163H2.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,flipY=TRUE,
                  bndtab1=rbind(r1_24_0159h1,r2_24_0159h1,r3_24_0159h1),bndtab2=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2))
adat<-fread("out_synteny_E240159H2_E240163H2.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,flipY=FALSE,
                  bndtab1=rbind(r1_24_0159h2,r2_24_0159h2,r3_24_0159h2),bndtab2=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2))
adat<-fread("out_synteny_E240161H1_E240163H2.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,
                  bndtab1=rbind(r1_24_0161h1,r2_24_0161h1,r3_24_0161h1),bndtab2=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2))
adat<-fread("out_synteny_E240163H1_E240163H2.psl",header=FALSE)
panelPlot(dat=adat,chtab=chtab,flipX=TRUE,flipY=TRUE,
                  bndtab1=rbind(r1_24_0163h1,r2_24_0163h1,r3_24_0163h1),bndtab2=rbind(r1_24_0163h2,r2_24_0163h2,r3_24_0163h2))
plot(0:1,0:1,type='n',xlab="",ylab="",axes=FALSE)
text(0.5,0.5,"24_0163 H2",cex=ct)
dev.off()

