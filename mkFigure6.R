## make figure 6, min max LD 
library(data.table)
library(scales)
library(RColorBrewer)

## graphical parameters
cl<-1.4;ca<-1.1;cm<-1.5

load("../gemma/LDchange.rdat")
##################################################################

### figure 6 manuscript #####
pdf("fig_LD_minmax.pdf",width=10.5,height=7)
layout(matrix(c(1,2,2,3,4,4,5,5),nrow=2,ncol=4,byrow=TRUE),widths=c(3,1.5,1.5,3),heights=c(3,3))
par(mar=c(4.5,5.5,2.5,1.5))
##########

lx<-min(c(rgRan[,2],rgMax[,2],pscoreRG))
ux<-max(c(rgRan[,2],rgMax[,2],pscoreRG))
ly<-min(c(gbRan[,2],gbMax[,2],pscoreGB))
uy<-max(c(gbRan[,2],gbMax[,2],pscoreGB))

plot(pscoreRG,pscoreGB,pch=19,cex.lab=cl,cex.axis=ca,col=alpha("black",.2),xlab="RG score",ylab="GB score",xlim=c(lx,ux),ylim=c(ly,uy))
points(rgRan[,2],gbRan[,2],pch=21,bg=alpha("#2166AC",.7))
title(main="(a) Polygenic scores low LD",cex.main=cm)

plot(pscoreRG,pscoreGB,pch=19,cex.lab=cl,cex.axis=ca,col=alpha("black",.7),xlab="RG score",ylab="GB score",xlim=c(lx,ux),ylim=c(ly,uy))
title(main="(b) Observed polygenic scores",cex.main=cm)

plot(pscoreRG,pscoreGB,pch=19,cex.lab=cl,cex.axis=ca,col=alpha("black",.2),xlab="RG score",ylab="GB score",xlim=c(lx,ux),ylim=c(ly,uy))
points(rgMax[,2],gbMax[,2],pch=21,bg=alpha("#B2182B",.7))
title(main="(c) Polygenic scores high LD",cex.main=cm)

plot(density(relf_ran_AC,adj=1.2),xlim=c(.6,2.25),ylim=c(0,5),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
legend(1.3,4.75,c("Obs. > low, pp 0.94","High > obs., pp 0.95"),fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(d) Effect of LD on A/C",cex.main=cm)

plot(density(relf_ran_MM,adj=1.2),xlim=c(.6,2.25),ylim=c(0,5),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
legend(1.3,4.75,c("Obs. > low, pp 0.81","High > obs., pp 0.82"),fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(d) Effect of LD on MM",cex.main=cm)

dev.off()
