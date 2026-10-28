## make relative fitness for different LD figures 
library(data.table)
library(scales)
library(RColorBrewer)

## graphical parameters
cl<-1.4;ca<-1.1;cm<-1.5
##################################################################



### green sort #####
pdf("sfig_LD_minmax_green_sort.pdf",width=10,height=9)
par(mfrow=c(3,2))

par(mar=c(4.5,5.5,2.5,1.5))
##########

xmax<-3.5
load("predictCh8_sort_green_0.05.rdat")
cat(length(hipip),"\n")## 31
plot(density(relf_ran_AC,adj=1.2),xlim=c(.6,xmax),ylim=c(0,4),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
pp1<-round(mean(relf_ran_AC > 1),2)
pp2<-round(mean(relf_max_AC > 1),2)
legend(1.3,3.8,c(paste("Obs. > low, ",pp1,sep=""),paste("High > obs., ",pp2,sep="")),
       fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(a) Effect of LD on A/C, PIP 0.05",cex.main=cm)

plot(density(relf_ran_MM,adj=1.2),xlim=c(.6,xmax),ylim=c(0,6),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
pp1<-round(mean(relf_ran_MM > 1),2)
pp2<-round(mean(relf_max_MM > 1),2)
legend(1.3,5.7,c(paste("Obs. > low, ",pp1,sep=""),paste("High > obs., ",pp2,sep="")),
       fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(b) Effect of LD on MM, PIP 0.05",cex.main=cm)

load("predictCh8_sort_green_0.075.rdat")
cat(length(hipip),"\n") ## 21
plot(density(relf_ran_AC,adj=1.2),xlim=c(.6,xmax),ylim=c(0,4),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
pp1<-round(mean(relf_ran_AC > 1),2)
pp2<-round(mean(relf_max_AC > 1),2)
legend(1.3,3.8,c(paste("Obs. > low, ",pp1,sep=""),paste("High > obs., ",pp2,sep="")),
       fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(c) Effect of LD on A/C, PIP 0.075",cex.main=cm)

plot(density(relf_ran_MM,adj=1.2),xlim=c(.6,xmax),ylim=c(0,6),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
pp1<-round(mean(relf_ran_MM > 1),2)
pp2<-round(mean(relf_max_MM > 1),2)
legend(1.3,5.7,c(paste("Obs. > low, ",pp1,sep=""),paste("High > obs., ",pp2,sep="")),
       fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(d) Effect of LD on MM, PIP 0.075",cex.main=cm)

load("predictCh8_sort_green_0.1.rdat")
cat(length(hipip),"\n") ## 14
plot(density(relf_ran_AC,adj=1.2),xlim=c(.6,xmax),ylim=c(0,4),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
pp1<-round(mean(relf_ran_AC > 1),2)
pp2<-round(mean(relf_max_AC > 1),2)
legend(1.3,3.8,c(paste("Obs. > low, ",pp1,sep=""),paste("High > obs., ",pp2,sep="")),
       fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(e) Effect of LD on A/C, PIP 0.1",cex.main=cm)

plot(density(relf_ran_MM,adj=1.2),xlim=c(.6,xmax),ylim=c(0,6),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
pp1<-round(mean(relf_ran_MM > 1),2)
pp2<-round(mean(relf_max_MM > 1),2)
legend(1.3,5.7,c(paste("Obs. > low, ",pp1,sep=""),paste("High > obs., ",pp2,sep="")),
       fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(f) Effect of LD on MM, PIP 0.1",cex.main=cm)
dev.off()

### red sort #####
pdf("sfig_LD_minmax_red_sort.pdf",width=10,height=9)
par(mfrow=c(3,2))

par(mar=c(4.5,5.5,2.5,1.5))
##########

xmax<-2.5
load("predictCh8_sort_red_0.05.rdat")
cat(length(hipip),"\n") ## 26
plot(density(relf_ran_AC,adj=1.2),xlim=c(.6,xmax),ylim=c(0,10),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
pp1<-round(mean(relf_ran_AC > 1),2)
pp2<-round(mean(relf_max_AC > 1),2)
legend(1.3,9.5,c(paste("Obs. > low, ",pp1,sep=""),paste("High > obs., ",pp2,sep="")),
       fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(a) Effect of LD on A/C, PIP 0.05",cex.main=cm)

plot(density(relf_ran_MM,adj=1.2),xlim=c(.6,xmax),ylim=c(0,10),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
pp1<-round(mean(relf_ran_MM > 1),2)
pp2<-round(mean(relf_max_MM > 1),2)
legend(1.3,9.5,c(paste("Obs. > low, ",pp1,sep=""),paste("High > obs., ",pp2,sep="")),
       fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(b) Effect of LD on MM, PIP 0.05",cex.main=cm)

load("predictCh8_sort_red_0.075.rdat") 
cat(length(hipip),"\n") ## 17
plot(density(relf_ran_AC,adj=1.2),xlim=c(.6,xmax),ylim=c(0,15),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
pp1<-round(mean(relf_ran_AC > 1),2)
pp2<-round(mean(relf_max_AC > 1),2)
legend(1.3,14.25,c(paste("Obs. > low, ",pp1,sep=""),paste("High > obs., ",pp2,sep="")),
       fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(c) Effect of LD on A/C, PIP 0.075",cex.main=cm)

plot(density(relf_ran_MM,adj=1.2),xlim=c(.6,xmax),ylim=c(0,15),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
pp1<-round(mean(relf_ran_MM > 1),2)
pp2<-round(mean(relf_max_MM > 1),2)
legend(1.3,14.25,c(paste("Obs. > low, ",pp1,sep=""),paste("High > obs., ",pp2,sep="")),
       fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(d) Effect of LD on MM, PIP 0.075",cex.main=cm)

load("predictCh8_sort_red_0.1.rdat")
cat(length(hipip),"\n") ## 9
plot(density(relf_ran_AC,adj=1.2),xlim=c(.6,xmax),ylim=c(0,25),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_AC,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
pp1<-round(mean(relf_ran_AC > 1),2)
pp2<-round(mean(relf_max_AC > 1),2)
legend(1.3,23.75,c(paste("Obs. > low, ",pp1,sep=""),paste("High > obs., ",pp2,sep="")),
       fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(e) Effect of LD on A/C, PIP 0.1",cex.main=cm)

plot(density(relf_ran_MM,adj=1.2),xlim=c(.6,xmax),ylim=c(0,25),xlab="Relative mean fitness",ylab="Density",main="",type='n',cex.lab=cl,cex.axis=ca)
dd<-density(relf_ran_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#2166AC",.5),border=NA)
dd<-density(relf_max_MM,adj=1.2)
polygon(c(dd$x,rev(dd$x)),c(rep(0,length(dd$y)),rev(dd$y)),col=alpha("#B2182B",.5),border=NA)
abline(v=1,lty=3,lwd=2)
pp1<-round(mean(relf_ran_MM > 1),2)
pp2<-round(mean(relf_max_MM > 1),2)
legend(1.3,23.75,c(paste("Obs. > low, ",pp1,sep=""),paste("High > obs., ",pp2,sep="")),
       fill=c(alpha("#2166AC",.5),alpha("#B2182B",.5)),bty='n',cex=1.5)
title(main="(f) Effect of LD on MM, PIP 0.1",cex.main=cm)
dev.off()
