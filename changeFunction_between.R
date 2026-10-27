## function to test for an increase in LD (here measured as PC correlations)
## beyond null expectations from random mortality or linear selection
## regions of interest for PCA are supplied by the user

## here only between region correlations are considered,
## this thus ignores correlations between PC1 and PC2 for the same region

## regions = matrix nrows = number of regions, ncols = 2 lower and upper bounds
## G = matrix of genotypes, nrows = SNPs, ncols = individuals
## pos = vector, position of SNPs, to cross reference with regions, same order as in G
## nreps = number of null sims
## trt = data frame with treatment information, e.g., sort_green_phfull_tchum2019.txt
LDChangePC_between<-function(regions=NA,G=NA,pos=NA,nreps=1000,
		     trt=NA){
	
	## PCS from regions
	nreg<-dim(regions)[1]
	regsnps<-vector("list",nreg) ## snps in regions
	pco<-vector("list",nreg) ## pca results by region
	pcx<-vector("list",nreg) ## first 2 PCs
	for(k in 1:nreg){
		regsnps[[k]]<-which(pos >= regions[k,1] & pos <= regions[k,2])
		pco[[k]]<-prcomp(t(G[regsnps[[k]],]),center=TRUE,scale=FALSE)
		pcx[[k]]<-pco[[k]]$x[,1:2]
	}
	pcs<-do.call(cbind,pcx) ## pcs across regions

	## get within region combinations
	ttt<-rep(LETTERS[1:nreg],each=2)
	nc<-nreg*2
	ttm<-matrix(NA,nrow=nr,ncol=nc)
	for(i in 1:nc){for(j in 1:nc){if(i!=j){ttm[i,j]<-as.numeric(ttt[i]==ttt[j])}}}
	## 1 = within region, 0 = between region

	## organize the relevant information from trt, including survival
	## same order as G
	t_ac<-which(trt$Treatment=="AC")
	t_mm<-which(trt$Treatment=="MM")

	t_ac_surv<-which(trt$Treatment=="AC" & trt$Survival==1)
	t_mm_surv<-which(trt$Treatment=="MM" & trt$Survival==1)

	## calculate observed correlations
	cors_ac<-abs(cor(pcs[t_ac,]))
	cors_mm<-abs(cor(pcs[t_mm,]))
	cors_sac<-abs(cor(pcs[t_ac_surv,]))
	cors_smm<-abs(cor(pcs[t_mm_surv,]))

	obs_ac<-mean(cors_ac[upper.tri(cors_ac)][ttm[upper.tri(ttm)]==0])
	obs_mm<-mean(cors_mm[upper.tri(cors_mm)][ttm[upper.tri(ttm)]==0])
	obs_acs<-mean(cors_sac[upper.tri(cors_sac)][ttm[upper.tri(ttm)]==0])
	obs_mms<-mean(cors_smm[upper.tri(cors_smm)][ttm[upper.tri(ttm)]==0])

	## relative to random mortality
	sran<-trt$Survival
	null_acs<-rep(NA,nreps)
	for(i in 1:nreps){
	        sran[t_ac]<-sample(sran[t_ac],length(t_ac),replace=FALSE)
	        r_ac_surv<-which(trt$Treatment=="AC" & sran==1)
	        rcors_acs<-abs(cor(pcs[r_ac_surv,]))
	        rcorsbar_acs<-mean(rcors_acs[upper.tri(rcors_acs)][ttm[upper.tri(ttm)]==0])
	        null_acs[i]<-rcorsbar_acs
	}
	pv_ran_ac<-(1+sum(null_acs >= obs_acs))/(nreps+1)
	ts_ac<-obs_acs-obs_ac
	null_ran_ac<-mean(null_acs-obs_ac)

	null_mms<-rep(NA,nreps)
	for(i in 1:nreps){
	        sran[t_mm]<-sample(sran[t_mm],length(t_mm),replace=FALSE)
	        r_mm_surv<-which(trt$Treatment=="MM" & sran==1)
	        rcors_mms<-abs(cor(pcs[r_mm_surv,]))
	        rcorsbar_mms<-mean(rcors_mms[upper.tri(rcors_mms)][ttm[upper.tri(ttm)]==0])
	        null_mms[i]<-rcorsbar_mms
	}
	pv_ran_mm<-(1+sum(null_mms >= obs_mms))/(nreps+1)
	ts_mm<-obs_mms-obs_mm
	null_ran_mm<-mean(null_mms-obs_mm)

	## relative to linear selection
	## AC
	null_acs_lrg<-rep(NA,nreps)
	null_acs_lgb<-rep(NA,nreps)

	o_ac_rg<-glm(trt$Survival[t_ac] ~ trt$RG[t_ac],family="binomial")
	o_ac_gb<-glm(trt$Survival[t_ac] ~ trt$GB[t_ac],family="binomial")

	p_ac_rg<-trt$Survival[t_ac]
	p_ac_rg[is.na(p_ac_rg)==FALSE]<-o_ac_rg$fitted.values
	p_ac_rg[is.na(p_ac_rg)==TRUE]<-mean(o_ac_rg$fitted.values,na.rm=TRUE)

	p_ac_gb<-trt$Survival[t_ac]
	p_ac_gb[is.na(p_ac_gb)==FALSE]<-o_ac_gb$fitted.values
	p_ac_gb[is.na(p_ac_gb)==TRUE]<-mean(o_ac_gb$fitted.values,na.rm=TRUE)

	for(i in 1:nreps){
	        sri<-sample(t_ac,size=sum(trt$Survival[t_ac]==1,na.rm=TRUE),replace=FALSE,prob=p_ac_rg)
	        sran<-rep(0,length(t_ac)+length(t_mm))
	        sran[sri]<-1
	        r_ac_surv<-which(sran==1)
	        rcors_acs<-abs(cor(pcs[r_ac_surv,]))
	        rcorsbar_acs<-mean(rcors_acs[upper.tri(rcors_acs)][ttm[upper.tri(ttm)]==0])
	        null_acs_lrg[i]<-rcorsbar_acs

	        sri<-sample(t_ac,size=sum(trt$Survival[t_ac]==1,na.rm=TRUE),replace=FALSE,prob=p_ac_gb)
	        sran<-rep(0,length(t_ac)+length(t_mm))
	        sran[sri]<-1
	        r_ac_surv<-which(sran==1)
	        rcors_acs<-abs(cor(pcs[r_ac_surv,]))
	        rcorsbar_acs<-mean(rcors_acs[upper.tri(rcors_acs)][ttm[upper.tri(ttm)]==0])
	        null_acs_lgb[i]<-rcorsbar_acs
	}

	## RG
	pv_ls_rg_ac<-(1+sum(null_acs_lrg >= obs_acs))/(nreps+1)
	null_ls_rg_ac<-mean(null_acs_lrg-obs_ac)
	## GB
	pv_ls_gb_ac<-(1+sum(null_acs_lgb >= obs_acs))/(nreps+1)
	null_ls_gb_ac<-mean(null_acs_lgb-obs_ac)

	### MM
	null_mms_lrg<-rep(NA,nreps)
	null_mms_lgb<-rep(NA,nreps)

	o_mm_rg<-glm(trt$Survival[t_mm] ~ trt$RG[t_mm],family="binomial")
	o_mm_gb<-glm(trt$Survival[t_mm] ~ trt$GB[t_mm],family="binomial")

	p_mm_rg<-trt$Survival[t_mm]
	p_mm_rg[is.na(p_mm_rg)==FALSE]<-o_mm_rg$fitted.values
	p_mm_rg[is.na(p_mm_rg)==TRUE]<-mean(o_mm_rg$fitted.values,na.rm=TRUE)

	p_mm_gb<-trt$Survival[t_mm]
	p_mm_gb[is.na(p_mm_gb)==FALSE]<-o_mm_gb$fitted.values
	p_mm_gb[is.na(p_mm_gb)==TRUE]<-mean(o_mm_gb$fitted.values,na.rm=TRUE)


	for(i in 1:nreps){
	        sri<-sample(t_mm,size=sum(trt$Survival[t_mm]==1,na.rm=TRUE),replace=FALSE,prob=p_mm_rg)
	        sran<-rep(0,length(t_ac)+length(t_mm))
	        sran[sri]<-1
	        r_mm_surv<-which(sran==1)
	        rcors_mms<-abs(cor(pcs[r_mm_surv,]))
	        rcorsbar_mms<-mean(rcors_mms[upper.tri(rcors_mms)][ttm[upper.tri(ttm)]==0])
	        null_mms_lrg[i]<-rcorsbar_mms

	        sri<-sample(t_mm,size=sum(trt$Survival[t_mm]==1,na.rm=TRUE),replace=FALSE,prob=p_mm_gb)
	        sran<-rep(0,length(t_ac)+length(t_mm))
	        sran[sri]<-1
	        r_mm_surv<-which(sran==1)
	        rcors_mms<-abs(cor(pcs[r_mm_surv,]))
	        rcorsbar_mms<-mean(rcors_mms[upper.tri(rcors_mms)][ttm[upper.tri(ttm)]==0])
	        null_mms_lgb[i]<-rcorsbar_mms
	}
	## RG
	pv_ls_rg_mm<-(1+sum(null_mms_lrg >= obs_mms))/(nreps+1)
	null_ls_rg_mm<-mean(null_mms_lrg-obs_mm)
	## GB
	pv_ls_gb_mm<-(1+sum(null_mms_lgb >= obs_mms))/(nreps+1)
	null_ls_gb_mm<-mean(null_mms_lgb-obs_mm)

	## collect outpus
	out<-c(ts_ac,null_ran_ac,pv_ran_ac,
	       null_ls_rg_ac,pv_ls_rg_ac,null_ls_gb_ac,pv_ls_gb_ac,
		ts_mm,null_ran_mm,pv_ran_mm,
	       null_ls_rg_mm,pv_ls_rg_mm,null_ls_gb_mm,pv_ls_gb_mm)
	return(out)
}
