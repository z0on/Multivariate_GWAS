if (length(commandArgs(trailingOnly=TRUE))<1) {
	stop("

Compiles RDA_GWAS.R results, for all chromosomes, for a series of hold-out replicates.

Arguments: 

in=[filename]       List of files each listing per-chromosome RDA_GWAS.R outputs (.RData file names) 
                     for the same hold-out replicate. 
                     
                     For example, if we have 100 hold-out replicates, and per-chromosome output files 
                     from RDA_GWAS.R have rep1_,rep2_,...rep10_ in their filenames, then
                     
                        for r in `seq 1 100`;do ls *rep${r}_*.RData >rep${r}_files; done
                        ls rep*_files >allreps
                     
                      and use in=allreps in the call to this script.

traits=[filename]    Tab-delimited table of trait(s), same as used for RDA_GWAS.R Should contain
                     all samples, including the hold-out ones. Predictions will be plotted against 
                     the first column in this table.
                     
minz=2               minimal SNP z-score cutoff 

runGLMnet=TRUE       Whether to re-run elastic net on compiled data for each replicate.

forcePred=FALSE      forces back-prediction of the same dataset (for sanity check)

forceAlpha=-1        if non-negative, sets alpha parameter for glmnet regression (0 - ridge, 1 - lasso)
                     otherwise alpha will be chosen automatically based on back-prediction accuracy.

outfile=[filename]   output file name
                   

Mikhail Matz, matz@utexas.edu, July 2020

")

infile=grep("in=",commandArgs())
if (length(infile)==0) { stop ("specify list of replicates (in=filename)\nRun script without arguments to see all options\n") }
infile=sub("in=","", commandArgs()[gtfile])

traits =grep("traits=",commandArgs())
if (length(traits)==0) { stop ("specify traits file (traits=filename)\nRun script without arguments to see all options\n") }
traits =sub("traits=","", commandArgs()[traits])

minz =grep("minz=",commandArgs())
if(length(minz)>0) { nsites=as.numeric(sub("minz=","", commandArgs()[minz])) } else { minz=2 }

outfile =grep("outfile=",commandArgs())
if (length(outfile)==0) { stop ("specify output file name (outfile=filename)\nRun script without arguments to see all options\n") }
outfile =sub("outfile=","", commandArgs()[outfile])

runGLMnet=TRUE
if(length(grep("runGLMnet=F",commandArgs()))>0) { runGLMnet=FALSE } 

forcePred=FALSE
if(length(grep("forcePred=T",commandArgs()))>0) { forcePred=TRUE } 

forceAlpha=-1
fa=grep("forceAlpha=",commandArgs())
if(length(fa)>0) { forceAlpha=as.numeric(sub("forceAlpha=","", commandArgs()[fa])) } 

require(vegan)
require(ggplot2)
require(glmnet)
require(scales)
require(RColorBrewer)

# function to test the accuracy of prediction for a specific alpha
gnets=function(dat,trait,alpha=0.5) {
	index=sample(1:nrow(dat),round(0.75*nrow(dat)))
	train=dat[index,]
	tr.tr=trait[index]
	test=dat[-index,]
	tr.tst=trait[-index]
	trains.CV = cv.glmnet(train, tr.tr, nfolds=10,alpha=alpha,family="gaussian")
	# The definition of the lambda parameter:
	lambda.trains = trains.CV$lambda.min
	# Fit the elastic net predictor to the training data
	trains = glmnet(train, tr.tr, family="gaussian", alpha=alpha, lambda=lambda.trains)
	# Arrive at an estimate of of DNAmAge
	preds=predict(trains,test,type="response",s=lambda.trains)
#	plot(preds~tr.tst,xlab="true",ylab="predicted")
#	abline(0,1,col="red")
	r2=summary(lm(preds~tr.tst))$adj.r.squared
#	mtext(paste("a=",alpha," n=",ncol(dat),"  R2=",round(r2,2)))
	return(list(r2,trains,preds))	
}

#------------- reading and aligning data

# for 10 replicates:
# mv ../chr*bl_rep*_*RData .
# for r in `seq 1 10`;do ls *bl_rep${r}_10.RData >rep${r}_bl; done
# ls rep*_bl >reps_bl

# setwd("~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS/PDruns")
     # infile = "reps_bl"
     # traits = "../bleach.traits"
     # outfile="bleach_rerunGlm.RData"
	 # minz=2
	 # forcePred=FALSE
	 # forceAlpha=-1
  	# runGLMnet=TRUE

reps=scan(infile,what="character")		
traits=read.table(traits,header=T,stringsAsFactors=F)
row.names(traits)=traits$sample
traits$sample=NULL
tnames=names(traits)

rr=reps[1]
preds=list()
for (r in 1:length(reps)){
	rr=reps[r]

# ----- reading all chromosome data
	infiles=scan(rr,what="character")
	gts.test=c();gts.use=c();outs=c()
	for (f in 1:length(infiles)) {
		ll=load(infiles[f])
		if(f==1) {
			goods.test=0
			goods.use=names(gt.s)
			if(length(gt.test)>1) { 
				goods.test=names(gt.test) 
				traits.test=data.frame(traits[goods.test,])
				traits.use=data.frame(traits[goods.use,])
				
			} else { 
				traits.test=traits.use=traits
			}
			row.names(traits.test)=goods.test
			colnames(traits.test)=tnames
		}
		if(goods.test[1]==0) { gt.test=gt.s	}	
		gts.test[[f]]=gt.test
		gts.use[[f]]=gt.s
		outs[[f]]=out
	}
	out=do.call(rbind,outs)
	gt.test=do.call(rbind,gts.test)
	gt.s=do.call(rbind,gts.use)

#------- running elastic net
		
	if (runGLMnet) {	
		# screening through alpha parameters
		if(forceAlpha>=0) { alpha=forceAlpha } else {
			message("rep ",r, ": glmnet: screening alphas...")
			as=c()
			for (a in seq(0,1,0.2)) {
				message("     ",a)
				r2s=c()
				for (i in 1:5) {
					p=gnets(t(gt.s),sample.scores,alpha=a)
					r2s=append(r2s,p[[1]])
				}
				as=data.frame(cbind(as,r2s))
			}
			meanr2=apply(as,2,mean)
			alpha=seq(0,1,0.2)[which(meanr2==max(meanr2))]
			message("rep ",r, ": alpha=",alpha)
		}

		net.CV = cv.glmnet(t(gt.s), sample.scores, nfolds=10,alpha=alpha,family="gaussian")
		lambda = net.CV$lambda.min
		model = glmnet(t(gt.s), sample.scores, family="gaussian", alpha=alpha, nlambda=100)
		betas=coef(model,s=lambda)[,1][-1]
		intercept=coef(model,s=lambda)[,1][1]
		length(betas)
		
		out$beta.rr=betas
		out$intercept.rr=intercept
	}
	save(out,gt.s,gt.test,sample.scores,file=paste(rr,outfile,sep="_"))

#------------- predicting 

	if(forcePred) { 
		gt.test=gt.s 
		traits.test=traits[colnames(gt.s),]
		}
	tt=traits.test[,1]
	tt[which(is.na(tt))]=mean(tt,na.rm=T)
	
	pree=c()
	for (i in 1:ncol(gt.test)) {
		pree[i]=sum(out$beta.rr[abs(out$zscore)>=minz]*gt.test[abs(out$zscore)>=minz,i])+out$intercept[1]
	}
		
	pree3=c()
	for (i in 1:ncol(gt.test)) {
		pree3[i]=sum(out$beta[abs(out$zscore)>=minz]*gt.test[abs(out$zscore)>=minz,i])
	}

	# par(mfrow=c(1,2))
	# plot(pree~traits.test[,1],main=paste(rr,"reg.reg. betas"))
	# mtext(round(cor(pree,tt),2))
	# plot(pree3~traits.test[,1],main=paste(rr,"simple betas"))
	# mtext(round(cor(pree3,tt),2))

	preds[[r]]=data.frame(cbind(rr=pree,lm=pree3,true=traits.test[,1]))
	row.names(preds[[r]])=row.names(traits.test)
	preds[[r]]$rep=rr

}

N=sum(abs(out$zscore)>=minz)
allpreds=do.call(rbind,preds)
jitter=0.01*max(allpreds$true)
allpreds$re.lm=rescale(allpreds$lm,range(allpreds$true))+rnorm(nrow(allpreds),0,jitter)
allpreds$re.rr=rescale(allpreds$rr,range(allpreds$true))+rnorm(nrow(allpreds),0,jitter)
allpreds$re.true=allpreds$true+rnorm(nrow(allpreds),0,jitter)

#cols=WGCNA::labels2colors(allpreds$rep)
getPal = colorRampPalette(brewer.pal(9, "Set1"))
repColors= getPal(length(reps))

pdf(paste(outfile,".pdf",sep=""),width=8,height=3.2)
par(mfrow=c(1,3))
plot(re.lm~re.true,allpreds,main=paste("z:",minz," Nsnps:",N," simple betas"),col= repColors,pch=16)
mtext(round(cor(allpreds$lm,allpreds$true),2))
plot(re.rr~re.true,allpreds,main=paste("z:",minz," Nsnps:",N," rr betas"),col= repColors,pch=16)
mtext(round(cor(allpreds$rr,allpreds$true),2))
plot(lm~rr,allpreds,col= repColors,main="prediction comparison",pch=16)
dev.off()
