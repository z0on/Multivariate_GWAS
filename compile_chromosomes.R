if (length(commandArgs(trailingOnly=TRUE))<1) {
	options(warning.length=8000)
	stop("

Compiles RDA_GWAS.R results for all chromosomes, for a single hold-out replicate (or full sample set).
Use compile_replicates.R after this script to summarize prediction accuracy in hold-out replicates.

Arguments: 

in=[filename]        List of per-chromosome RDA_GWAS.R outputs (.RData file names) 
                     for the same hold-out replicate. 
                     
                     For example, if we have 50 hold-out replicates, and per-chromosome output files 
                     from RDA_GWAS.R have rep1_,rep2_,...rep10_ in their filenames, then
                     
                       >comp
                       for r in `seq 1 50`; do
                       ls *rep${r}_*.RData >rep${r}_files; 
                       echo \"Rscript compile_chromosomes.R in=rep${r}_files\">>comp;
                       done
 
                       and then run everything in comp.

runGLMnet=TRUE       Whether to re-run elastic net on compiled data.

plotManhattan=TRUE   Whether to plot the whole-genome Manhattan plot (for all SNPs with z-score > 3).

forceAlpha=-1        if non-negative, sets alpha parameter for glmnet regression (0 - ridge, 1 - lasso)
                     otherwise alpha will be chosen automatically based on back-prediction accuracy.

Output: RData bundle containing data frame \"out\" (zscores, pvalues, betas.lm, betas.rr), genotypes for \'use\' set (gt.use) 
genotypes for held-out \'test\' set (gt.test), and manhattan data for all sites (\'manh\')                 

Mikhail Matz, matz@utexas.edu, July 2020

") }

infile=grep("in=",commandArgs())
if (length(infile)==0) { stop ("specify list of replicates (in=filename)\nRun script without arguments to see all options\n") }
infile=sub("in=","", commandArgs()[infile])

runGLMnet=TRUE
if(length(grep("runGLMnet=F",commandArgs()))>0) { runGLMnet=FALSE } 

plotManhattan=TRUE
if(length(grep("plotManhattan=F",commandArgs()))>0) { plotManhattan=FALSE } 

forceAlpha=-1
fa=grep("forceAlpha=",commandArgs())
if(length(fa)>0) { forceAlpha=as.numeric(sub("forceAlpha=","", commandArgs()[fa])) } 

require(glmnet)
require(ggplot2)

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


# setwd("~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS/")
  # infile = "rep2_rf"
 # # forceAlpha=-1
  # runGLMnet=FALSE
# plotManhattan=TRUE

# ----- reading all chromosome data

infiles =scan(infile,what="character")
gts.test=list()
gts.use=list()
outs=list()
manhs=list()

for (f in 1:length(infiles)) {
	ll=load(infiles[f])
	gts.test[[f]]=gt.test
	gts.use[[f]]=gt.s
	outs[[f]]=out
	manhs[[f]]=manh
}
out=do.call(rbind,outs)
gt.test=do.call(rbind,gts.test)
gt.s=do.call(rbind,gts.use)
manh=do.call(rbind,manhs)

#------- running elastic net on combined data
		
if (runGLMnet) {	
	# screening through alpha parameters
	if(forceAlpha>=0) { alpha=forceAlpha } else {
		message("glmnet: screening alphas...")
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
		message("chosen alpha=",alpha)
	}

	net.CV = cv.glmnet(t(gt.s), sample.scores, nfolds=10,alpha=alpha,family="gaussian")
	lambda = net.CV$lambda.min
	model = glmnet(t(gt.s), sample.scores, family="gaussian", alpha=alpha, nlambda=100)
	betas=coef(model,s=lambda)[,1][-1]
	intercept=coef(model,s=lambda)[,1][1]
	
	out$beta.rr=betas
	out$intercept.rr=intercept
}

save(out,gt.s,gt.test,manh,file=paste(infile,".RData",sep=""))

#--------------- Manhattan plot (only for "chr" contigs)

if(plotManhattan) {
	pdf(paste(infile,"_manhattan.pdf",sep=""),width=15,height=3.5)
	manh$pos.Mb=manh$pos/1e+3
	mh=manh[abs(manh$zscore)>3,]
	mh=mh[grep("chr",mh$chrom),]
	mh$chrN=gsub("\\D","",mh$chrom)
	mh$chrom=factor(mh$chrom,levels=unique(mh$chrom)[order(as.numeric(unique(mh$chrN)))])
	sign=as.numeric(mh$zscore>0)
	sign[sign==0]=-1
	mh$slogp=mh$logp*sign
	mh$logp.adj[mh$logp.adj>2]=2
	mh$lpa2=(mh$logp.adj)^3
	pp=ggplot(mh,aes(pos.Mb,abs(zscore)))+
		geom_point(shape = 21, colour = "grey20", aes(size=lpa2,fill=slogp))+
	    scale_size_continuous(limits=c(0,8),breaks=c(0.3,0.6,1,1.3,2)^3,labels=c(0.5,0.25,0.1,0.05,1e-2))+
	  #  scale_size_continuous(limits=c(0,3),breaks=c(1,1.3,2),labels=c(0.1,0.05,1e-2))+
		scale_fill_gradient(low="cyan3",high="coral")+
		theme_bw() + labs(size = "p.adj")+theme(axis.text.x=element_text(angle=45, hjust=1))+
		ylim(min(abs(mh$zscore)),max(c(10,max(abs(mh$zscore)))))+
		ggtitle(paste(infile,"z > 3"))+
		xlab("position,Mb")+facet_grid(~chrom,scale="free_x",space="free_x")
	plot(pp)
	dev.off()
}
