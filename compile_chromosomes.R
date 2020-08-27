if (length(commandArgs(trailingOnly=TRUE))<1) {
	options(warning.length=8000)
	stop("

Compiles RDA_GWAS.R results for all chromosomes, for a single hold-out replicate (or full sample set).
Use compile_replicates.R after this script to summarize prediction accuracy in hold-out replicates.

Arguments: 

in=[filename]        List of per-chromosome RDA_GWAS.R '_gwas.RData' outputs
                     for the same hold-out replicate. 

traits=[filename]    RData bundle 'traits_etc_*.RData' containing stuff for the hold-out replicate 
                     This is the second file saved by RDA_GWAS.R (same for all chromosomes)
                     
gt=[filename]        list of genotype filenames, corresponding to chromosomes (same as used by RDA_GWAS.R).
					 If left unspecified, elastic net regression (GLMnet) will NOT be rerun.                     

gt.samples=[filename]   single-column list of sample names corresponding to the genotype files.
                           (same as for RDA_GWAS run)


runGLMnet=TRUE       Whether to re-run elastic net regression on compiled data. If not, betas from 
                     per-chromosome elastic net regressions will be used.

plots=TRUE           Whether to plot the whole-genome Manhattan plot (for all SNPs with z-score > 2) and predictions.

forceAlpha=-1        if non-negative, sets alpha parameter for elastic net regression (0 - ridge, 1 - lasso)
                     otherwise alpha will be chosen automatically based on prediction accuracy of the hold-out samples.

Output: RData bundle containing data frame \"gwas\" (zscores, pvalues, betas.lm, betas.rr), genotypes for \'use\' set (gt.use) 
genotypes for held-out \'test\' set (gt.test), and manhattan data for all sites (\'manh\')                 

Mikhail Matz, matz@utexas.edu, July 2020

") }

infile=grep("in=",commandArgs())
if (length(infile)==0) { stop ("specify list of replicates (in=filename)\nRun script without arguments to see all options\n") }
infile=sub("in=","", commandArgs()[infile])

bams =grep("gt.samples=",commandArgs())
if (length(bams)==0) { stop ("specify file listing sample names (gt.samples=filename)\nRun script without arguments to see all options\n") }
bams =sub("gt.samples=","", commandArgs()[bams])

runGLMnet=TRUE
traits =grep("traits=",commandArgs())
if (length(traits)==0) { 
	stop("specify traits etc. file\nRun script without arguments to see all options\n") 
	}
traits =sub("traits=","", commandArgs()[traits]) 

gts=grep("gts=",commandArgs())
if (length(gts)==0) { 
	print("no list of files containing per-chromosome posterior allele counts supplied\n elastic net will NOT be rerun \nRun script without arguments to see all options\n") 
	runGLMnet=FALSE
	} else { gts =sub("in=","", commandArgs()[gts]) }

if(length(grep("runGLMnet=F",commandArgs()))>0) { runGLMnet=FALSE } 

plotManhattan=TRUE
if(length(grep("plotManhattan=F",commandArgs()))>0) { plotManhattan=FALSE } 

forceAlpha=-1
fa=grep("forceAlpha=",commandArgs())
if(length(fa)>0) { forceAlpha=as.numeric(sub("forceAlpha=","", commandArgs()[fa])) } 

require(dplyr)
require(glmnet)
require(ggplot2)
require(data.table)
require(scales)
options(datatable.fread.datatable=FALSE)

source("RDA_GWAS_functions.R")

#----------- debug params

# setwd("~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS/")
  # infile = "gwass0"
  # gts="gts"
  # traits="traits_etc_0.RData"
  # forceAlpha=-1
  # runGLMnet=TRUE
  # plots=TRUE
  # bams = "bams.qc"

# ----- reading all chromosome data

infiles =scan(infile,what="character")
bams=scan(bams,what="character")
#removing path
bams=sub(".+/","",bams)
#removing extension
bams=sub("\\..+","",bams)

if(runGLMnet) {
 gtss=scan(gts,what="character")
 ll=load(traits)
}

gts.test=list()
gts.use=list()
gwas.c=list()
message("reading genotypes...")
for (f in 1:length(infiles)) {
	ll=load(infiles[f])
	gwas.c[[f]]=gwas
	if(runGLMnet) {
		gt=fread(gtss[f],nThread=4)
		# removing possibly duplicated sites
		gt=distinct(gt,paste(gt[,1],gt[,2],sep=":"),.keep_all=T)
		row.names(gt)=paste(gt[,1],gt[,2],sep=":")
		gt=gt[,-c(1,2,ncol(gt))]
		colnames(gt)=bams
		gts.test[[f]]=gt[,goods.test]
		gts.use[[f]]=gt[,goods.use]
	}
}
gt.test=do.call(rbind,gts.test)
gt=do.call(rbind,gts.use)
gwas=do.call(rbind,gwas.c)

# recalcuating logp.adjusted
gwas$logp.adj=-log(p.adjust(10^(-gwas$logp),method="BH"),10)

chosen=which(gwas$ldpruned==1)

#---- manhattan

if(plots) {
    pdf(paste(infile,"_plots.pdf",sep=""))
	plot(manhat(gwas[abs(gwas$zscore)>2,]))
}

#------- running elastic net on combined data
		
if(runGLMnet) {
	if(forceAlpha>=0) { alpha=forceAlpha } else {
		message("glmnet: screening alphas...",appendLF=FALSE)
		r2s=c()
		for (a in seq(0,1,0.1)) {
			message("     ",a)
			p=gnets(scale(t(gt[chosen,]),scale=F),scale(t(gt.test[chosen,]),scale=F),sample.scores,traits.test[,1],alpha=a)
			r2s=append(r2s,p[[1]])
		}
		alpha=seq(0,1,0.1)[which(r2s==max(r2s))]
	}
	message("alpha ",alpha)
	
	net.CV = cv.glmnet(t(gt[chosen,]), sample.scores, nfolds=10,alpha=alpha,family="gaussian")
	#plot(net.CV)
	lambda = net.CV$lambda.min
	model = glmnet(scale(t(gt[chosen,]),scale=F), sample.scores, family="gaussian", alpha= alpha, nlambda=100)
	# plot(model)
	# print(model)
	# str(coef(model))
	
	gwas$beta.rr=0
	preds.rr=predict(model,t(gt.test[chosen,]),type="response",s=lambda)
	betas=as.matrix(coef(model,s=lambda))[-1]
	gwas$beta.rr[chosen]=betas
}

#--------------- predictions of hold-out samples (or back-prediction of the same)

goodst=row.names(traits.test)[which(!is.na(traits.test[,1]))]
if(is.character(traits.test[goodst[1],1])) { 
	tt=as.numeric(as.factor(traits.test[goodst,1])) 
} else { tt=as.numeric(traits.test[goodst,1]) }
gt.test=gt.test[,goodst]

# scanning numbers of good SNPs, looking for best prediction
snps=row.names(gwas[chosen,])[order(abs(gwas[chosen,"zscore"]),decreasing=T)]
pred.rr=c();pred.lm=c()
# head(out)
# plot(beta.rr~beta,out)
message("scanning for best number of SNPs...")
ns=unique(round(10^(seq(0.1,log(length(snps),10),length.out=30))))
pb=txtProgressBar(0,length(ns))
zr=vector("list", length = length(ns));cp=0;maxcp=0
for (j in 1:length(ns)) {
	s=ns[j]
	chosen2=snps[1:s]
	betas=gwas[chosen2,"beta"];
	gtt=gt.test[chosen2,]
	for (i in 1:ncol(gtt)) {
		pred.lm[i]=sum(betas*gtt[,i])
	}	
	preds=data.frame(cbind(lm=pred.lm,true=tt))
	row.names(preds)=goodst
	zr[[s]]=preds
	setTxtProgressBar(pb,j)
}

zr2=c();Ns=c()
for (s in ns) {
	 Ns=append(Ns,s)
	 allpreds=zr[[s]]
	 zr2=append(zr2,cor(allpreds$lm,allpreds$true))
}
if(plots) {
 	plot(zr2~Ns,xlab="N(SNPs)",ylab="prediction R",log="x")
	lines(zr2~Ns)
}

best=which(zr2==max(zr2))[1]
allpreds=zr[[ns[best]]]
allpreds$rr=preds.rr[,1]
allpreds$rr=rescale(allpreds$rr,range(allpreds$true))
allpreds$lm=rescale(allpreds$lm,range(allpreds$true))

if(plots) {
	N=ns[best]
	plot(lm~true,allpreds,main=paste("Nsnps:",N," lm"),ylab="predicted",xlab="observed")
	mtext(paste("R2 =",round(cor(allpreds$lm,allpreds$true)^2,2)))
	N=sum(gwas[chosen,"beta.rr"]!=0)
	plot(rr~true,allpreds,main=paste("Nsnps:",N," glmnet"),ylab="predicted",xlab="observed")
	mtext(paste("R2 =",round(cor(allpreds$rr,allpreds$true)^2,2)))	
	dev.off()
}

save(gwas,file=paste(infile,"_gwas.RData",sep=""))
save(allpreds,file=paste(infile,"_predictions.RData",sep=""))
