#------- empirical pvalue function

getEmpP <- function(Obs, Null,nq=10000){
#Null=nulls;Obs=snp.scores
	pb=txtProgressBar(0,length(Obs))
	Null=abs(Null)
	Obs=abs(Obs)
	LN=length(Null)
	logq=1+log(seq(1/nq,1,1/nq),nq+1)
	qn=quantile(Null,logq)
	perq=length(Null)/nq
	lastnull=Null[Null>qn[nq-1]]
 	pvals=c()
	for (i in 1:length(Obs)) {
		qo=sum(qn<Obs[i])
		if(qo==0) { 
			pval=1
		} else {
			if (qo==nq) {
				pval=0.5/LN
			} else { 
				if (qo<nq-1) {
					pval = 1-(logq[qo])
				} else {
					pval = sum(lastnull>=Obs[i])/LN
				}
			}
		}
	    pvals[i]=pval
		setTxtProgressBar(pb,i)
	}
    return(pvals)
}

#--------- manhattan plot function

manhat=function(mh) {
	mh=mh[grep("chr",mh$chrom),]
	mh$chrN=gsub("\\D","",mh$chrom)
	mh$chrom=factor(mh$chrom,levels=unique(mh$chrom)[order(as.numeric(unique(mh$chrN)))])
	sign=as.numeric(mh$zscore>0)
	sign[sign==0]=-1
	mh$pos.Mb=mh$pos/1e+6
	mh$signed.logp=mh$logp*sign
	mh$logp.adj[mh$logp.adj>2]=2
	mh$lpa2=(mh$logp.adj)^3
	pp=ggplot(mh,aes(pos.Mb,abs(zscore)))+
		geom_point(shape = 21, colour = "grey20", aes(size=lpa2,fill=signed.logp))+
	    scale_size_continuous(limits=c(0,8),breaks=c(0.3,0.6,1,1.3,2)^3,labels=c(0.5,0.25,0.1,0.05,1e-2))+
	  #  scale_size_continuous(limits=c(0,3),breaks=c(1,1.3,2),labels=c(0.1,0.05,1e-2))+
		scale_fill_gradient(low="cyan3",high="coral")+
		theme_bw() + labs(size = "p.adj")+theme(axis.text.x=element_text(angle=45, hjust=1))+
		ylim(min(abs(mh$zscore)),max(c(5,max(abs(mh$zscore)))))+
#		ggtitle(paste(gtfile,"z > 3"))+
		xlab("position,Mb")+facet_grid(~chrom,scale="free_x",space="free_x")
	return(pp)
}

#--------- selecting best alpha in glmnet

# cros-validation using provided test set
gnets=function(train,test,tr.tr,tr.tst,alpha=0.5) {
	trains.CV = cv.glmnet(train, tr.tr, nfolds=10,alpha=alpha,family="gaussian")
	# The definition of the lambda parameter:
	lambda.trains = trains.CV$lambda.min
	# Fit the elastic net predictor to the training data
	trains = glmnet(train, tr.tr, family="gaussian", alpha=alpha, lambda=lambda.trains)
	# predict trait in test set
	preds=predict(trains,test,type="response",s=lambda.trains)
#	plot(preds~tr.tst,xlab="true",ylab="predicted",main=alpha)
#	abline(0,1,col="red")
	r2=summary(lm(preds~tr.tst))$adj.r.squared
#	mtext(paste("a=",alpha,"  R2=",round(r2,2)))
	return(list(r2,trains,preds))	
}

#========================


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
                     
gts=[filename]        list of genotype filenames, corresponding to chromosomes (same as used by RDA_GWAS.R)                  

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

if(length(grep("plots=F",commandArgs()))>0) { plots=FALSE } else { plots=TRUE }

traitfile=grep("traits=",commandArgs())
if (length(traitfile)==0) { 
	stop("specify traits etc. file\nRun script without arguments to see all options\n") 
	}
traitfile=sub("traits=","", commandArgs()[traitfile]) 

runGLMnet=TRUE
if(length(grep("runGLMnet=F",commandArgs()))>0) { runGLMnet=FALSE } 

gts=grep("gts=",commandArgs())
if (length(gts)==0) { 
	stop("no list of files containing per-chromosome posterior allele counts supplied\nRun script without arguments to see all options\n") 
	runGLMnet=FALSE
	} else { gts =sub("gts=","", commandArgs()[gts]) }

bams =grep("gt.samples=",commandArgs())
if (length(bams)==0) { stop ("specify file listing sample names (gt.samples=filename)\nRun script without arguments to see all options\n") }
bams =sub("gt.samples=","", commandArgs()[bams])

if(runGLMnet){
  forceAlpha=-1
  fa=grep("forceAlpha=",commandArgs())
  if(length(fa)>0) { forceAlpha=as.numeric(sub("forceAlpha=","", commandArgs()[fa])) } 
}

require(dplyr)
require(glmnet)
require(ggplot2)
require(data.table)
require(scales)
options(datatable.fread.datatable=FALSE)

# source("RDA_GWAS_functions.R")

#----------- debug params

# setwd("/work/01211/cmonstr/impute_gwas")
# infile = "gws_rep100_20"
# gts="gts"
# traitfile="traits_etc_rep100_20.RData"
# forceAlpha=-1
# runGLMnet=TRUE
# plots=TRUE
# bams = "samples"

# ----- reading all chromosome data

infiles =scan(infile,what="character")

if(runGLMnet) {
  bams=scan(bams,what="character")
  #removing path 
  bams=sub(".+/","",bams)
  #removing extension
  bams=sub("\\..+","",bams)

  gtss=scan(gts,what="character")
  gts.test=list()
  gts.use=list()
}

ll=load(traitfile)
gwas.c=list()
message("reading data...")
for (f in 1:length(infiles)) {
	message(infiles[f])
	ll=load(infiles[f])
	gwas.c[[f]]=gwas
	if(runGLMnet) {

	  #---- loading genotypes
	  gt=fread(gtss[f],nThread=4)
	  row.names(gt)=paste(gt[,1],gt[,2],sep=":")
	  gt=gt[,-c(1,2)]
	  colnames(gt)=bams
	  goods=intersect(row.names(gt), row.names(gwas))
	  gt=gt[goods,]
	  gwas=gwas[goods,]
	  
	  # message("done")
	  # # removing possibly duplicated sites
	  # gt=distinct(gt,paste(gt[,1],gt[,2],sep=":"),.keep_all=T)
	  # row.names(gt)=paste(gt[,1],gt[,2],sep=":")
	  # if(grep("paste",names(gt)[ncol(gt)])==1) { gt=gt[,-ncol(gt)]}
	  # gt=gt[,-c(1,2)]
	  # colnames(gt)=bams
	  # 
	  # # removing sites with NAs and invariable sites
	  # message("removing sites with NA entries and invariant sites...")
	  # gt=na.omit(gt)
	  # sds=apply(gt,1,sd)
	  # gt=gt[sds>0,]
	  # 
	  # message("removing low-freq (maf<0.05) sites...")
	  # af=apply(gt,1,sum)
	  # af=af/(2*ncol(gt))
	  # gt=gt[af>0.05,]
	  # 
	  # # ----- reading bad sites, removing them
	  # 
	  # if(is.character(badsites)) {
	  #   badsi=scan(badsites,what="character")
	  #   gt=gt[!(row.names(gt) %in% badsi),]
	  # }
	  # 
		gts.test[[f]]=gt[,goods.test]
		gts.use[[f]]=gt[,goods.use]
	}
}
if(runGLMnet) {
 gt.test=do.call(rbind,gts.test)
 gt=do.call(rbind,gts.use)
}

gwas=do.call(rbind,gwas.c)

# recalcuating logp.adjusted
gwas$logp.adj=-log(p.adjust(10^(-gwas$logp),method="BH"),10)

chosen=row.names(gwas)[which(gwas$ldpruned==1)]

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
	gwas[chosen,"beta.rr"]=betas
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
	 if(sd(allpreds$lm)==0) { 
	   zr2=append(zr2,0) 
	 } else { 
	   zr2=append(zr2,cor(allpreds$lm,allpreds$true))
	 }
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
