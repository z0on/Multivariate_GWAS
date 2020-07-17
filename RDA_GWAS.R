if (length(commandArgs(trailingOnly=TRUE))<1) {
	options(warning.length=8000)
	stop("
	
Calculates SNP scores, empirical p-values, betas, and r-squares for RDA-based GWAS. 
Any number of covariates can be removed wihtout loss of power.

Per-SNP values are computed against the first constrained ordination axis. 
Several correlated traits can be used to define that axis.
 
arguments (all tables are space-delimited, can be .gz): 

gt=[filename]        genotypes: table of minor allele counts 
                     (rows - loci, columns - samples) the first two columns must be chromosome, 
                     position header line must be present (chr, pos, sample names)
				
covars=[filename]    table of covariates to use (rows - samples, columns - covariates).
                     First column must be sample names. Header line must be present (sample, covariates). 
                     May not fully match the genotype table. Rows containing NA will be removed.

traits=[filename]    table of trait(s). First column must be sample names. 
                     There must be at least 2 columns (samples, 1 trait). Header line 
                     must be present (sample, traits). May not fully match the genotype table. 
                     Rows containing NA will be removed.

gdist=[filename]     matrix of genetic distances between samples listed in the genotype 
                     file (e.g. IBS matrix from angsd). No header line or other non-numeric columns.
					
gdist.samples=[filename]   single-column list of sample names corresponding to the genotype 
                           AND genetic distance matrix. Could be filenames with leading path and 
                           trailing extension (these will be removed).
							
hold.out=[filename]  list of samples to hold out from the whole analysis for testing the predictors.

outfile=[filename]    output file name

plots=TRUE           whether to plot diagnostic plots ([out]_plots.pdf)

nsites=5500000       number of sites to compute FDR (for Manhattan plot)

prune.dist=50000     pruning distance (selected SNPs must be at least that far apart)

Output:              An RData bundle containing results table for pruned SNPs (out) with zscores, pvalues, 
                     betas and R2, their genotypes (gt.s), genotypes of hold-out samples (if any) at the selected
                     SNPs (gt.test), manhattan plot data for all sites (manh), and sample scores 
                     in the ordination (sample.scores).
   
Mikhail Matz, matz@utexas.edu, July 2020

")
}


gtfile =grep("gt=",commandArgs())
if (length(gtfile)==0) { stop ("specify genotype file (gt=filename)\nRun script without arguments to see all options\n") }
gtfile=sub("gt=","", commandArgs()[gtfile])

covars =grep("covars=",commandArgs())
if (length(covars)==0) { stop ("specify covariates file (covars=filename)\nRun script without arguments to see all options\n") }
covars =sub("covars=","", commandArgs()[covars])

traits =grep("traits=",commandArgs())
if (length(traits)==0) { stop ("specify traits file (traits=filename)\nRun script without arguments to see all options\n") }
traits =sub("traits=","", commandArgs()[traits])

ibs =grep("gdist=",commandArgs())
if (length(ibs)==0) { stop ("specify genetic distance file (gdist=filename)\nRun script without arguments to see all options\n") }
ibs =sub("gdist=","", commandArgs()[ibs])

bams =grep("gdist.samples=",commandArgs())
if (length(bams)==0) { stop ("specify file listing samples for genetic distances (gdist.samples=filename)\nRun script without arguments to see all options\n") }
bams =sub("gdist.samples=","", commandArgs()[bams])

outfile =grep("outfile=",commandArgs())
if (length(outfile)==0) { stop ("specify output file name (outfile=filename)\nRun script without arguments to see all options\n") }
outfile =sub("outfile=","", commandArgs()[outfile])

if(length(grep("plots=FALSE",commandArgs()))>0) { plots=FALSE } else { plots=TRUE }
nsites =grep("nsites=",commandArgs())
if(length(nsites)>0) { nsites=as.numeric(sub("nsites=","", commandArgs()[nsites])) } else { nsites=5500000 }
prune.dist =grep("prune.dist=",commandArgs())
if(length(prune.dist)>0) { prune.dist=as.numeric(sub("prune.dist=","", commandArgs()[prune.dist])) } else { prune.dist=50000 }

hold.out =grep("hold.out=",commandArgs())
if(length(hold.out)>0) { hold.out=sub("hold.out=","", commandArgs()[hold.out]) } else { hold.out=0 }

require(dplyr)
require(vegan)
require(ggplot2)
require(glmnet)
require(data.table)
require(R.utils)
require(scales)
options(datatable.fread.datatable=FALSE)


#---- reading and aligning data

 # setwd("~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS")
        # gtfile = "chr7.postAlleles.gz"
        # covars = "mds2_14_25"
        # traits = "pd.traits"
        # bams = "bams.qc"
        # ibs="zz8.ibsMat"
        # outfile="c7pd.RData"
        # plots=TRUE
        # nsites=5500000
        # prune.dist=50000
        # hold.out="rep14_25"
 #       hold.out=0

bams=scan(bams,what="character")
#removing path
bams=sub(".+/","",bams)
#removing extension
bams=sub("\\..+","",bams)

ibs = as.matrix(read.table(ibs))
if(length(bams)!=nrow(ibs)) { stop("genetic distance matrix and gdist.sample file don't seem to match\n")}
dimnames(ibs)=list(bams,bams)
goods.ibs=row.names(na.omit(ibs))

traits=read.table(traits,header=T,stringsAsFactors=F)
row.names(traits)=traits$sample
traits$sample=NULL
tnames=names(traits)
goods.traits=row.names(na.omit(traits))

covars=read.table(covars,header=T,stringsAsFactors=F)
row.names(covars)=covars$sample
covars$sample=NULL
goods.covars=row.names(na.omit(covars))

goods=intersect(intersect(goods.ibs,goods.traits),goods.covars)
goods.test=0
if (hold.out!=0) {
	goods.test=scan(hold.out,what="character")
	goods.use=goods[!(goods %in% goods.test)]
	traits.test=data.frame(traits[goods.test,])
	row.names(traits.test)=goods.test
	colnames(traits.test)=tnames
} else { goods.use=goods}

ibs=ibs[goods.use,goods.use]
covars=data.frame(covars[goods.use,])

# dummifying covariates
covs=c()
for (ci in 1:ncol(covars)) {
#	if(is.integer(covars[,ci]) | is.character(covars[,ci])) { covars[,ci]=as.factor(covars[,ci]) }
	co=covars[,ci]
	if(is.integer(co) | is.character(co)) { co=as.factor(co) }
	covs=cbind(covs,model.matrix(~0+co))
}

traits=data.frame(traits[goods.use,])
row.names(traits)=goods.use
colnames(traits)=tnames

message(nrow(traits)," samples used, ",sum(goods.test!=0)," held out")
message("reading genotypes...",appendLF=FALSE)
gt=fread(gtfile,nThread=4)
message("done")
# removing possibly duplicated sites
gt=distinct(gt,paste(gt[,1],gt[,2],sep=":"),.keep_all=T)
row.names(gt)=paste(gt[,1],gt[,2],sep=":")
gt=gt[,-c(1,2,ncol(gt))]
colnames(gt)=bams
if (hold.out>0) { gt.test=gt[,goods.test] } else { gt.test=0 }
gt=gt[,goods.use]

#------- computing RDA and SNP scores against CAP1

if (plots) { pdf(paste(outfile,"_plots.pdf",sep="")) }

message("computing ordination and SNP scores...",appendLF=FALSE)
#cap=capscale(ibs~.+Condition(covs),data=traits)
cap=capscale(ibs~.+Condition(covs),data=traits,comm=t(gt))
message("done")

if(plots) { plot(cap,scaling=3,display=c("wa","cn"),mgp=c(2.3,1,0), main="sample ordination") }

# correcting sign to have positive scores for positive traits[,1] values
flip=1
if(coef(lm(cap$CCA$u[,1]~traits[,1]))[2]<0) { flip=(-1)}

snp.scores=as.vector(scale(cap$CCA$v[,1]))*flip
names(snp.scores)=row.names(cap$CCA$v)
sample.scores=as.vector(cap$CCA$u[,1])*flip
cap$CCA$biplot[,1]=cap$CCA$biplot[,1]*flip

#------- q-q plot to see if there is signal

if(plots) {
	
	nulls=data.frame(cap$CA$v[,107:110]) 
	nulls=data.frame(apply(nulls,2,scale))
	nulls=stack(nulls)$values
	
	qqplot(nulls,scale(snp.scores),cex=0.7,main="q-q",mgp=c(2.3,1,0))
	abline(0,1,col="red")	
}

#----- "caviar plots", colored rin: >5 SD from 0

if(plots) {
	
	pscores=data.frame(cbind(snp.scores,scale(cap$CA$v[,100])))
	names(pscores)=c("CAP1","MDS100")
	
	biplot=data.frame(cap$CCA$biplot)
	
	if(ncol(traits)<2) { biplot$CAP2=0 }
	biplot$x1=0
	biplot$y1=0
	biplot=biplot/(max(abs(biplot)/max(abs(pscores[,1]))))
	distfrom0=apply(pscores[,1:2],1,function(x){sqrt(x[1]^2+x[2]^2)})
	pscores$z=0
	pscores$z[distfrom0>2]=2
	pscores$z[distfrom0>3]=3
	pscores$z[distfrom0>4]=4
	pscores$z[distfrom0>5]=5
	pscores$z=factor(pscores$z)
	ggplot()+
		geom_point(data=pscores,aes(CAP1,MDS100,fill=z,color=z),shape = 21, colour = "grey20")+
		geom_segment(data=biplot,aes(x=biplot$x1,y=biplot$y1,xend=biplot$CAP1,yend=biplot$CAP2,color="cyan3"),arrow = arrow(length = unit(0.3, "cm")))+
		theme_bw()+coord_equal()+theme(legend.position="n")
}

#------- empirical pvalue function

getEmpP <- function(Obs, Null,nq=10000){
#Null=nulls;Obs=snp.scores
	pb=txtProgressBar(0,length(Obs))
	Null=abs(Null)
	Obs=abs(Obs)
	LN=length(Null)
	logq=1+log(seq(1/nq,1,1/nq),nq+1)
	qn=quantile(Null,logq)
	length(qn)
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

#----------computing p-values

nulls=data.frame(cap$CA$v[,20:120])
nulls=data.frame(apply(nulls,2,scale))
nulls=stack(nulls)$values
message("calculating pvalues...")
pvals=getEmpP(snp.scores,nulls)

#--------- manhattan plot (by chromosome)

logp=-log(pvals,10)
padj=p.adjust(pvals,method="BH",n=nsites)
logp.adj=-log(padj,10)
chrom=sub(":.+","",row.names(gt))
pos=sub(".+:","",row.names(gt))
manh=data.frame(cbind(chrom,pos,"zscore"=snp.scores,logp,logp.adj),stringsAsFactors=F)
manh$zscore =as.numeric(manh$zscore)
manh$logp=as.numeric(manh$logp)
manh$logp.adj=as.numeric(manh$logp.adj)
manh$pos.kb=as.numeric(manh$pos)/1000
manh$pos.Mb=manh$pos.kb/1000
manh$chrom=as.factor(manh$chrom)

if(plots) {
	for (chr in levels(manh$chrom)) { 
		mc=subset(manh,chrom==chr)
		pp=ggplot(mc,aes(pos.Mb,logp))+
			geom_point(shape = 21, colour = "grey20", aes(size=logp.adj,fill=zscore))+
		    scale_size_continuous(limits=c(0,3),breaks=c(0.3,0.6,1,1.3,2),labels=c(0.5,0.25,0.1,0.05,1e-2))+
			scale_fill_gradient(low="cyan3",high="coral")+
			theme_bw() + labs(size = "p.adj")+
			ylim(0,max(c(7,max(manh$logp))))+ggtitle(paste(chr,"raw"))+
			xlab("position,Mb")
		plot(pp)
	}
}

#----------- pruning SNPs by z-scores and distance, computing simple betas and r2s (by chromosome)

blips=which(abs(snp.scores)>2)
message("\n",length(blips)," blips; pruning...")
tops=gt[blips,]
tops.scores=snp.scores[blips]
chs=row.names(tops)
chroms=as.factor(sub(":.+","",chs))
poss=as.numeric(sub(".+:","",chs))
r2s=c();bs=c()
for (chr in levels(chroms)) {
	message("\n     chromosome ",chr)
	tsc=tops.scores[chroms==chr]
	tops.c=tops[chroms==chr,]
	sorted=order(abs(tsc),decreasing=T)
message("      ",length(sorted)," blips")
	pb=txtProgressBar(0,length(sorted))
	pos=poss[chroms==chr][sorted]
	r2=vector(mode="numeric",length=length(sorted));b=vector(mode="numeric",length=length(sorted))
	plus=minus=c(-2*prune.dist);s=100
	for (s in 1:length(sorted)) {
		i=sorted[s]
		setTxtProgressBar(pb,s)
	# skipping sites if closer than prune.dist to any of the already profiled ones
		ss=tsc[i]
#		if(abs(ss)<2) { message("\n===> low score: ",ss," s:",s," i:",i)}
		if (ss>=0) {
			if(min(abs(plus-pos[s]))<prune.dist) {
				b[i]=0
				r2[i]=0
				next
			} else { plus=append(plus,pos[s]) }
		} else {
			if(min(abs(minus-pos[s]))<prune.dist) {
				b[i]=0
				r2[i]=0
				next
			} else { minus=append(minus,pos[s]) }
		}
		mm=lm(sample.scores~t(tops.c[i,]))
		r2[i]=round(summary(mm)$adj.r.squared,4)
		b[i]=coef(mm)[2]
	}
	r2s=append(r2s,r2)
	bs=append(bs,b)
}

out=data.frame(cbind(zscore=tops.scores,pvals=pvals[blips],padj=padj[blips],beta=bs,r2=r2s),stringsAsFactors=F)
row.names(out)=row.names(gt[blips,])
out$chr=as.factor(sub(":.+","",row.names(gt[blips,])))
out$pos=as.numeric(sub(".+:","",row.names(gt[blips,])))
out=out[!(bs==0),]
out$logp=-log(out$pvals,10)
out$logp.adj=-log(out$padj,10)
message("\n",nrow(out)," independent blips collected")
gt.s=gt[row.names(out),]
if(hold.out!=0) { gt.test=gt.test[row.names(out),] }

#------ replotting pruned manhattan plot

if (plots) { 
	out$pos.Mb=out$pos/1e+6
	for (ch in levels(out$chr)) { 
		mc=out[out$chr==ch,]
		pp=ggplot(mc,aes(pos.Mb,logp))+
			geom_point(shape = 21, colour = "grey20", aes(size=logp.adj,fill=zscore))+
		    scale_size_continuous(limits=c(0,3),breaks=c(0.3,0.6,1,1.3,2),labels=c(0.5,0.25,0.1,0.05,1e-2))+
			scale_fill_gradient(low="cyan3",high="coral")+
			theme_bw() + labs(size = "p.adj")+
			ylim(min(mc$logp),max(c(7,max(mc$logp))))+ggtitle(paste(ch,"pruned"))+
			xlab("position,Mb")
		plot(pp)
	}
}

#------- regularized regression

# determining best alpha 
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

# screening through alpha parameters
message("glmnet: screening alphas...",appendLF=FALSE)
as=c()
for (a in seq(0,1,0.2)) {
	message("     ",a)
	r2s=c()
	for (i in 1:10) {
		p=gnets(t(gt.s),sample.scores,alpha=a)
		r2s=append(r2s,p[[1]])
	}
	as=data.frame(cbind(as,r2s))
}
meanr2=apply(as,2,mean)
alpha=seq(0,1,0.2)[which(meanr2==max(meanr2))]
message(" ",alpha)

net.CV = cv.glmnet(t(gt.s), sample.scores, nfolds=10,alpha=0,family="gaussian")
lambda = net.CV$lambda.min
model = glmnet(t(gt.s), sample.scores, family="gaussian", alpha=alpha, nlambda=100)
betas=coef(model,s=lambda)@x[-1]
intercept=coef(model,s=lambda)@x[1]

out$beta.rr=betas
out$intercept.rr=intercept

save(out,gt.s,gt.test,sample.scores,manh,file=outfile)

#---------------- predictng test set (if not specified, predict same set)

if(plots){
	if (hold.out==0) { 
		gt.test=gt.s
		traits.test=traits
		}

	goodst=row.names(traits.test)[which(!is.na(traits.test[,1]))]
	if(is.character(traits.test[goodst[1],1])) { 
		tt=as.numeric(as.factor(traits.test[goodst,1])) 
	} else { tt=as.numeric(traits.test[goodst,1]) }
	gt.test=gt.test[,goodst]

# scanning through z-score cutoffs, looking for best prediction
	zscan=seq(2,5,0.1);zr=vector("list", length = length(zscan))
	pred.rr=c();pred.lm=c()
	for (z in 1:length(zscan)) {
		minz=zscan[z]
		for (i in 1:ncol(gt.test)) {
			pred.rr[i]=sum(out$beta.rr[abs(out$zscore)>= minz]*gt.test[abs(out$zscore)>= minz,i])+out$intercept[1]
			pred.lm[i]=sum(out$beta[abs(out$zscore)>= minz]*gt.test[abs(out$zscore)>= minz,i])
		}	
		preds=data.frame(cbind(rr=pred.rr,lm=pred.lm,true=tt))
		row.names(preds)=goodst
		zr[[z]]=preds
	}
	zr2=c();Ns=c()
	for (z in 1:length(zscan)) {
		Ns=append(Ns,sum(abs(out$zscore)>=zscan[z]))
		allpreds=zr[[z]]
		zr2=append(zr2,cor(allpreds$lm,allpreds$true)^2)
	}
	plot(zr2~Ns,xlab="N(SNPs)",ylab="prediction R2")
	lines(zr2~Ns)

	best=which(zr2==max(zr2))
	allpreds=zr[[best]]

	N=Ns[best]
	jitter=0.01*max(allpreds$true)
	allpreds$re.lm=rescale(allpreds$lm,range(allpreds$true))+rnorm(nrow(allpreds),0,jitter)
	allpreds$re.rr=rescale(allpreds$rr,range(allpreds$true))+rnorm(nrow(allpreds),0,jitter)
	allpreds$re.true=allpreds$true+rnorm(nrow(allpreds),0,jitter)

	plot(re.lm~re.true,allpreds,main=paste("z:",zscan[best]," Nsnps:",N," simple betas"),ylab="predicted",xlab="observed")
	mtext(round(cor(allpreds$lm,allpreds$true)^2,2))
	plot(re.rr~re.true,allpreds,main=paste("z:",zscan[best]," Nsnps:",N," rr betas"),ylab="predicted",xlab="observed")
	mtext(round(cor(allpreds$rr,allpreds$true)^2,2))

	dev.off() 
}

