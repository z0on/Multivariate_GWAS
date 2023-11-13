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
  require(ggplot2)
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
	
Calculates SNP scores, empirical p-values, betas, and r-squares for RDA-based GWAS. 
Any number of covariates can be removed wihtout loss of power.

Per-SNP values are computed against the first constrained ordination axis. 
Several correlated traits can be used to define that axis.
 
arguments (all tables are space-delimited, can be .gz): 

gt=[filename]        genotypes: table of minor allele counts 
                     (rows - loci, columns - samples) the first two columns must be chromosome, 
                     position. Header line must be present (chr, pos, sample names)
				
covars.e=[filename]      table of NON-GENETIC covariates (e.g. sampling site, age of individual)
                         (rows - samples, columns - covariates).
                         These will be regressed out of traits.
                         First column must be sample names. Header line must be present (sample, covariates). 
                         May not fully match the genotype table. Rows containing NA will be removed.

covars.g=[filename]     table of GENETIC covariates (e.g. sequencing batch, read depth, first couple of genetic PCs)
                        (rows - samples, columns - covariates).
                        These will be regressed out of genotypes.
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
							
badsites=[filename]  list of sites to exclude, one per line in format chr1:12345

hold.out=[filename]  list of sample names to hold out from the whole analysis for testing the predictors.

plots=TRUE           whether to plot diagnostic plots ([out]_plots.pdf)

nsites=5500000       number of sites to compute FDR (for Manhattan plot)

prune.dist=50000     pruning distance (selected SNPs must be at least that far apart). Alternatively, a RData bundle containing
                     object rdlm, output of LDq.R script - distance to R2 dropoff below 0.1.

Output:              An RData bundle containing results table for pruned SNPs (gwas) with zscores, pvalues, 
                     betas and R2, their genotypes (gt.s), genotypes of hold-out samples (if any) at the selected
                     SNPs (gt.test), manhattan plot data for all sites (manh), and sample scores 
                     in the ordination (sample.scores).
   
Mikhail Matz, matz@utexas.edu, July 2020

")
}


gtfile =grep("gt=",commandArgs())
if (length(gtfile)==0) { stop ("specify genotype file (gt=filename)\nRun script without arguments to see all options\n") }
gtfile=sub("gt=","", commandArgs()[gtfile])

covs.e =grep("covars.e=",commandArgs())
if (length(covs.e)==0) { message ("running without environmental covariates\n");covs.e=0 } else { covs.e =sub("covars.e=","", commandArgs()[covs.e]) }

covs.g =grep("covars.g=",commandArgs())
if (length(covs.g)==0) { message ("running without genetic covariates\n");covs.g=0 } else { covs.g =sub("covars.g=","", commandArgs()[covs.g]) }

traits =grep("traits=",commandArgs())
if (length(traits)==0) { stop ("specify traits file (traits=filename)\nRun script without arguments to see all options\n") }
traits =sub("traits=","", commandArgs()[traits])

ibs =grep("gdist=",commandArgs())
if (length(ibs)==0) { stop ("specify genetic distance file (gdist=filename)\nRun script without arguments to see all options\n") }
ibs =sub("gdist=","", commandArgs()[ibs])

gdist.samples =grep("gdist.samples=",commandArgs())
if (length(gdist.samples)==0) { stop ("specify file listing samples for genetic distances (gdist.samples=filename)\nRun script without arguments to see all options\n") }
gdist.samples =sub("gdist.samples=","", commandArgs()[gdist.samples])

badss =grep("badsites=",commandArgs())
if(length(badss)>0) { badsites=as.character(sub("badsites=","", commandArgs()[badss])) } else { badsites=0 }

if(length(grep("plots=F",commandArgs()))>0) { plots=FALSE } else { plots=TRUE }
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
# source("RDA_GWAS_functions.R")
options(datatable.fread.datatable=FALSE)

#---- reading and aligning data

      # setwd("~/Dropbox/quigley_rda_gwas_rerun_jan2023/")
      #   gtfile = "cleanGT_ten.tab"
      #   covs.g = "nas_pops.tab"
      #   covs.e = 0
      #   traits = "prepost_shuffled.tab"
      #   gdist.samples = "cors.bams"
      #   ibs="corMatrix.tab"
      #   plots=T
      #   nsites=5500000
      #   prune.dist=50
      #   hold.out="holdout12.shuffled"
      #   badsites=0

outfile=paste(sub("\\..+","",gtfile),sub("\\..+","",traits),sub("\\..+","",hold.out),sep=".")

bams=scan(gdist.samples,what="character")
#removing path
bams=sub(".+/","",bams)
#removing extension
bams=sub("\\..+","",bams)

ibs = as.matrix(read.table(ibs))
if(length(bams)!=nrow(ibs)) { stop("genetic distance matrix and gdist.sample file don't seem to match\n")}
dimnames(ibs)=list(bams,bams)
goods.ibs=row.names(na.omit(ibs))

traits=fread(traits)
row.names(traits)=traits$sample
traits$sample=NULL
tnames=names(traits)
goods.traits=row.names(na.omit(traits))
traits=data.frame(scale(traits,scale=FALSE))
trcols=colnames(traits)

if(covs.e!=0) {
	covars.e=fread(covs.e)
	row.names(covars.e)=covars.e$sample
	covars.e$sample=NULL
	goods.covars.e=row.names(na.omit(covars.e))
}

if(covs.g!=0) {
	covars.g=fread(covs.g)
	row.names(covars.g)=covars.g$sample
	covars.g$sample=NULL
	goods.covars.g=row.names(na.omit(covars.g))
}

goods.covars=goods.traits
if(covs.e!=0 & covs.g!=0 ) { 
  goods.covars=intersect(goods.covars.e,goods.covars.g) 
  } else { 
    if(covs.e!=0 ) { 
      goods.covars=goods.covars.e 
    } else { 
        if (covs.g!=0) {
          goods.covars=goods.covars.g 
          } 
    }
}

#---- loading genotypes

message("reading genotypes...",appendLF=FALSE)
gt=as.data.frame(fread(gtfile,nThread=4))

message("done")
# removing possibly duplicated sites
gt=distinct(gt,gt[,1:2],.keep_all=T)
row.names(gt)=paste(t(gt[,1]),t(gt[,2]),sep=":")
if(length(grep("paste",names(gt)[ncol(gt)]))==1) { gt=gt[,-ncol(gt)]}
gt=gt[,-c(1,2)]
colnames(gt)=bams

# removing sites with NAs and invariable sites
message("removing sites with NA entries and invariant sites...")
gt=na.omit(gt)
sds=apply(gt,1,sd)
gt=gt[sds>0,]

message("removing low-freq (maf<0.05) sites...")
af=apply(gt,1,sum)
af=af/(2*ncol(gt))
gt=gt[af>0.05,]

# ----- reading bad sites, removing them

if(is.character(badsites)) {
  message("removing bad sites...\n")
  badsi=scan(badsites,what="character")
  gt=gt[!(row.names(gt) %in% badsi),]
}

# ---- aligning all data

goods=intersect(intersect(goods.covars,goods.traits),colnames(gt))
#colnames(gt)[!(colnames(gt) %in% goods.covars)]
#goods.covars[!(goods.covars %in% colnames(gt))]

if(covs.g!=0) { covars.g=data.frame(covars.g[goods,]); row.names(covars.g)=goods }
if(covs.e!=0) { covars.e=data.frame(covars.e[goods,]); row.names(covars.e)=goods }

gt=gt[,goods]
traits=data.frame(traits[goods,])
row.names(traits)=goods
colnames(traits)=tnames

# ------ dummifying genetic covariates

if(covs.g!=0) { 
	 covs=c()
	 for (ci in 1:ncol(covars.g)) {
	   if(is.factor(covars.g[,ci]) | is.integer(covars.g[,ci]) | is.character(covars.g[,ci])) { 
	     co=as.factor(covars.g[,ci])
	     dum=data.frame(model.matrix(~0+co))
	     dum=dum[,-ncol(dum)]
	     colnames(dum)=paste(colnames(covars.g)[ci],colnames(dum),sep="_")
	     colnames(dum)=sub("_co","_",colnames(dum))
	     colnames(dum)=sub("^co","",colnames(dum))
	     colnames(dum)=sub("_model.+","",colnames(dum))
	     covs=cbind(covs,dum)
	     
	   }  else { 
	     covs=cbind(covs,covars.g[,ci])
	     colnames(covs)[ncol(covs)]=colnames(covars.g)[ci]
	   }
	 }	 
}
row.names(covs)=row.names(covars.g)

# ---- splitting into train and test sets

goods.test=0
if (hold.out!=0) {
	goods.test=scan(hold.out,what="character")
	#removing path
	goods.test=sub(".+/","",goods.test)
	#removing extension
	goods.test=sub("\\..+","",goods.test)
	goods.use=goods[!(goods %in% goods.test)]
	traits.test=data.frame(traits[goods.test,])
	row.names(traits.test)=goods.test
	colnames(traits.test)=tnames
} else { goods.use=goods}

ibs=ibs[goods.use,goods.use]
if(covs.g!=0) { 
  covs=data.frame(covs[goods.use,])
#  colnames(covs)=covcols 
}
traits=data.frame(traits[goods.use,])
rownames(traits)=goods.use
colnames(traits)=trcols

message(nrow(traits)," samples used, ",sum(goods.test!=0)," held out")

if (hold.out>0) { 
	gt.test=gt[,goods.test]
} else { 
	gt.test=0 
	goods.test=goods.use
}
gt=gt[,goods.use]

# --- Regressing environmental covariates out of traits

trcols=colnames(traits)
trrows=row.names(traits)
if(covs.e!=0) { 
   	cov.formula="~";pl=""
	for (ci in 1:ncol(covars.e)) {
	   if(is.integer(covars.e[,ci]) | is.character(covars.e[,ci])) { 
			covars.e[,ci]=as.factor(covars.e[,ci])
		}
		cov.formula=paste(cov.formula,pl,"covars.e[goods.use,",ci,"]",sep="")
		pl="+"
	}
	for (t in ncol(traits)) {
		cov.formt=paste("traits[,",t,"]",cov.formula,sep="")
		traits[,t]=residuals(lm(formula(cov.formt)))
	}
}

# hist(af,breaks=40)

#------- computing RDA and SNP scores against CAP1

message("computing ordination and SNP scores...",appendLF=FALSE)

# using scaled genotypes
if(covs.g!=0) { 
	cap=capscale(ibs~.+Condition(as.matrix(covs)),data=traits,comm=scale(t(gt)))
} else {
	cap=capscale(ibs~.,data=traits,comm=scale(t(gt)))
} 

nullpc=round(ncol(gt)*0.75)

if(plots) { 
  png(paste(outfile,"_ordination.png",sep=""),height=3000,width=2800,res=400)
  # crosses - unadjusted (measured) per-sample trait values, points - adjusted based on genotypes (used for GWAS)
  plot(cap,scaling=1, choices=c(1,nullpc),display=c("wa","cn"),mgp=c(2.3,1,0), main="sample ordination",type="n") 
  points(cap,scaling=1, choices=c(1,nullpc),display="wa",pch=21,col="grey10",bg="skyblue")
  points(cap,scaling=1, choices=c(1,nullpc),display="lc",pch=4,col="coral")
  ordispider(cap,scaling=1, choices=c(1,nullpc),col="coral")
  text(cap, scaling=1, display="bp", col="firebrick", cex=1, lwd=3,choices=c(1,2))
  invisible(dev.off())
}

# correcting sign to have positive scores for positive traits[,1] values
flip=1
if(coef(lm(cap$CCA$u[,1]~traits[,1]))[2]<0) { flip=(-1)}

snp.scores=scale(scores(cap,display="species",choices=1))*flip
sample.scores=scores(cap,display="sites",choices=1)*flip
names(snp.scores)=row.names(cap$CCA$v)
sample.scores=rescale(sample.scores,range(traits[,1]))
cap$CCA$biplot[,1]=cap$CCA$biplot[,1]*flip

#------- q-q plot to see if there is signal

if(plots) {	
  pdf(paste(outfile,"_qq.pdf",sep=""),height=4.7,width=4)
#  png(paste(outfile,"_qq.png",sep=""),height=1900,width=1600,res=400)
  nulls=data.frame(scale(scores(cap,display="species",choices=c(nullpc:(nullpc+3))))) 
	nulls=stack(nulls)$values
	qqplot(nulls,snp.scores,cex=0.7,main="q-q",mgp=c(2.3,1,0))
	abline(0,1,col="red")	
	invisible(dev.off())
}

#----- "bullseye plot", colored rings: >2,3,4,5 SD from 0

if(plots) {
  pdf(paste(outfile,"_snpScores.pdf",sep=""),height=4,width=4)
#  png(paste(outfile,"_snpScores.png",sep=""),height=1600,width=1600,res=400)
	pscores=data.frame(scale(flip*scores(cap,display="species",choices=c(1,nullpc))))
	names(pscores)=c("CAP1","MDS120")
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
	pp=ggplot()+
		geom_point(data=pscores,aes(CAP1,MDS120,fill=z,color=z),shape = 21, colour = "grey20")+
		geom_segment(data=biplot,aes(x=x1,y=y1,xend=CAP1,yend=CAP2,color="cyan3"),arrow = arrow(length = unit(0.3, "cm")))+
		theme_bw()+coord_equal()+theme(legend.position="n")
	plot(pp)
	invisible(dev.off())
}

#----------computing p-values

nulls=data.frame(scale(scores(cap,display="species",choices=c((nullpc-round(2*nullpc/3)):nullpc))))
nulls=stack(nulls)$values
message("calculating pvalues...")
pvals=getEmpP(snp.scores,nulls)
names(pvals)=row.names(gt)

#--------- assembling results table

logp=-log(pvals,10)
padj=p.adjust(pvals,method="BH")
logp.adj=-log(padj,10)
chrom=sub(":.+","",row.names(gt))
pos=sub(".+:","",row.names(gt))
gwas=data.frame(cbind(chrom,pos,snp.scores,logp,logp.adj),stringsAsFactors=F)
names(gwas)[3]="zscore"
gwas$zscore =as.numeric(gwas$zscore)
gwas$logp=as.numeric(gwas$logp)
gwas$logp.adj=as.numeric(gwas$logp.adj)
gwas$pos=as.numeric(gwas$pos)
gwas$chrom=as.factor(gwas$chrom)

#---- manhattan

if(plots) {
  pdf(paste(outfile,"_manhattan.pdf",sep=""),height=3.5,width=4)
#  png(paste(outfile,"_manhattan.png",sep=""),height=1600,width=4800,res=400)
	plot(manhat(gwas[abs(gwas$zscore)>1,]))
	invisible(dev.off())
}

#----------- pruning SNPs by z-scores and distance, computing simple betas and r2s (by chromosome)
gwas$ldpruned=0
lds=FALSE
if(length(grep("RData",prune.dist))>0) { 
	 lds=TRUE
	 load(prune.dist)
	 lops=c()
	 if(!("chrom" %in% names(rdlm))) { rdlm$chrom=gwas$chrom[1]}
	 for (chr in levels(rdlm$chrom)) {
 	 	rdlms=subset(rdlm,chrom==chr)
		 yy=rdlms$rolldrop
		 xx=rdlms$rollpos
		 lo=loess(yy~xx,span=0.005)
		 lop=predict(lo,newdata=data.frame(xx=rdlm$rollpos))
		 lop[lop<250]=250
		 lopd=data.frame(cbind(chrom=chr,pos=rdlms$rollpos,ldist=lop*1.5))
		 lopd$pos=as.numeric(as.character(lopd$pos))
		 lopd$ldist=as.numeric(as.character(lopd$ldist))
	 }
	 lops=data.frame(rbind(lopd))
}

blips=which(abs(snp.scores)>2)
message("\n",length(blips)," blips; pruning...")
tops=gt[blips,]
tops.scores=snp.scores[blips]
chs=row.names(tops)
chroms=as.factor(sub(":.+","",chs))
poss=as.numeric(sub(".+:","",chs))

for (chr in levels(chroms)) {
	ps=c()
	message("\n     chromosome ",chr)
	tsc=tops.scores[chroms==chr]
	tops.c=tops[chroms==chr,]
	sorted=order(abs(tsc),decreasing=T)
message("      ",length(sorted)," blips")
	n0=length(sorted)
	pb=txtProgressBar(0,length(sorted))
	pos=poss[chroms==chr]
	if(lds) { 
		lops.c=subset(lops, chrom==chr)
		prdist=approx(x=lops.c$pos,y=lops.c$ldist,xout=pos)$y
# filling in possible NAs at tails
		tail1=prdist[min(which(!is.na(prdist)))]
		tail2=prdist[max(which(!is.na(prdist)))]
		prdist[is.na(prdist)]=tail1
		prdist[(max(which(!is.na(prdist)))):length(prdist)]=tail1
	} else { 
		prdist=rep(prune.dist,length(pos)) 
	}
	r2=c();b=c();p=c()
	while(length(pos)>0) {
		i=sorted[1]
		p=c(p,pos[i])
		if(is.na(prdist[i])) { 
			message("WARN: unpredictable position in ",chr,": ",pos[i]," - skipping")
			toprune=pos[i]
			pos=pos[-toprune]
			tsc=tsc[-toprune]
			tops.c=tops.c[-toprune,]
			sorted=order(abs(tsc),decreasing=T)
	#		message(" ", length(sorted)," remaining")
			setTxtProgressBar(pb,n0-length(sorted))
			next
			}
		posdelta=abs(pos-pos[i])
		toprune=which(posdelta<=prdist[i] & tsc*tsc[i]>0)
#		message(pos[i]," z:",round(tsc[i],1)," prdist:",round(prdist[i],0)," prune:",min(pos[toprune]), "-",max(pos[toprune])," (",max(pos[toprune])-min(pos[toprune]),")"," zs:", paste(round(tsc[toprune],1),collapse=" "), appendLF=FALSE)
		pos=pos[-toprune]
		tsc=tsc[-toprune]
		tops.c=tops.c[-toprune,]
		sorted=order(abs(tsc),decreasing=T)
#		message(" ", length(sorted)," remaining")
		setTxtProgressBar(pb,n0-length(sorted))
	}
	ps=p[order(p)]
#	chrs=append(chrs,rep(chr,length(p)))
	gwas$ldpruned[gwas$chrom==chr & gwas$pos %in% p]=1
}

chosen=which(gwas$ldpruned==1)
message(length(chosen)," blips left after pruning")

# replotting pruned manhattan plot
if (plots) { 
  pdf(paste(outfile,"_manhattanPruned.pdf",sep=""),height=3.5,width=4)
#  png(paste(outfile,"_manhattanPruned.png",sep=""),height=1600,width=4800,res=400)
  plot(manhat(gwas[chosen,]))
	invisible(dev.off())
}

#----- computing betas and r2s

message("\ncalculating betas...")
beta=function(x,y) { return(coef(lm(y~x))[2]) }
b=apply(t(gt[chosen,]),2,beta,y=sample.scores)
message("\ncalculating R2s...")
r2=apply(t(gt[chosen,]),2,cor,y=sample.scores)
r2=r2^2
gwas$beta=0
gwas$r2=0
gwas$beta[chosen]=b
gwas$r2[chosen]=r2

#------- regularized regression

if (hold.out==0) { 
	gt.test=gt
	traits.test=traits
}

save(traits.test,traits,sample.scores,goods.use,goods.test,file=paste("traits_etc_",hold.out,".RData",sep=""))

gwas$beta.rr=0
message("glmnet: screening alphas...",appendLF=FALSE)
r2s=c()
for (a in seq(0,1,0.1)) {
	message("     ",a)
	p=gnets(scale(t(gt[chosen,]),scale=F),scale(t(gt.test[chosen,]),scale=F),sample.scores,traits.test[,1],alpha=a)
	r2s=append(r2s,p[[1]])
}
alpha=seq(0,1,0.1)[which(r2s==max(r2s))[1]]
message(" ",alpha)

net.CV = cv.glmnet(t(gt[chosen,]), sample.scores, nfolds=10,alpha=alpha,family="gaussian")
#plot(net.CV)
lambda = net.CV$lambda.min
model = glmnet(scale(t(gt[chosen,]),scale=F), sample.scores, family="gaussian", alpha= alpha, nlambda=100)
# plot(model)
# print(model)
# str(coef(model))

preds=predict(model,t(gt.test[chosen,]),type="response",s=lambda)
betas=as.matrix(coef(model,s=lambda))[-1]
# Nb=sum(betas!=0)
# plot(preds~traits.test[,1],xlab="true",ylab="predicted",main=paste("glmnet predict: ",Nb," SNPs"))
# mtext(round(cor(preds,traits.test[,1]),2),cex=0.8)

gwas[chosen,"beta.rr"]=betas

save(gwas,file=paste(outfile,"_gwas.RData",sep=""))

#---------------- predicting test set (if not specified, predict same set)

if(plots){
  pdf(paste(outfile,"_predictions.pdf",sep=""),height=4,width=8)
#  png(paste(outfile,"_predictions.png",sep=""),height=1200,width=3000,res=400)
  par(mfrow=c(1,3))
  goodst=row.names(traits.test)[which(!is.na(traits.test[,1]))]
  if(is.character(traits.test[goodst[1],1])) { 
    tt=as.numeric(as.factor(traits.test[goodst,1])) 
  } else { tt=as.numeric(traits.test[goodst,1]) }
  if(sd(tt)==0) { stop("zero variation in trait in holdout set")}
  gt.test=gt.test[,goodst]
  
  # scanning numbers of good SNPs, looking for best prediction
  snps=row.names(gwas[chosen,])[order(abs(gwas[chosen,"zscore"]),decreasing=T)]
  pred.rr=c();pred.lm=c()
  # head(out)
  # plot(beta.rr~beta,out)
  message("scanning for best number of SNPs...")
  ns=unique(round(10^(seq(0.1,log(length(snps),10),length.out=30))))
  pb=txtProgressBar(0,length(ns))
  zr=vector("list", length = length(ns));cp=0;maxcp=0;j=3
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
  zr[[4]]
  zr2=c();Ns=c();s=1
  for (s in ns) {
    Ns=append(Ns,s)
    allpreds=zr[[s]]
    if(sd(allpreds$lm)==0) { 
      zr2=append(zr2,0) 
      } else { 
        zr2=append(zr2,cor(allpreds$lm,allpreds$true))
      }
  }
  plot(zr2~Ns,xlab="N(SNPs)",ylab="prediction R",log="x",mgp=c(2.3,1,0))
  lines(zr2~Ns)
  
  best=which(zr2==max(zr2))[1]
  allpreds=zr[[ns[best]]]
  N=ns[best]
  
  
  message("\n------------\nSimple lm prediction:")
  message("    N SNPs: ",N)
  message("    R2: ",round(max(zr2),2))
  print(head(gwas[snps[1:N],c("zscore","beta","beta.rr","r2")]))
  
  
  jitter=0.1*sd(allpreds$true)
  allpreds$re.lm=rescale(allpreds$lm,range(allpreds$true))+rnorm(nrow(allpreds),0,jitter)
  allpreds$re.true=allpreds$true+rnorm(nrow(allpreds),0,jitter)
  
  plot(re.lm~re.true,allpreds,main=paste("Nsnps:",N," lm"),ylab="predicted",xlab="observed",mgp=c(2.3,1,0))
  mtext(paste("R2 =",round(cor(allpreds$lm,allpreds$true)^2,2)),cex=0.6)
  
  # regularized prediction
  rr=c()
  for (i in 1:ncol(gt.test)) {
    rr[i]=sum(gwas[chosen,"beta.rr"]*gt.test[chosen,i])
  }
  allpreds=data.frame(cbind(rr,true=tt))
  N=sum(gwas[chosen,"beta.rr"]!=0)
  #	head(allpreds)
  allpreds$re.rr=rescale(allpreds$rr,range(allpreds$true))+rnorm(nrow(allpreds),0,jitter)
  allpreds$re.true=allpreds$true+rnorm(nrow(allpreds),0,jitter)
  
  plot(re.rr~re.true,allpreds,main=paste("Nsnps:",N," glmnet"),ylab="predicted",xlab="observed",mgp=c(2.3,1,0))
  mtext(paste("R2 =",round(cor(allpreds$rr,allpreds$true)^2,2)),cex=0.6)
  invisible(dev.off())
  
  message("\n------------\nRegularized regression prediction:")
  message("    N snps: ",N)
  message("    R2: ",round(cor(allpreds$rr,allpreds$true)^2,2))
  message("    top ten SNPs:")
  print(head(gwas[chosen,c("zscore","beta","beta.rr","r2")][order(gwas[chosen,"beta.rr"],decreasing=T),]))
  
}

