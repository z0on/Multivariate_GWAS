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

#--------- manhattan for arbitrary measure

manhat.y=function(mh,y) {
	mh=mh[grep("chr",mh$chrom),]
	mh$chrN=gsub("\\D","",mh$chrom)
	mh$chrom=factor(mh$chrom,levels=unique(mh$chrom)[order(as.numeric(unique(mh$chrN)))])
	mh$zscore=mh[,y]
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
#		ylim(min(abs(mh$zscore)),max(c(5,max(abs(mh$zscore)))))+
#		ggtitle(paste(gtfile,"z > 3"))+
		xlab("position,Mb")+facet_grid(~chrom,scale="free_x",space="free_x")
	return(pp)
}

#--------- manhattan for weighted pvals

manhat.y.wt=function(mh,y) {
	mh=mh[grep("chr",mh$chrom),]
	mh$chrN=gsub("\\D","",mh$chrom)
	mh$chrom=factor(mh$chrom,levels=unique(mh$chrom)[order(as.numeric(unique(mh$chrN)))])
	mh$zscore=mh[,y]
	mh$logp=mh$logp.wt
	mh$logp.adj=mh$logp.adj.wt
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
#		ylim(min(abs(mh$zscore)),max(c(5,max(abs(mh$zscore)))))+
#		ggtitle(paste(gtfile,"z > 3"))+
		xlab("position,Mb")+facet_grid(~chrom,scale="free_x",space="free_x")
	return(pp)
}

#--------- selecting best alpha in glmnet

# cros-validation by same-data resampling
gnets0=function(dat,trait,alpha=0.5) {
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


