if (length(commandArgs(trailingOnly=TRUE))<1) {
	options(warning.length=8000)
	stop("

Compiles hold-out replicates using results of compile_chromosomes.R

Arguments: 

reps=[filename]       List of replicates (.RData files output by compile_chromosomes.R)                     

traits=[filename]    Tab-delimited table of trait(s), same as used for RDA_GWAS.R 
                                          
forcePred=FALSE      forces back-prediction of the same dataset (sanity check)

Output: pdf file [input filename]_predictions.pdf showing
        panel 1: scan through z-score cutoffs for best predictions (using lm betas)
        panel 2: predictions for hold-out samples based on lm betas
        panel 3: predictions based on regularized betas (glmnet)
        panel 4: comparison of simple and regularized predictions

Mikhail Matz, matz@utexas.edu, July 2020

") }


require(scales)
require(RColorBrewer)

infile=grep("reps=",commandArgs())
if (length(infile)==0) { stop ("specify list of replicates (.RData files output by compile_chromosomes.R)\nRun script without arguments to see all options\n") }
infile=sub("reps=","", commandArgs()[infile])

forcePred=FALSE
if(length(grep("forcePred=T",commandArgs()))>0) { forcePred=TRUE } 

traits =grep("traits=",commandArgs())
if (length(traits)==0) { stop ("specify traits file (traits=filename)\nRun script without arguments to see all options\n") }
traits =sub("traits=","", commandArgs()[traits])


 # infile="rf19"
 # traits="rf.traits"
 # forcePred=FALSE

traits=read.table(traits,header=T,stringsAsFactors=F)
row.names(traits)=traits$sample
tnames=names(traits)

outfile=paste(infile,"_predictions.pdf",sep="")
reps=scan(infile,what="character")		

zscan=seq(2,5,0.1);zr=vector("list", length = length(zscan))
for (r in 1:length(reps)){
	rr=reps[r]
	ll=load(rr)

	if(forcePred) { 
		gt.test=gt.s 
	}
	
	traits.test=traits[colnames(gt.test),]
	traits.test$sample=NULL

	goodst=row.names(traits.test)[which(!is.na(traits.test[,1]))]
	if(is.character(traits.test[goodst[1],1])) { 
		tt=as.numeric(as.factor(traits.test[goodst,1])) 
	} else { tt=as.numeric(traits.test[goodst,1]) }
	
	gt.test=gt.test[,goodst]
	
	pred.rr=c();pred.lm=c()
	for (z in 1:length(zscan)) {
		minz=zscan[z]
		for (i in 1:ncol(gt.test)) {
			pred.rr[i]=sum(out$beta.rr[abs(out$zscore)>= minz]*gt.test[abs(out$zscore)>= minz,i])+out$intercept[1]
			pred.lm[i]=sum(out$beta[abs(out$zscore)>= minz]*gt.test[abs(out$zscore)>= minz,i])
		}	
		preds=data.frame(cbind(rr=pred.rr,lm=pred.lm,true=tt))
		row.names(preds)=goodst
		preds$rep=rr
		zr[[z]][[r]]=preds
	}
}
pdf(outfile,width=11,height=3.2)
par(mfrow=c(1,4))
zr2=c();Ns=c()
for (z in 1:length(zscan)) {
	Ns=append(Ns,sum(abs(out$zscore)>=zscan[z]))
	allpreds=do.call(rbind,zr[[z]])
	zr2=append(zr2,cor(allpreds$lm,allpreds$true)^2)
}
#plot(zr2~zscan)
plot(zr2~Ns,xlab="N(SNPs)",ylab="prediction R2")
lines(zr2~Ns)

best=which(zr2==max(zr2))
allpreds=do.call(rbind,zr[[best]])
	
getPal = colorRampPalette(brewer.pal(9, "Set1"))
repColors= getPal(length(reps))
N=Ns[best]
jitter=0.01*max(allpreds$true)
allpreds$re.lm=rescale(allpreds$lm,range(allpreds$true))+rnorm(nrow(allpreds),0,jitter)
allpreds$re.rr=rescale(allpreds$rr,range(allpreds$true))+rnorm(nrow(allpreds),0,jitter)
allpreds$re.true=allpreds$true+rnorm(nrow(allpreds),0,jitter)

plot(re.lm~re.true,allpreds,main=paste("z:",zscan[best]," Nsnps:",N," simple betas"),col= repColors,pch=16,ylab="predicted",xlab="observed")
mtext(round(cor(allpreds$lm,allpreds$true)^2,2))
plot(re.rr~re.true,allpreds,main=paste("z:",zscan[best]," Nsnps:",N," rr betas"),col= repColors,pch=16,ylab="predicted",xlab="observed")
mtext(round(cor(allpreds$rr,allpreds$true)^2,2))
plot(lm~rr,allpreds,col= repColors,main="prediction comparison",pch=16)
dev.off()

