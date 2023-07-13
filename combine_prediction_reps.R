setwd("~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS/")
reps=list()
for (r in 1:25) {
	load(paste("gws_rep",r,"_25_predictions.RData",sep=""))
	reps[[r]]=allpreds
}
Reps=do.call(rbind,reps)

# sign prediction accuracy
tt=table((Reps$rr>0 & Reps$true>0)| (Reps$rr<0 & Reps$true<0))
acc.rr=round(tt[2]/sum(Reps$true>0 | Reps$true<0),2)
tt=table((Reps$lm>0 & Reps$true>0)| (Reps$lm<0 & Reps$true<0))
acc.lm=round(tt[2]/sum(Reps$true>0 | Reps$true<0),2)

Reps$j.true=Reps$true+rnorm(nrow(Reps),0,sd(Reps$true)/5)
Reps$j.lm=Reps$lm+rnorm(nrow(Reps),0,sd(Reps$lm)/10)
pdf("overall_prediction_accuracy.pdf",height=4,width=6.7)
par(mfrow=c(1,2))
plot(j.lm~j.true,Reps,main="LM",xlab="true",ylab="predicted",mgp=c(2.3,1,0))
abline(h=0,lty=3,col="red")
abline(v=0,lty=3,col="red")
abline(lm(lm~true,Reps))
mtext(paste("R2 =",round(cor(Reps$lm,Reps$true)^2,2),"; class acc.=",acc.lm))

# LM sign prediction accuracy
tt=table((Reps$lm>0 & Reps$true>0)| (Reps$lm<0 & Reps$true<0))
acc.lm=tt[2]/sum(Reps$true>0 | Reps$true<0)

Reps$j.rr=Reps$rr+rnorm(nrow(Reps),0,sd(Reps$rr)/10)
plot(j.rr~j.true,Reps,main="Elastic Net",xlab="true",ylab="predicted",mgp=c(2.3,1,0))
abline(h=0,lty=3,col="red")
abline(v=0,lty=3,col="red")
abline(lm(rr~true,Reps))
mtext(paste("R2 =",round(cor(Reps$rr,Reps$true)^2,2),"; class acc.=",acc.rr))
dev.off()
