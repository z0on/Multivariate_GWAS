setwd("~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS/")
reps=list()
for (r in 1:25) {
	load(paste("gws_rep",r,"_25_predictions.RData",sep=""))
	reps[[r]]=allpreds
}
Reps=do.call(rbind,reps)
head(Reps)
par(mfrow=c(1,2))
plot(lm~true,Reps,main="LM")
abline(h=0,lty=3,col="red")
abline(v=0,lty=3,col="red")
mtext(paste("R2 =",round(cor(Reps$lm,Reps$true)^2,2)))

# sign prediction accuracy
tt=table((Reps$lm>0.2 & Reps$true>0.2)| (Reps$lm<(-0.2) & Reps$true<(-0.2)))
acc=tt[2]/sum(Reps$true>0.2 | Reps$true<(-0.2))


plot(rr~true,Reps,main="Elastic Net")
abline(h=0,lty=3,col="red")
abline(v=0,lty=3,col="red")
mtext(paste("R2 =",round(cor(Reps$rr,Reps$true)^2,2)))
