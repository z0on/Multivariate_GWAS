setwd('~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS')
load('~/Dropbox/amil_RDA_association_jun2020/Amil2020_bigtraits.RData')
head(bigtraits)


propsym=read.table("~/Dropbox/amil_RDA_association_jun2020/prop_symb.plink.txt",header=T)
head(propsym)
pd=propsym[,c(1,5)]
names(pd)=c("sample","propD")
write.table(pd,file="pd.traits",quote=F,row.names=F)

bleach=bigtraits[,c(2,6,4,5)]
head(bleach)
write.table(bleach,file="bleach.traits",quote=F,row.names=F)
bleach=bigtraits[,c(2,6)]
head(bleach)
write.table(bleach,file="blscore.traits",quote=F,row.names=F)


load("~/Dropbox/amil_RDA_association_jun2020/rda_covariates.RData")
pp0=capscale(ibs~1)
mds=pp0$CA$u[,1:2]
mds=cbind(sample=row.names(mds),mds)
head(mds)
write.table(mds,file="mds2",quote=F,row.names=F)

pdmds2=merge(pd,mds,by="sample")
write.table(pdmds2,file="pdmds2",quote=F,row.names=F)


covars=read.table("~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS/simple.covars",header=T,stringsAsFactors=F)
str(covars)
write.table(covars[,c(1,3,4)],file="~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS/technical.covars",quote=F,row.names=F)
write.table(covars[,c(1,2)],file="~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS/reefsites",quote=F,row.names=F)
