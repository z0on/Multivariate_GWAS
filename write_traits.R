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


load("rda_covariates.RData")
pp0=capscale(ibs~1)
mds=pp0$CA$u[,1:7]
mds=cbind(sample=row.names(mds),mds)
head(mds)
write.table(mds,file="mds7",quote=F,row.names=F)


covars=read.table("~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS/simple.covars",header=T,stringsAsFactors=F,sep="\t")
write.table(covars,file="~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS/simple.covars",quote=F,row.names=F)
