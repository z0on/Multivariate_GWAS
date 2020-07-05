  traits = "pd.traits"
  bams = "bams.qc"
  ibs="zz8.ibsMat"

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

goods=intersect(goods.ibs,goods.traits)
ibs=ibs[goods,goods]
dim(ibs)

library(vegan)

nsamp=10;r=1
for (r in 1:100) {
	test=sample(goods,nsamp)
	goods.use=goods[!(goods %in% test)]
	ibs.r=ibs[goods.use,goods.use]
	pp0=capscale(ibs.r~1)
	mds=cbind(sample=row.names(pp0$CA$u),pp0$CA$u[,1:7])
	write.table(mds,file=paste("mds7_",r,"_",nsamp,sep=""),quote=F,row.names=F)
	mds=cbind(sample=row.names(pp0$CA$u),pp0$CA$u[,1:2])
	write.table(mds,file=paste("mds2_",r,"_",nsamp,sep=""),quote=F,row.names=F)
	writeLines(test,paste("rep",r,"_",nsamp,sep=""))
}

nsamp=25;r=1
for (r in 1:100) {
	test=sample(goods,nsamp)
	goods.use=goods[!(goods %in% test)]
	ibs.r=ibs[goods.use,goods.use]
	pp0=capscale(ibs.r~1)
	mds=cbind(sample=row.names(pp0$CA$u),pp0$CA$u[,1:7])
	write.table(mds,file=paste("mds7_",r,"_",nsamp,sep=""),quote=F,row.names=F)
	mds=cbind(sample=row.names(pp0$CA$u),pp0$CA$u[,1:2])
	write.table(mds,file=paste("mds2_",r,"_",nsamp,sep=""),quote=F,row.names=F)
	writeLines(test,paste("rep",r,"_",nsamp,sep=""))
}