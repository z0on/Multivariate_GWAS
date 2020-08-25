require(fields)
require(zoo)
require(data.table)
require(R.utils)

#setwd("~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS")
infile=commandArgs(T)

wd=10;threshold=0.1; step=5; r2limit=0.6

message("reading input...",appendLF=FALSE)
# infile="bighead.gz"
ld=fread(infile,nThread=4)
ld=ld[,-c(4:6)]
message("done")
names(ld)=c("pos1","pos2","dis","r2em")
setkey(ld,pos1,pos2)
all1=unique(ld$pos1)

#snps=unique(c(ld$pos1,ld$pos2))
snps=unique(ld$pos1)
pos=as.numeric(as.character(sub(".+:","",snps)))
snpso=snps[order(pos)]


# x=snpso[100:109];LD=ld
# start=which(snpso=="chr1:56")
# x=snpso[start:(start+9)]
ld.drop=function(x,LD=ld,th=0.1,qs=log(c(1,50,100,200,400,800,1600,3200,6400,10000,15000,20000,25000)), plots=F,r2=0.6) {
	lds1=LD[x]
	lds2=LD[pos2 %in% x]
	lds=rbind(lds1,lds2)
    message(x[1]," N:",nrow(lds))
#	setkey(lds,dis)
	lds$quant=cut(log(lds$dis+1), qs) 
#	table(lds$quant)
	qr2=aggregate(lds$r2em,list(lds$quant),mean)
	qd=aggregate(lds$dis,list(lds$quant),mean)
#	plot(qr2$x~qd$x)
	xx=log(qr2$x+1e-3)
	yy=log(qd$x)
	lo=lm(yy~xx) #,span=1.1)
	if(summary(lo)$adj.r.squared<0.6 | coef(lo)[2]>0) { return (NA) }
	log.edge=predict(lo,newdata=data.frame(xx=log(threshold)))
	if(plots) {
			plot(yy~xx,main=paste(x[1],round(exp(log.edge),0)))
			xs=seq(min(xx),max(xx),length=100)
			lop=predict(lo,newdata=data.frame(xx=xs))
			lines(xs,lop,col="purple",lwd=3)
			abline(h=log.edge,col="grey70")
			abline(v=log(threshold),col="grey70")
	}
	return(exp(log.edge))
}

meanpos=function(x) {
	pos=as.numeric(as.character(sub(".+:","",x)))
	return(round(mean(pos),0))
}

# cuts=c(1,100);f=1.5
# for (i in 1:100) {
	# nn=cuts[length(cuts)]*f
	# if(nn>25000) {break}
	# cuts=append(cuts,nn)
# }
# cuts=append(cuts,25000)
cuts=c(1,50,100,200,400,800,1600,3200,6400,10000,15000,20000,25000)
# LD$quant=cut(log(LD$dis+1), log(cuts)) 


# # #	table(lds$quant)
	# qr2=aggregate(LD$r2em,list(LD$quant),mean)
	# qd=aggregate(LD$dis,list(LD$quant),mean)
	# plot(qr2$x~qd$x)
	# xx=log(qr2$x+1e-3)
	# yy=log(qd$x)
	# lo=rlm(yy~xx) #,span=1.1)
	# ?rlm
	# log.edge=predict(lo,newdata=data.frame(xx=log(threshold)))
	# plot(yy~xx,main=paste(x[1],round(exp(log.edge),0)))
	# xs=seq(min(xx),max(xx),length=100)
	# lop=predict(lo,newdata=data.frame(xx=xs))
	# lines(xs,lop,col="purple",lwd=3)
	# abline(h=log.edge,col="grey70")
	# abline(v=log(threshold),col="grey70")
# exp(log.edge)


rolldrop=rollapply(snpso,wd,FUN=function(x) ld.drop(x,th=threshold,plots=F,r2=r2limit),by=step)
#table(rolldrop<250)
message("goods : NAs")
 table(is.na(rolldrop))

rolldrop[rolldrop<50]=50
#rolldrop[rolldrop>25000]=25000
rollpos=rollapply(snpso,wd,meanpos,by=step)
rollpos=rollpos[1:length(rolldrop)]

rolldrop=rolldrop[!(is.na(rolldrop))]
rollpos=rollpos[!(is.na(rolldrop))]

rolldrop=c(rolldrop[1],rolldrop,rolldrop[length(rolldrop)])
rollpos=c(1,rollpos,max(pos))
rdlm=data.frame(cbind(rollpos,rolldrop))

tail(rdlm)
head(rdlm)
pdf(paste(infile,".ldlm01.pdf",sep=""),width=20, height=4)
	 yy=rolldrop
	 xx=rollpos
	lo=loess(yy~xx,span=0.001)
	lop=predict(lo,newdata=data.frame(xx=rdlm$rollpos))
	lop[lop<250]=250
	plot(rdlm,type="l",lty=3,ylim=c(0,max(lop)*1.5),ylab=paste("distance ro R2=",threshold,sep=""),xlab="position")
	# lines(rdlm,type="l",lty=3,col="red")
	 points(rdlm,cex=0.5,pch=16,col="grey50")
	# points(lo,pch="*",col="red")
	lines(rdlm$rollpos,lop,col="coral",lwd=2)
	lines(rdlm$rollpos,lop*1.5,col="cyan3",lwd=2)
	legend("topleft",lty=1, lwd=2,col=c("coral","cyan3"),legend=c("loess","loess*1.5"),bty="n")
dev.off()

save(rdlm,lop,file=paste(infile,".ldlm01.RData",sep=""))
