require(ggplot2 )

infile=commandArgs(T)
#infile="bleachfull.RData"
load(infile)
pdf(paste(infile,"_manhattan.pdf",sep=""),width=15,height=3.5)
manh$pos.Mb=manh$pos/1e+3
mh=manh[abs(manh$zscore)>3,]
mh=mh[grep("chr",mh$chrom),]
mh$chrN=gsub("\\D","",mh$chrom)
mh$chrom=factor(mh$chrom,levels=unique(mh$chrom)[order(as.numeric(unique(mh$chrN)))])
sign=as.numeric(mh$zscore>0)
sign[sign==0]=-1
mh$slogp=mh$logp*sign
mh$logp.adj[mh$logp.adj>2]=2
mh$lpa2=(mh$logp.adj)^3
pp=ggplot(mh,aes(pos.Mb,abs(zscore)))+
	geom_point(shape = 21, colour = "grey20", aes(size=lpa2,fill=slogp))+
    scale_size_continuous(limits=c(0,8),breaks=c(0.3,0.6,1,1.3,2)^3,labels=c(0.5,0.25,0.1,0.05,1e-2))+
  #  scale_size_continuous(limits=c(0,3),breaks=c(1,1.3,2),labels=c(0.1,0.05,1e-2))+
	scale_fill_gradient(low="cyan3",high="coral")+
	theme_bw() + labs(size = "p.adj")+theme(axis.text.x=element_text(angle=45, hjust=1))+
	ylim(min(abs(mh$zscore)),max(c(10,max(abs(mh$zscore)))))+
	ggtitle(paste(infile,"z > 3"))+
	xlab("position,Mb")+facet_grid(~chrom,scale="free_x",space="free_x")
plot(pp)
dev.off()
