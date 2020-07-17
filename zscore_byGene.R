if (length(commandArgs(trailingOnly=TRUE))<1) {
	options(warning.length=8000)
	stop("

Writes down max or mean z-score per gene, for GO_MWU analysis.

Arguments: 

infile=[filename]    .RData file output by RDA_GWAS.R or compile_chromosomes.R

genes=[filename]     table of gene coordinates with four columns: 
                     chromosome, start, end, gene name.
                     Gene names should be the same as used in the GO annotations table.

margin=0             Length of additional regions on both ends of gene

stat=max             Statistic to calculate: max or mean z-score per gene

Output: csv file [infile]_[stat]_2go, suitable for GO_MWU                 

Mikhail Matz, matz@utexas.edu, July 2020

") }

#To generate "genes" table from a typical gff3 file:
# grep gene mygenome.gff3 | awk '$3==\"gene\"' | cut -f 1,4,5,9 | perl -pe 's/ID=(Prefix[\d\w]+).+/$1/' >mygenome_genes.gff3
#where \"Prefix\" is the common start of all gene IDs.

gff3 =grep("genes=",commandArgs())
if (length(gff3)==0) { stop ("specify genes file (genes=filename)\n") }
gff3=sub("genes=","", commandArgs()[gff3])

infile =grep("infile=",commandArgs())
if (length(infile)==0) { stop ("specify input file, produced by RDA_GWAS.R or compile_chromosomes.R (infile=filename)\n") }
infile =sub("infile=","", commandArgs()[infile])

stat="max"
if(length(grep("stat=mean",commandArgs()))>0) { stat="mean" } 

margin =grep("margin=",commandArgs())
if(length(margin)>0) { margin=as.numeric(sub("margin=","", commandArgs()[margin])) } else { margin=0 }


# infile="pdfull.RData"
# gff3="amil_genes.gff3"
# margin=0
# stat="max"

gtab=read.table(gff3)
names(gtab)=c("chr","start","end","gname")
load(infile) # supposed to contain out (selected sites) and manh (all sites) objects 

genes2go=c();i=1;nn=c();maxz=c();gene=c()
for (i in 1:length(levels(manh$chrom))) {
	ch=levels(manh$chrom)[i]
	message("\n",ch)
	gtc=subset(gtab,chr == ch)
	if (nrow(gtc)==0) { next }
	zzc=subset(manh,chrom==ch)
	zzc$pos=zzc$pos*1000
	pb=txtProgressBar(0,nrow(gtc))
	for ( g in 1:nrow(gtc)) {
		z=subset(zzc,pos>=(gtc[g,"start"]-margin) & pos<(gtc[g,"end"]+ margin))
		nn=append(nn,nrow(z))
		if(nrow(z)==0) { maxz=append(maxz,0) } else { 
			if (stat=="max") { 	
				maxz=append(maxz,z$zscore[which(abs(z$zscore)==max(abs(z$zscore)))]) 
			} else {
				maxz=append(maxz,mean(z$zscore)) 				
			}
		 }
		 gene=append(gene,as.character(gtc$gname[g]))
		setTxtProgressBar(pb,g)
	}
}
#table(nn>=5)
maxz=maxz[nn>=5]
gene=gene[nn>=5]
nn=nn[nn>=5]

genes2go=data.frame(cbind(gene=gene,z=maxz,nn=nn),stringsAsFactors=F)
genes2go$z=as.numeric(genes2go$z)
genes2go$nn=as.numeric(genes2go$nn)
#which(is.na(genes2go$z))

# ------ linear model-based correction for log(number of SNPs per gene)
maxz=genes2go$z
nn=genes2go$nn

 lnn.p=lm(abs(maxz)~log(nn))
 lnn.n=lm(-abs(maxz)~log(nn))
 maxz.c=maxz
 maxz.c[maxz>0]=resid(lnn.p)[maxz>0]
 maxz.c[maxz<0]=resid(lnn.n)[maxz<0]

pdf(paste(infile,"_lmscaled_z_",stat,".pdf",sep=""),width=10,height=10)
 plot(maxz~log(nn))
 abline(lnn.p)
 abline(lnn.n)
 points(maxz.c~log(nn),col="red")
dev.off()

genes2go$zc=maxz.c
g2g=genes2go[,c(1,4)]
write.csv(g2g,file=paste(infile,"_lmscaled_z_",stat,".csv",sep=""),quote=F,row.names=F)

# ------ same, via quantile rescaling

qq=scales::rescale(seq(0,1,0.05),range(log(nn)))
rmz=c()
for (i in 2:length(qq)) {
	qm=genes2go[log(genes2go$nn)<=qq[i] & log(genes2go$nn)>=qq[i-1],]
	qm$z=scale(qm$z)
	rmz=rbind(rmz,qm)
}
rmz=rmz[order(rmz$gene),]
genes2go=genes2go[order(genes2go$gene),]

pdf(paste(infile,"_quantscaled_z_",stat,".pdf",sep=""),width=10,height=10)
 plot(maxz~log(nn))
 points(rmz$z~log(rmz$nn),col="red")
dev.off()
g2g=rmz[,1:2]
write.csv(g2g,file=paste(infile,"_qscaled_z_",stat,".csv",sep=""),quote=F,row.names=F)

