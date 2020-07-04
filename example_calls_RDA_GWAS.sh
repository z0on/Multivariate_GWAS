# for bleaching traits (3 traits together)
Rscript RDA_GWAS.R gt=chr7.postAlleles covars=mds7 traits=bleach.traits gdist.samples=bams.qc gdist=zz8.ibsMat hold.out=rep10_25 outfile=chr7.bl.mds7.rep10.RData
# for proportion of D symbionts
Rscript RDA_GWAS.R gt=chr7.postAlleles covars=mds7 traits=pd.traits gdist.samples=bams.qc gdist=zz8.ibsMat hold.out=rep10_25 outfile=chr7.pd.mds7.rep10.RData
# prop D with different covariates
Rscript RDA_GWAS.R gt=chr7.postAlleles covars=simple.covars traits=pd.traits gdist.samples=bams.qc gdist=zz8.ibsMat hold.out=rep10_25 outfile=chr7.pd.simcov.rep10.RData
