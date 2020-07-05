# Multivariate GWAS 
## based on constrained ordination

### Advantages:
- leverages existing SNP covariance structure to detect polygenic signals
- any number of covariates can be removed without loss of power
- naturally generates empirical null distribution to detect true signal and calculate p-values
- several correlated traits can be used together as a "compound trait"

The key script here is **RDA_GWAS.R**. Below is the "help" page it would print if run without any arguments. The genotype file needed to run these, *chr14.postAlleles.gz*, is here: https://www.dropbox.com/s/12oi4dmfep7meup/chr14.postAlleles.gz . See *example_calls_RDA_GWAS.sh* for test calls, or just copy-paste from here:
```bash
# for bleaching traits (3 traits together)
Rscript RDA_GWAS.R gt=chr14.postAlleles covars=mds2_10_25 traits=bleach.traits gdist.samples=bams.qc gdist=zz8.ibsMat hold.out=rep10_25 outfile=chr14.bl.mds2.rep10.RData
# for proportion of D symbionts
Rscript RDA_GWAS.R gt=chr14.postAlleles.gz covars=mds2_10_25 traits=pd.traits gdist.samples=bams.qc gdist=zz8.ibsMat hold.out=rep10_25 outfile=chr14.pd.mds2.rep10.RData
# prop D with different covariates
Rscript RDA_GWAS.R gt=chr14.postAlleles.gz covars=technical.covars traits=pd.traits gdist.samples=bams.qc gdist=zz8.ibsMat hold.out=rep10_25 outfile=chr14.pd.techcov.rep10.RData
```
This project is based on the idea of using constrained ordination to look for genotype-environment associations, presented in papers by Brenna R. Forester et al: 
https://doi.org/10.1111/mec.13476
https://doi.org/10.1111/mec.14584

#### RDA_GWAS.R: Arguments
> **Note:** all tables must be space-delimited, and can be compressed .gz files.

**gt=[filename]** Genotypes: table of minor allele counts (rows - loci, columns - samples) the first two columns must be chromosome, position header line must be present (chr, pos, sample names).

**covars=[filename]**  Table of covariates to use (rows - samples, columns - covariates). First column must be sample names. Header line must be present (sample, covariates). May not fully match the genotype table. Rows containing NA will be removed.

**traits=[filename]** Table of trait(s). First column must be sample names. There must be at least 2 columns (samples, 1 trait). Header line must be present (sample, traits). May not fully match the genotype table. Rows containing NAs will be removed.

**gdist=[filename]** Matrix of genetic distances between samples listed in the genotype file (e.g. IBS matrix from angsd). No header line or other non-numeric columns.

**gdist.samples=[filename]** Single-column list of sample names corresponding to the genotype AND genetic distance matrix. Could be filenames with leading path and trailing extension (these will be removed).

**hold.out=[filename]**  List of samples to hold out from the whole analysis for testing the predictors.

**outfile=[filename]**  Output file name.

**plots=TRUE** Whether to plot fun plots ([outfile]_plots.pdf).

**nsites=5500000** Number of sites to compute FDR (for Manhattan plot).

**prune.dist=50000** Pruning distance (chosen SNPs must be at least that far apart).

**Output:**   RData bundle containing results table for pruned SNPs (*out*) with zscores, pvalues, betas and R2s, their genotypes (*gt.s*), and Manhattan plot data for all sites (*manh*).

Mikhail Matz, matz@utexas.edu, July 2020

