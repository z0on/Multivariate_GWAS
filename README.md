# Multivariate GWAS 
## based on constrained ordination

### Advantages:
- leverages existing SNP covariance structure to detect polygenic signals
- any number of covariates can be removed without loss of power
- naturally generates empirical null distribution to detect true signal and calculate p-values
- several correlated traits can be used together as a "compound trait"

The key script here is **RDA_GWAS.R**. Below is the "help" page it would print if run without any arguments. The genotype file needed to run these, *chr14.postAlleles.gz*, is here: https://www.dropbox.com/s/12oi4dmfep7meup/chr14.postAlleles.gz . 

This project is based on the idea of using constrained ordination to look for genotype-environment associations, presented in papers by Brenna R. Forester et al: 
https://doi.org/10.1111/mec.13476
https://doi.org/10.1111/mec.14584

#### RDA_GWAS.R: Arguments (things we need to run this method)
> **Note:** all tables must be space-delimited, and can be compressed .gz files.

**gt=[filename]** Genotypes: table of minor allele counts (rows - loci, columns - samples) the first two columns must be chromosome, position header line must be present (chr, pos, sample names). Better have multiple genotype tables, split by chromosome, to run things with less memory and in parallel.

**covars=[filename]**  Table of covariates to use (rows - samples, columns - covariates). First column must be sample names. Header line must be present (sample, covariates). May not fully match the genotype table. Rows containing NA will be removed.

**traits=[filename]** Table of trait(s). First column must be sample names. There must be at least 2 columns (samples, 1 trait). Header line must be present (sample, traits). May not fully match the genotype table. Rows containing NAs will be removed.

**gdist=[filename]** Matrix of genetic distances between samples listed in the genotype file (e.g. IBS matrix from angsd). No header line or other non-numeric columns.

**gdist.samples=[filename]** Single-column list of sample names *exactly corresponding* to the genotype AND genetic distances matrix. Could be filenames with leading path and trailing extension (these will be removed) - basically use the same file that was used for -b argument in angsd to obtain IBS matrix and genotypes.

**hold.out=[filename]**  File listing sample names to hold out from the whole analysis for subsequent testing. Omit this if working with all samples

#### other RDA_GWAS.R arguments

**outfile=[filename]**  Output file name.

**plots=TRUE** Whether to plot fun plots ([outfile]_plots.pdf).

**nsites=5500000** Number of sites to compute FDR (for Manhattan plot).

**prune.dist=50000** Pruning distance (chosen SNPs must be at least that far apart).

**Output or RDA_GWAS.R:**   RData bundle containing results table for pruned SNPs (*out*) with zscores, pvalues, betas and R2s, their genotypes (*gt.s*), genotypes of the test sample set (*gt.test*), sample scores for the trait (*sample.scores* - not including the hold-out samples), and Manhattan plot data for all sites (*manh*).

## Simple run, for a whole dataset (without hold-out samples) ## 
 
```bash
# Assuming we have multiple *.postAlleles.gz files with genotypes, one file per chromosome

>allchroms
for CHR in `ls *postAlleles.gz`; do
# remove everything from the file name after the first non-alphanumeric character, to get neater output names
OUTN=`echo $CHR | perl -pe 's/^([\w\d]+)\\..+/$1/'`;
echo "Rscript RDA_GWAS.R gt=$CHR covars=mds2 traits=pd.traits gdist.samples=bams.qc gdist=zz8.ibsMat outfile=${OUTN}_pd.RData">>allchroms;
done

# examine per-chomosome *pdf files - they would contain sample and SNP ordination plots, q-q plot to show if there is any signal, manhattan plots of all and chosen SNPs (after distance-pruning)
# execute all commands in allchroms

# to compile all chromosomes together and plot genome-wide manhattan plot (only uses contigs with "chr" in the name!):
# list all per-chromosome output files in a text file:
ls *_pd.RData >pds
Rscript compile_chromosomes.R in=pds

```
Running this with hold-out samples is a bit more involved. 

Mikhail Matz, matz@utexas.edu, July 2020

