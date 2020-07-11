# Multivariate GWAS 
## based on constrained ordination

### Advantages:
- leverages existing SNP covariance structure to detect polygenic signals
- any number of covariates can be removed without loss of power
- naturally generates empirical null distribution to detect true signal and calculate p-values
- several correlated traits can be used together as a "compound trait"

The key script here is **RDA_GWAS.R**. The genotype file needed to run the examples, *chr14.postAlleles.gz*, is here: https://www.dropbox.com/s/12oi4dmfep7meup/chr14.postAlleles.gz . 

#### RDA_GWAS.R: Arguments (things we need to run this method)
> **Note:** all tables must be space-delimited, and can be compressed .gz files.

**gt=[filename]** Genotypes: table of minor allele counts (rows - loci, columns - samples) the first two columns must be chromosome, position header line must be present (chr, pos, sample names). Better have multiple genotype tables, split by chromosome, to run things with less memory and in parallel.

**covars=[filename]**  Table of covariates to use (rows - samples, columns - covariates). First column must be sample names. Header line must be present (sample, covariates). May not fully match the genotype table. Rows containing NA will be removed.

**traits=[filename]** Table of trait(s). First column must be sample names. There must be at least 2 columns (samples, 1 trait). Header line must be present (sample, traits). May not fully match the genotype table. Rows containing NAs will be removed.

**gdist=[filename]** Matrix of genetic distances between samples listed in the genotype file (e.g. IBS matrix from angsd). No header line or other non-numeric columns.

**gdist.samples=[filename]** Single-column list of sample names *exactly corresponding* to the genotype AND genetic distances matrix. Could be filenames with leading path and trailing extension (these will be removed) - basically use the same file that was used for -b argument in angsd to obtain IBS matrix and genotypes.

**hold.out=[filename]**  File listing sample names to hold out from the whole analysis for subsequent testing. Omit this if working with all samples

### other RDA_GWAS.R arguments

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

# execute all commands in allchroms

# examine per-chomosome *pdf files - they would contain:
#   - sample and SNP ordination plots, 
#   - q-q plot to show if there is any signal, 
#   -manhattan plots of all and chosen SNPs (after distance-pruning)


# to compile all chromosomes together and plot genome-wide manhattan plot (only uses contigs with "chr" in the name!):
ls *_pd.RData >pds
Rscript compile_chromosomes.R in=pds

```
## Run with hold-out samples ##
This is a bit more involved. First, we need to list hold-out sample names in a file. We might wish to make many such files listing randomly picked hold-out samples. So we will have a bunch of sample-listing files named, for example, *rep10_25* - which would be 10th replicate of witholding 25 samples. We might also need to make replicate-specific tables of covariates, especially if they include unconstrained MDSes - those we are supposed to compute based on the dataset *without the hold-out samples*. See *write_holdout_reps.R* for example R code (spagetty warning). 

Then we basically need to run the above code for each of the replicates, with additional hold.out=[filename] argument to *RDA_GWAS.R* . Here is one way to do this with bash looping. Note that in this case we are setting up a run with 50 hold-out replicates, with replicate-specific covariate files named like *mds2_10_25* (to correspond with the hold-out samples filenames like *rep10_25*). Note that we do NOT need to subset our genotypes, genetic distances, or traits tables - this will happen automatically, just supply full files for all samples.

```bash
>pdd
for R in `seq 1 50`; do
REP=rep${R}_10;
MDS=mds2_${R}_10;
for CHR in `ls *postAlleles.gz`; do
OUTN=`echo $CHR | perl -pe 's/.+(chr\d+).+/$1/'`;
echo "Rscript RDA_GWAS.R gt=$CHR covars=$MDS traits=pd.traits gdist.samples=bams.qc plots=FALSE gdist=zz8.ibsMat hold.out=$REP outfile=${OUTN}_pd_${REP}.RData">>pdd;
done;
done
```
Executing all commands in *pdd* (much preferrably in parallel!) will give us 50 hold-out runs per each chromosome. To sort out this mess, first we need to compile results for all chromosomes for each hold-out replicate, using *compile_chromosomes.R*:
```bash
>compd
for r in `seq 1 50`; do
ls *pd_rep${r}_*.RData >rep${r}_pd; 
echo "Rscript compile_chromosomes.R in=rep${r}_pd plotManhattan=FALSE">>compd;
done
```
Executing all commands in *compd* gives us per-replicate, whole-genome results. Each output file includes genotypes for hold-out samples, not included in the main analysis. The last stage is to see how well the trait values in these samples are predicted based on their polygenic scores, and summarize the results of all replicates in a nice plot: 
```bash
ls rep*_pd.RData >reps50
compile_replicates.R in=reps50 traits=pd.traits
```
Note that we need to supply argument *traits* again, pointing to the same table as we were using for *RDA_GWAS.R*

This will generate a plot *reps50.pdf* looking somewhat like this:

![predictions](pd_predictions.png)

Where:
* panel 1: scan through z-score cutoffs for best predictions (using lm betas)
* panel 2: predictions for hold-out samples based on lm betas
* panel 3: predictions based on regularized betas (glmnet)
* panel 4: comparison of simple and regularized predictions

### Where does it come from?
This project is based on the idea of using constrained ordination to look for genotype-environment associations, presented in papers by Brenna R. Forester et al: 
https://doi.org/10.1111/mec.13476
https://doi.org/10.1111/mec.14584


Mikhail Matz, matz@utexas.edu, July 2020

