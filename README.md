# Multivariate GWAS 
## based on constrained ordination

### Advantages:
- leverages natural SNP covariance structure to detect polygenic signals;
- any number of nuisance covariates can be removed without loss of power;
- naturally generates empirical null distribution to detect true signal and calculate p-values;
- several correlated traits can be used together as a "compound trait".

The key script here is **`RDA_GWAS.R`**, which is designed for command-line usage (`Rscript RDA_GWAS.R [arguments]`). 
The genotype file needed to run exmple code below, `chr14.postAlleles.gz`, is here: https://www.dropbox.com/s/12oi4dmfep7meup/chr14.postAlleles.gz . 

### *RDA_GWAS.R*: Main arguments (things we need to run this method)
> **Note:** all tables must be space-delimited, and can be compressed .gz files.

`gt=[filename]` Genotypes: table of minor allele counts in each sample (rows - loci, columns - samples). The first two columns must be chromosome, position. Header line must be present (chr, pos, sample names). I recommend running the method on `gt` files for individual chromosomes, to use less memory and to run it in parallel. (see **Appendix** about how to get this from **`angsd`**)

`covars=[filename]`  Table of covariates (rows - samples, columns - covariates). First column must be sample names. Header line must be present (sample, names of covariates). May not fully match the genotype table - the script will match them using the `sample` column. Rows containing NA will be removed.

`traits=[filename]` Table of trait(s). First column must be sample names. There must be at least 2 columns (samples, 1 trait). Header line must be present (sample, names of traits). Just as *covars*, this table may not fully match the genotype table; rows containing NAs will be removed.

`gdist=[filename]` Matrix of genetic distances between samples listed in the genotype file (e.g. IBS matrix from **`angsd`**, see **Appendix**). Note: there must be no header line or other non-numeric columns.

`gdist.samples=[filename]` Single-column list of sample names *exactly corresponding* to the genotype AND genetic distances matrix. Could be filenames with leading path and trailing extension (these will be removed) - basically use the same file that was used for `-b` argument in **`angsd`** to obtain IBS matrix and genotypes (see **Appendix**).

`hold.out=[filename]`  File listing sample names to hold out from the whole analysis for subsequent testing of the polygenic score's prediction accuracy. May be omitted.

### Other *RDA_GWAS.R* arguments

`outfile=[filename]`  Output file name.

`plots=TRUE` Whether to plot diagnostic plots (`[outfile]_plots.pdf`).

`nsites=5500000` Total number of sites *acros the whole genome* that are being analyzed - this is to compute genome-wide FDR (Benjamini-Hochberg method).

`prune.dist=50000` Pruning distance (chosen SNPs must be at least that far apart).

## Simple run, for a whole dataset (without hold-out samples) ## 
Assuming we have multiple `*.postAlleles.gz` files with genotypes, one file per chromosome:
```bash
>allchroms
for CHR in `ls *postAlleles.gz`; do
# removing everything from the file name after the first non-alphanumeric character, to get neater output names
OUTN=`echo $CHR | perl -pe 's/^([\w\d]+)\\..+/$1/'`;
# writing down a list of calls to RDA_GWAS.R, for each chromosome
echo "Rscript RDA_GWAS.R gt=$CHR covars=mds2 traits=pd.traits gdist.samples=bams.qc gdist=zz8.ibsMat outfile=${OUTN}_pd.RData">>allchroms;
done
```
Execute all commands in `allchroms` (preferably in parallel)

This will generate RData bundles, one for each chromosome, containing the following R objects:
* `out` : results table for pruned SNPs containing zscores, pvalues, betas (for simple linear model and elastic net regression), and r-squares for lm regressions;
* `gt.s` : genotypes of chosen SNPs;
* `gt.test` : genotypes of hold-out samples (if any) at the selected SNPs;
* `manh` : manhattan plot data (zscores, pvalues) for ALL analyzed sites;
* `sample.scores` : sample scores along the first constrained ordination axis.

Also, unless `plots=FALSE` option was given, there will be `[outfile]_plots.pdf` files generated for each chromosome, containing the following plots:

![sample ordination](sample_ordination.png)
* constrained ordination plot for samples, and the trait(s) vector(s). The analysis uses sample scores along the first constrained axis, CAP1, but multiple correlated traits can be used to define it. 
> Note: the metod always flips the CAP1 axis so that increase in the trait value (specifically, first column in the *traits* table) corresponds with increase of CAP1 score. Since CAP1 orientation is arbitrary, this does not change anything except making the results easier for a human to comprehend. The plot above is a raw plot, before flipping - the SNP ordination plot below (colored dots in rings) will be a flipped version of this one. 

![qq plot](qqplot.png)
* q-q plot of SNP scores along CAP1 compared to SNP scores along a very high-order MDS representing noise. Departure upwards from the red line at the top right corner indicates positive signal, departure downwards in the lower left corner - negative signal. In this case these is definitely some positive signal, but little or no negative signal.

![snp scores](snp_ordination.png)
* SNP scores in the same ordination space: CAP1 (trait) vs MDS100 (noise). Colored rings are increasing z-scores of distance from 0, the outmost ring is z > 5. The idea is to check if the cloud is more extended / has more outliers along CAP1 compared to MDS100.

![raw manhattan](raw_manhattan.png)
* Manhattan plot of all analyzed sites. Adjusted p-values are supposed to be genome-wide, if the total number of analyzed SNPs (across the whole genome) was supplied to *RDA_GWAS.R* as *nsites=1234567* argument.

![pruned manhattan](pruned_manhattan.png)
* Manhattan plot for distace-pruned top-zscore SNPs. Pruning follows the same procedure as LD-clumping but with physical distances instead of LD (LD stuff is currently in the works). In short, out of SNPs with z-score exceeding 2, we choose the best-zscore SNP, remove all SNPs within `prune.dist` from it that have the same z-score sign, repeat for next-best remaining SNP, and so on. `prune.dist` is the argument to `RDA_GWAS.R`, default is 50000.

![zscan](back_predict_zscan.png)  
* Change in trait prediction accuracy when relaxing z-score cutoff (from just the few top-zscore SNPs to a whole lot of SNPs including those with much worse zscore)  

![bp](back_predict_simple_betas.png)  
* Predicted vs observed trait, using prediction from simple lm model.  

![bp](back_predict_rr_betas.png)
* Predicted vs observed trait, using prediction from regularized regression model.  
>**NOTE:** In this run, the last three plots are more like a sanity check rather than actual test of prediction power, because the trait is predicted **in the same samples** that were used to derive predictions. This always works too well! To check our prediction power properly, we must predict the trait in "hold-out" samples, that were NOT used for deriving predictions (see next section).

To combine all chromosomes together and plot genome-wide manhattan plot (the plot only uses contigs with "chr" in the name!):
```bash
ls *_pd.RData >pds
Rscript compile_chromosomes.R in=pds
```
## Run with hold-out samples ##
The idea is to withold some samples from the analysis and use them later to test whether we can predict the trait in them from their genotypes. This is a bit more involved. First, we need to list hold-out sample names in a file. We might wish to make many such files listing randomly picked hold-out samples. So we will have a bunch of sample-listing files named, for example, `rep10_25` - which would be 10th replicate of witholding 25 samples. We might also need to make replicate-specific tables of covariates, especially if they include unconstrained MDSes - those we are supposed to compute based on the dataset *without the hold-out samples*. See 'write_holdout_reps.R' for example R code (*spaghetti warning...*). 

Then we basically need to run the above code for each of the replicates, with additional `hold.out=[filename]` argument to `RDA_GWAS.R` . Here is one way to do this with bash looping. Note that in this case we are setting up a run with 50 hold-out replicates, with replicate-specific covariate files named like `mds2_10_25` (to correspond with the hold-out samples filenames like `rep10_25`). Note that we do NOT need to subset our genotypes, genetic distances, or traits tables - this will happen automatically, just supply full files for all samples. Also there is no need to have replicate-specific covariates if the covariates are not computed from the dataset itself (for example, sequencing batch, population designations of each sample, sequencing quality metrics).

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
Executing all commands in `pdd` (much preferrably in parallel!) will give us 50 hold-out runs per each chromosome. To sort out this mess, first we need to compile results for all chromosomes for each hold-out replicate, using `compile_chromosomes.R`:
```bash
>compd
for r in `seq 1 50`; do
ls *pd_rep${r}_*.RData >rep${r}_pd; 
echo "Rscript compile_chromosomes.R in=rep${r}_pd plotManhattan=FALSE">>compd;
done
```
Executing all commands in `compd` gives us per-replicate, whole-genome results. Each output file includes genotypes for hold-out samples, not included in the main analysis. The last stage is to see how well the trait values in these samples are predicted based on their polygenic scores, and summarize the results of all replicates in a nice plot: 
```bash
ls rep*_pd.RData >reps50
compile_replicates.R reps=reps50 traits=pd.traits
```
Note that we need to supply argument `traits` again, pointing to the same table as we were using for `RDA_GWAS.R`

This will generate a plot `reps50.pdf` looking somewhat like this:

![predictions](pd_predictions.png)

Where:
* panel 1: scan through z-score cutoffs for best predictions (using lm betas). Basically it shows whether preriction improves if we keep including more and more SNPs with worse z-scores. 
* panel 2: predictions for hold-out samples based on lm betas, for the z-score cutoff giving maximal prediction accuracy (listed above the panel). Colors are different replicates, just for sanity check: all the well-predicted samples must not belong to a single replicate. 
* panel 3: same as panel 2, but based on regularized betas (from elastic net regression)
* panel 4: comparison of simple and regularized predictions

## Where does it come from?
This project is based on the idea of using constrained ordination to look for genotype-environment associations, presented in papers by Brenna R. Forester et al: 
https://doi.org/10.1111/mec.13476 ; https://doi.org/10.1111/mec.14584


# Appendix
#### How to get genotypes (posterior minor allele counts) and genetic distance matrix (IBS) from ANGSD

Assume we have a file `bams.qc` listing our (indexed) bam files, from which we have already tossed all the samples that are severely under-sequenced, clonal, wrong species, or just look weird a PCoA plot. We have already decided on the genotyping rate cutoff (`-minInd` argument to `angsd`), which is the number of individuals in which a locus must be represented by at least one read (idealy it shoudl be set to 75-80% of total number of samples). We are going after variants of minor allele 0.05 and higher (`-minMaf 0.05`): 
```bash
FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 30 -minQ 20 -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 152 -snp_pval 1e-5 -minMaf 0.05 "
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doGeno 8 -doPost 1"
angsd -b bams.qc -GL 1 $FILTERS $TODO -P 12 -out zz8
```
>If this runs out of memory, try reducing `-P` , all the way to `-P 1`.

The output file `zz8.ibsMat` is the genetic distances matrix (identity-by-state) that we can use for GWAS here. The file `zz8.geno.gz` contains posterior genotype probabilities, it needs to be massaged a bit before we can use it.

First, let's unarchive it and split by chromosome:
```bash
zcat zz8.geno.gz | awk -F, 'BEGIN { FS = "\t" } ; {print > $1".split.geno"}'
```
If you have some short contigs in addition to chromosomes, you might wish to concatenate them together into a separate `unplaced.split.geno` file before proceeding. If your genome is highly fragmented, pre-concatenate it into "fake chromosomes" before mapping (see [`concat_fasta.pl`](https://github.com/z0on/2bRAD_denovo/blob/master/concatFasta.pl).

Then, we calculate posterior number of minor alleles:
```bash
>bychrom
for GF in *.split.geno; do
echo "awk '{ printf \$1\"\\t\"\$2; for(i=4; i<=NF-1; i=i+3) { i2=i+1; printf \"\\t\"\$i+2*\$i2} ; printf \"\\n\";}' $GF > ${GF/.split.geno/}.postAlleles" >>bychrom
done
```
Execute all lines in `bychrom`, and we got ourselves genotype data tables suitable for `RDA_GWAS.R`.

You might wish to compress them, for tidyness, although it is not necessary for `RDA_GWAS.R`:
```bash
for F in *.postAlleles;do gzip -f $F;done
```




Mikhail Matz, matz@utexas.edu, July 2020

