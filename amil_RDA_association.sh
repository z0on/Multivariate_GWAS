
# the project is archived at 
ranch.tacc.utexas.edu:/stornext/ranch_01/ranch/users/01211/cmonstr/rda_gwas_july2020
$ARCHIVER:/stornext/ranch_01/ranch/users/01211/cmonstr/rda_gwas_july2020

>cpp
for F in `ls /corral-repl/utexas/tagmap/zach_shared/*.bam`; do
echo "cp $F ." >>cpp;
done
ls5_launcher_creator.py -j cpp -n cpp -a tagmap -t 4:00:00 -w 24 -e matz@utexas.edu -q normal
sbatch cpp.slurm


# indexing huge bam files, to make csi indices
>baii
for F in `ls *.bam`; do
echo "samtools index -m 14 $F" >>baii;
done
ls5_launcher_creator.py -j baii -n baii -a tagmap -t 4:00:00 -w 4 -e matz@utexas.edu -q normal
sbatch baii.slurm


# extracting cladeC
>s3
for F in `ls *es.bam`; do
echo "samtools view $F symbiont3 -bF 0x400 -O bam -o ${F/.marked_duplicates.bam/}.sym3.bam && samtools index -m 14 ${F/.marked_duplicates.bam/}.sym3.bam" >>s3;
done
ls5_launcher_creator.py -j s3 -n s3 -a tagmap -t 2:00:00 -w 4 -e matz@utexas.edu -q normal
sbatch s3.slurm


ls *sym3.bam > s3bams
>s3i
for F in `ls *sym3.bam`; do
echo "samtools index -m 14 $F" >>s3i;
done
ls5_launcher_creator.py -j s3i -n s3i -a tagmap -t 0:30:00 -w 24 -e matz@utexas.edu
sbatch s3i.slurm

ls *sym4.bam > s4bams
>s4i
for F in `ls *sym4.bam`; do
echo "samtools index -m 14 $F" >>s4i;
done
ls5_launcher_creator.py -j s4i -n s4i -a tagmap -t 0:30:00 -w 24 -e matz@utexas.edu
sbatch s4i.slurm


FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1000 -checkBamHeaders 0"
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
angsd -b s3bams -r symbiont3:1-1000000 -GL 1 $FILTERS $TODO -P 12 -out dd3 


# extracting cladeD
>s4
for F in `ls *es.bam`; do
echo "samtools view $F symbiont4 -bF 0x400 -O bam -o ${F/.marked_duplicates.bam/}.sym4.bam && samtools index -m 14 ${F/.marked_duplicates.bam/}.sym4.bam" >>s4;
done
ls5_launcher_creator.py -j s4 -n s4 -a tagmap -t 2:00:00 -w 4 -e matz@utexas.edu -q normal
sbatch s4.slurm


# ----- subsampling high-coverage data to 20% of reads, fixing headers

cp /corral-repl/utexas/tagmap/zach_shared/Amil.symb.cat.genome.fasta* .
samtools faidx $REF

REF=Amil.symb.cat.genome.fasta
>subsamp2
for F in `ls /corral-repl/utexas/tagmap/zach_shared/*.bam`; do
outname=`echo $F | perl -pe 's/\/corral-repl\/utexas\/tagmap\/zach_shared\/|\.marked_duplicates\.bam//g'`
echo "samtools view -F 260 -s 0.2 $F | samtools view -bS -t ${REF}.fai -O bam -o ${outname}.subsampled.bam">>subsamp2;
done
ls5_launcher_creator.py -j subsamp2 -n subsamp2 -a mega2014 -t 4:00:00 -w 24 -e matz@utexas.edu -q normal
sbatch subsamp2.slurm


ls /corral-repl/utexas/tagmap/zach_shared/low_coverage/*.bam >>bams

# removing low-coverage files that are among high-coverage, also removing high-coverage files that are not among ba
# used script hc_lc_bams.R
ls *subsampled.bam
cat hcbams bams.lc >bams.combo

# running mpileup (to remove overlapping read-pairs)
REF=Amil.symb.cat.genome.fasta
echo "samtools mpileup -b bams.combo -f $REF --output-MQ >bams.qc.pileup">pile
ls5_launcher_creator.py -j pile -n pile -a mega2014 -t 12:00:00 -w 1 -e matz@utexas.edu -q normal
sbatch pile.slurm

# Calculate fasta lengths, making bed file per-chromosome
cat Amil.symb.cat.genome.fasta | awk '$0 ~ ">" {if (NR > 1) {print 0"\t"(c-1);} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print 0"\t"(c-1); }' > amil.bed

# writing mpileup commands per chromosome
REF=Amil.symb.cat.genome.fasta
>pile.chroms
for C in `seq 1 14`; do
head -$C amil.bed | tail -1 >chr${C}.bed;
echo "samtools mpileup -b bams.combo -l chr${C}.bed -f $REF --output-MQ >chr${C}.pileup">>pile.chroms
done
ls5_launcher_creator.py -j pile.chroms -n pile.chroms -a mega2014 -t 12:00:00 -w 8 -e matz@utexas.edu -q normal
sbatch pile.chroms.slurm

# for -minMaf 0.05, all good bams (incl those without traits)
export FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 152 -snp_pval 1e-5 -minMaf 0.05 "
export TODO8="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"
echo '#!/bin/bash
#SBATCH -J zz8
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p largemem512GB
#SBATCH -o zz8.o%j
#SBATCH -e zz8.e%j
#SBATCH -t 24:00:00
#SBATCH -A tagmap
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
angsd -b bams.qc -GL 1 $FILTERS $TODO8 -P 12 -out zz8' > zz8
zz8job=$(sbatch zz8 | grep "Submitted batch job" | perl -pe 's/\D//g')





FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1000 -checkBamHeaders 0"
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
angsd -b bams -r chr2:1000-1001000 -GL 1 $FILTERS $TODO -P 12 -out dd 

Rscript ~/bin/plotQC.R dd >qranks


# rerunning for -doGeno 8 (ngsLD and GS in WGCNA)
module load TACC-largemem

# for -minMaf 0.05, all good bams (incl those without traits)
export FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 152 -snp_pval 1e-5 -minMaf 0.05 "
export TODO8="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"
echo '#!/bin/bash
#SBATCH -J zz8
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p largemem512GB
#SBATCH -o zz8.o%j
#SBATCH -e zz8.e%j
#SBATCH -t 24:00:00
#SBATCH -A tagmap
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
angsd -b bams.qc -GL 1 $FILTERS $TODO8 -P 12 -out zz8' > zz8
zz8job=$(sbatch zz8 | grep "Submitted batch job" | perl -pe 's/\D//g')

export FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 152 -snp_pval 1e-5 -minMaf 0.05 "
export TODO8="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"
echo '#!/bin/bash
#SBATCH -J zz9
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p largemem512GB
#SBATCH -o zz8.o%j
#SBATCH -e zz8.e%j
#SBATCH -t 12:00:00
#SBATCH -A tagmap
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
angsd -b bams.qc -r chr1 -GL 1 $FILTERS $TODO8 -P 12 -out zz9' > zz9
zz8job=$(sbatch zz9 | grep "Submitted batch job" | perl -pe 's/\D//g')

zcat zz9.geno.gz | awk '{ printf $1"\t"$2; for(i=4; i<=NF-1; i=i+3) { i2=i+1; printf "\t"$i+2*$i2} ; printf "\n";}' | gzip > zz9.postAlleles.gz

# saving read counts for a single SNP
# chr1_10235367
export TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -dumpCounts 4"
angsd -b bams.qc -r chr1:10235367 -GL 1 $TODO -P 12 -out chr1_10235367

zcat zz9.geno.gz | head -1000 > genotest
cat zz9.postAlleles | head -1000 >patest

# to make BCF file (not needed typically, skip this)
export FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 152 -snp_pval 1e-5 -minMaf 0.05 "
export TODO1="-gl 1 -dopost 1 -domajorminor 1 -domaf 1 -dobcf 1 --ignore-RG 0 -dogeno 1 -docounts 1"
echo '#!/bin/bash
#SBATCH -J zz1
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p development
#SBATCH -o zz1.o%j
#SBATCH -e zz1.e%j
#SBATCH -t 2:00:00
#SBATCH -A tagmap
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
angsd -b bams.qc -r chr1 $FILTERS $TODO1 -P 6 -out zz1' > zz1
sbatch zz1


# splitting zz8,geno.gz by chromosome
zcat zz8.geno.gz | awk -F, 'BEGIN { FS = "\t" } ; {print > $1".split.geno"}'  
# concatenating small contigs
for i in `seq 0 9`;do
cat x*Sc0000${i}*geno >xSc0000${i}.pool.split.geno;
done
cat x*Sc0001*geno >> xSc00009.pool.split.geno
cat x*Sc0000[4378].pool.*geno >> xSc00009.pool.split.geno && rm x*Sc0000[4378].pool.*geno
rm x*Sc*[0-9].split*
rm Sc*[0-9].split*


# calculating posterior number of derived alleles
>bychrom
for GF in *.split.geno; do
echo "awk '{ printf \$1\"\\t\"\$2; for(i=4; i<=NF-1; i=i+3) { i2=i+1; printf \"\\t\"\$i+2*\$i2} ; printf \"\\n\";}' $GF > ${GF/.split.geno/}.postAlleles" >>bychrom
done


# removing extra columns
>r22
for F in *split.geno.LD; do
#echo "cat $F | awk '\$7>=0.05' | cut -f 1,2,3,7 >${F/.*/}.r2">>r22;
echo "cat $F | cut -f 1,2,3,7 >${F/.*/}.full.r2">>r22;
done
ls5_launcher_creator.py -j r22 -n r22 -a tagmap -t 0:10:00 -w 24 -e matz@utexas.edu
sbatch r22.slurm


# compressing
> gzz
for F in *.postAlleles;do echo "gzip -f $F " >>gzz;done
ls5_launcher_creator.py -j gzz -n gzz -a tagmap -t 0:10:00 -w 24 -e matz@utexas.edu
sbatch gzz.slurm

# on laptop, run write_traits.R to make tables of tratis:
# pd.traits  : proportion of D
# bleach.traits : 2 bleaching traits, blscore, logsymb, and logchl
# blscore.traits : just the bleaching score
# mt.traits : mitochondria type
# rf.traits : reef type, "I" (0) vs "M" (1)
# zz8.ibsMat
# bams.qc
 
# ------ calculating LD and distance where R2 drops below 0.1 (skip this if you want to simply use fixed LD clumping distance)

# subsetting geno files for maf>=0.1
MAF=0.1
>tomaf
for GF in *.split.geno.gz; do
echo "zcat $GF | awk -v maf=$MAF '{ sum=0; for(i=4; i<=NF-1; i=i+3) { i2=i+1; sum+=(\$i+2*\$i2) } ; af=sum/(2*(NF-2)/3); if (af>=maf) { print } }' | gzip > ${GF/.split.geno.gz/}.maf01.geno.gz" >>tomaf
done
ls5_launcher_creator.py -j tomaf -n tomaf -a tagmap -t 0:10:00 -w 24 -e matz@utexas.edu
sbatch tomaf.slurm

# calculating LD within 25kb distance
module load gsl
NB=`cat bams.qc | wc -l`
>ld
KBDIST=25
for GF in *.maf01.geno.gz; do echo "zcat $GF | cut -f 1,2 > ${GF}.sites && NS=\`wc -l ${GF}.sites\` && ngsLD --geno $GF --probs 1 --n_ind $NB --pos ${GF}.sites --n_sites \$NS --max_kb_dist $KBDIST --out ${GF}.LD --n_threads 4 ">> ld; done
ls5_launcher_creator.py -j ld -n ld -a tagmap -t 2:00:00 -w 6 -e matz@utexas.edu
sbatch ld.slurm

# obtaining LD loess (distance at which R2 drops below 0.05)
>ldl
for GF in chr*.maf01.geno.gz.LD; do 
echo "Rscript LDq.R $GF">> ldl; 
done
ls5_launcher_creator.py -j ldl -n ldl -a tagmap -t 2:00:00 -w 4 -e matz@utexas.edu
sbatch ldl.slurm


#-----------
cd 
rm -rf Multivariate_GWAS/
git clone https://github.com/z0on/Multivariate_GWAS.git
cd -

# get glmnet packade for R below 3.6 here: https://packages.debian.org/source/stable/r-cran-glmnet 

# running 50 hold-out replicates on all chromosomes for proportion of D
>pdd
for R in `seq 1 25`; do
REP=rep${R}_25;
for CHR in `ls chr*postAlleles.gz`; do
OUTN=`echo $CHR | perl -pe 's/.+(chr\d+).+/$1/'`;
echo "Rscript ~/Multivariate_GWAS/RDA_GWAS.R gt=$CHR covars.e=reefsites covars.g=technical.covars traits=pd.traits gdist.samples=bams.qc plots=FALSE gdist=zz8.ibsMat hold.out=$REP">>pdd;
done;
done
ls5_launcher_creator.py -j pdd -n pdd -a tagmap -e matz@utexas.edu -t 0:10:00 -w 6 -q normal
sbatch pdd.slurm


>compchrom
ls chr*postAlleles.gz >gts
for R in `seq 1 25`; do
REP=rep${R}_25;
ls chr*${REP}_gwas.RData >gws_${REP};
echo "Rscript ~/Multivariate_GWAS/compile_chromosomes.R in=gws_${REP} gts=gts gt.samples=bams.qc traits=traits_etc_${REP}.RData">>compchrom;
done
ls5_launcher_creator.py -j compchrom -n compchrom -a tagmap -t 2:00:00 -w 1 -e matz@utexas.edu -q normal
sbatch compchrom.slurm

>compchrom0
ls chr*postAlleles.gz >gts
for R in `seq 1 25`; do
REP=rep${R}_25;
ls chr*${REP}_gwas.RData >gws_${REP};
echo "Rscript ~/Multivariate_GWAS/compile_chromosomes.R in=gws_${REP} traits=traits_etc_${REP}.RData">>compchrom0;
done
ls5_launcher_creator.py -j compchrom0 -n compchrom0 -a tagmap -t 2:00:00 -w 24 -e matz@utexas.edu -q normal
sbatch compchrom0.slurm



# running 50 hold-out replicates on all chromosomes for bleaching traits
>bll
for R in `seq 1 50`; do
REP=rep_pdmds2_${R}_10;
MDS=pdmds2_${R}_10;
for CHR in `ls chr*postAlleles.gz`; do
OUTN=`echo $CHR | perl -pe 's/.+(chr\d+).+/$1/'`;
echo "Rscript RDA_GWAS.R gt=$CHR covars=$MDS traits=bleach.traits gdist.samples=bams.qc plots=FALSE gdist=zz8.ibsMat hold.out=$REP outfile=${OUTN}_bl_${REP}.RData">>bll;
done;
done
ls5_launcher_creator.py -j bll -n bll -a tagmap -e matz@utexas.edu -t 0:15:00 -w 12 -q normal
sbatch bll.slurm

# running 50 hold-out replicates on all chromosomes for mt type
>mtt
for R in `seq 1 50`; do
REP=rep${R}_10;
MDS=mds2_${R}_10;
for CHR in `ls chr*postAlleles.gz`; do
OUTN=`echo $CHR | perl -pe 's/.+(chr\d+).+/$1/'`;
echo "Rscript RDA_GWAS.R gt=$CHR covars=$MDS traits=mt.traits gdist.samples=bams.qc plots=FALSE gdist=zz8.ibsMat hold.out=$REP outfile=${OUTN}_mt_${REP}.RData">>mtt;
done;
done
ls5_launcher_creator.py -j mtt -n mtt -a tagmap -e matz@utexas.edu -t 0:10:00 -w 6 -q normal
sbatch mtt.slurm

# running 50 hold-out replicates on all chromosomes for cross-shelf association
>rff
for R in `seq 1 50`; do
REP=rep${R}_10;
MDS=mds2_${R}_10;
for CHR in `ls chr*postAlleles.gz`; do
OUTN=`echo $CHR | perl -pe 's/.+(chr\d+).+/$1/'`;
echo "Rscript RDA_GWAS.R gt=$CHR covars=$MDS traits=rf.traits gdist.samples=bams.qc plots=FALSE gdist=zz8.ibsMat hold.out=$REP outfile=${OUTN}_rf_${REP}.RData">>rff;
done;
done
ls5_launcher_creator.py -j rff -n rff -a tagmap -e matz@utexas.edu -t 0:10:00 -w 6 -q normal
sbatch rff.slurm


# compiling all chromosomes for each replicate
>comprf
for r in `seq 1 50`; do
ls *rf_rep${r}_*.RData >rep${r}_rf; 
echo "Rscript compile_chromosomes.R in=rep${r}_rf plotManhattan=FALSE">>comprf;
done
ls5_launcher_creator.py -j comprf -n comprf -a tagmap -e matz@utexas.edu -t 0:10:00 -w 24 -q normal
sbatch comprf.slurm

>compd
for r in `seq 1 50`; do
ls *pd_rep${r}_*.RData >rep${r}_pd; 
echo "Rscript compile_chromosomes.R in=rep${r}_pd plotManhattan=FALSE">>compd;
done
ls5_launcher_creator.py -j compd -n compd -a tagmap -e matz@utexas.edu -t 0:10:00 -w 24 -q normal
sbatch compd.slurm

>compmt
for r in `seq 1 50`; do
ls *pd_rep${r}_*.RData >rep${r}_pd; 
echo "Rscript compile_chromosomes.R in=rep${r}_pd plotManhattan=FALSE">>compmt;
done
ls5_launcher_creator.py -j compmt -n compmt -a tagmap -e matz@utexas.edu -t 0:10:00 -w 24 -q normal
sbatch compmt.slurm

>compbl
for r in `seq 1 50`; do
ls *bl_rep_pdmds2_${r}_*.RData >rep${r}_bl; 
echo "Rscript compile_chromosomes.R in=rep${r}_bl plotManhattan=FALSE">>compbl;
done
ls5_launcher_creator.py -j compbl -n compbl -a tagmap -e matz@utexas.edu -t 0:10:00 -w 24
sbatch compbl.slurm


#------ FULL RUN:

>bigpd
for CHR in `ls *postAlleles.gz`; do
OUTN=`echo $CHR | perl -pe 's/^([\w\d]+)\\..+/$1/'`;
echo "Rscript RDA_GWAS.R gt=$CHR covars=mds2 traits=pd.traits gdist.samples=bams.qc gdist=zz8.ibsMat outfile=${OUTN}_pd.RData">>bigpd;
done
ls5_launcher_creator.py -j bigpd -n bigpd -a tagmap -e matz@utexas.edu -t 0:15:00 -w 6 -q normal
sbatch bigpd.slurm

>bigrf
for CHR in `ls *postAlleles.gz`; do
OUTN=`echo $CHR | perl -pe 's/^([\w\d]+)\\..+/$1/'`;
echo "Rscript RDA_GWAS.R gt=$CHR covars=mds2 traits=rf.traits gdist.samples=bams.qc gdist=zz8.ibsMat outfile=${OUTN}_rf.RData">>bigrf;
done
ls5_launcher_creator.py -j bigrf -n bigrf -a tagmap -e matz@utexas.edu -t 0:15:00 -w 6 -q normal
sbatch bigrf.slurm

>bigmt
for CHR in `ls *postAlleles.gz`; do
OUTN=`echo $CHR | perl -pe 's/^([\w\d]+)\\..+/$1/'`;
echo "Rscript RDA_GWAS.R gt=$CHR covars=mds2 traits=mt.traits gdist.samples=bams.qc gdist=zz8.ibsMat outfile=${OUTN}_mt.RData">>bigmt;
done
ls5_launcher_creator.py -j bigmt -n bigmt -a tagmap -e matz@utexas.edu -t 0:15:00 -w 6 -q normal
sbatch bigmt.slurm

>bigbl
for CHR in `ls *postAlleles.gz`; do
OUTN=`echo $CHR | perl -pe 's/^([\w\d]+)\\..+/$1/'`;
echo "Rscript RDA_GWAS.R gt=$CHR covars=pdmds2 traits=blscore.traits gdist.samples=bams.qc gdist=zz8.ibsMat outfile=${OUTN}_bl.RData">>bigbl;
done
ls5_launcher_creator.py -j bigbl -n bigbl -a tagmap -e matz@utexas.edu -t 0:15:00 -w 6 -q normal
sbatch bigbl.slurm


# compiling full run data
ls *_mt.RData >mtfull
echo "Rscript compile_chromosomes.R in=mtfull">fullc
ls *_pd.RData >pdfull
echo "Rscript compile_chromosomes.R in=pdfull ">>fullc
ls *_rf.RData >rffull
echo "Rscript compile_chromosomes.R in=rffull ">>fullc
ls *_bl.RData >bleachfull
echo "Rscript compile_chromosomes.R in=bleachfull ">>fullc
ls5_launcher_creator.py -j fullc -n fullc -a tagmap -e matz@utexas.edu -t 0:30:00 -w 4 
sbatch fullc.slurm





#---- older work


# compiling results
for r in `seq 1 50`;do ls *pd_rep${r}_*.RData >rep${r}_pd; done
ls rep*_pd >reps_pd
for r in `seq 1 50`;do ls *mt_rep${r}_*.RData >rep${r}_mt; done
ls rep*_mt >reps_mt

for r in `seq 1 50`;do ls *bl_rep_pdmds2_${r}_*.RData >rep${r}_bl2; done
ls rep[0-9]*_bl2 >reps_bl2
idev -t 1:30:00
Rscript compile_chromosomes.R in=reps_bl2 traits=bleach.traits outfile=bleach_pdmds_50 

for r in `seq 1 50`;do ls *rf_rep${r}_*.RData >rep${r}_rf; done
ls rep[0-9]*_rf >reps_rf
idev -t 1:30:00
Rscript compile_chromosomes.R in=reps_rf traits=rf.traits outfile=rf50 


Rscript compile_chromosomes.R in=reps_pd traits=pd.traits outfile=pd50 
Rscript compile_chromosomes.R in=reps_bl traits=bleach.traits outfile=bl50 
Rscript compile_chromosomes.R in=reps_mt traits=mt.traits outfile=mt50 
ls5_launcher_creator.py -j compiles -n compiles -a tagmap -e matz@utexas.edu -t 2:00:00 -w 4 
sbatch compiles.slurm

scp *_full.* $ARCHIVER:/stornext/ranch_01/ranch/users/01211/cmonstr/rda_gwas_july2020







# ------on laptop:

cd ~/Dropbox/amil_RDA_association_jun2020/RDA_GWAS/
>bll
for R in `seq 1 10`; do
REP=rep${R}_10;
MDS=mds2_${R}_10;
for CHR in `ls ../chr*postAlleles.gz`; do
OUTN=`echo $CHR | perl -pe 's/.+(chr\d+).+/$1/'`;
echo "Rscript RDA_GWAS.R gt=$CHR covars=$MDS traits=bleach.traits gdist.samples=bams.qc plots=FALSE gdist=zz8.ibsMat hold.out=$REP outfile=${OUTN}_bl_${REP}.RData">>bll;
done;
done

>mtt
for R in `seq 1 20`; do
REP=rep${R}_10;
MDS=mds2_${R}_10;
for CHR in `ls ../chr*postAlleles.gz`; do
OUTN=`echo $CHR | perl -pe 's/.+(chr\d+).+/$1/'`;
echo "Rscript RDA_GWAS.R gt=$CHR covars=$MDS traits=mt.traits gdist.samples=bams.qc plots=FALSE gdist=zz8.ibsMat hold.out=$REP outfile=${OUTN}_mt_${REP}.RData">>mtt;
done;
done

>rff
for R in `seq 1 20`; do
REP=rep${R}_10;
MDS=mds2_${R}_10;
for CHR in `ls ../chr*postAlleles.gz`; do
OUTN=`echo $CHR | perl -pe 's/.+(chr\d+).+/$1/'`;
echo "Rscript RDA_GWAS.R gt=$CHR covars=$MDS traits=rf.traits gdist.samples=bams.qc plots=FALSE gdist=zz8.ibsMat hold.out=$REP outfile=${OUTN}_rf_${REP}.RData">>rff;
done;
done

for r in `seq 1 7`;do ls *rf_rep${r}_*.RData >rep${r}_rf; done
ls rep*_rf >reps_rf
Rscript compile_chromosomes.R in=reps_rf traits=rf.traits outfile=rf7


# full runs




#----- older work
>sc
for PA in *.postAlleles; do
echo "Rscript rda_scores.R $PA">>sc;
done
ls5_launcher_creator.py -j sc -n sc -t 1:00:00 -e matz@utexas.edu -w 4 -a tagmap
sbatch sc.slurm

>rda
for PA in *.postAlleles; do
echo "Rscript rda_manhattans_auto.R $PA">>rda;
done
ls5_launcher_creator.py -j rda -n rda -t 1:00:00 -e matz@utexas.edu -w 12 -a tagmap
sbatch rda.slurm

>rda
for PA in chr*.postAlleles; do
echo "Rscript twoRDAs_mds.R $PA">>rda;
done
ls5_launcher_creator.py -j rda -n rda -t 1:00:00 -e matz@utexas.edu -w 12 -a tagmap
sbatch rda.slurm


>rda100
for PA in chr*.postAlleles; do
echo "Rscript 100gwasTests.R $PA">>rda100;
done
ls5_launcher_creator.py -j rda100 -n rda100 -t 4:00:00 -e matz@utexas.edu -w 6 -a tagmap -q normal
sbatch rda100.slurm

>coll
for REP in `seq 1 50`;do
echo "Rscript polygenic_scores_allchroms.R $REP">>coll;
done
ls5_launcher_creator.py -j coll -n coll -t 2:00:00 -e matz@utexas.edu -w 6 -a tagmap 
sbatch coll.slurm

>rda100m
for PA in chr*.postAlleles; do
echo "Rscript 100gwasTests_mtrf_1_50.R $PA">>rda100m;
done
ls5_launcher_creator.py -j rda100m -n rda100m -t 2:00:00 -e matz@utexas.edu -w 6 -a tagmap
sbatch rda100m.slurm


>rda100b
for PA in chr*.postAlleles; do
echo "Rscript 100gwasTests_mtrf_50_100.R $PA">>rda100b;
done
ls5_launcher_creator.py -j rda100b -n rda100b -t 2:00:00 -e matz@utexas.edu -w 6 -a tagmap
sbatch rda100b.slurm

>coll
for REP in `seq 50 86`;do
echo "Rscript polygenic_scores_allchroms_rf.R $REP
Rscript polygenic_scores_allchroms_mt.R $REP">>coll;
done
ls5_launcher_creator.py -j coll -n coll -t 2:00:00 -e matz@utexas.edu -w 12 -a tagmap 
sbatch coll.slurm

>rdall
for PA in *.postAlleles; do
echo "Rscript 1gwasTest.R $PA">>rdall;
done
ls5_launcher_creator.py -j rdall -n rdall -t 1:00:00 -e matz@utexas.edu -w 6 -a tagmap
sbatch rdall.slurm

# max z-scores per gene (correncted by N SNPs per gene)
grep gene /Users/c-monstr/Dropbox/Amil_v2.01_annotated/Amil.coding.gff3 | awk '$3=="gene"' | cut -f 1,4,5,9 | perl -pe 's/ID=(Amillepora\d+).+/$1/' >amil_genes.gff3
# scp amil_genes.gff3 and gwas_rda_zscore_byGene.R to tacc
ls *d_full.RData >infiles
echo "Rscript gwas_rda_max_zscore_byGene.R
Rscript gwas_rda_max_logP_byGene.R" >byg
ls5_launcher_creator.py -j byg -n byg -t 0:15:00 -e matz@utexas.edu -w 2 -a tagmap
sbatch byg.slurm


awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$6 }' Amillepora_euk.emapper.annotations | grep GO | perl -pe 's/,/;/g' >amil_gene2go.tab
#  KOG classes (single-letter):
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$12 }' Amillepora_euk.emapper.annotations | grep -Ev "[,#S]" >amil_gene2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt amil_gene2kogClass1.tab > amil_gene2kogClass.tab



# for -minMaf 0.05
export FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 152 -snp_pval 1e-5 -minMaf 0.05 "
export TODO8="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1"
echo '#!/bin/bash
#SBATCH -J zz8
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p largemem512GB
#SBATCH -o zz8.o%j
#SBATCH -e zz8.e%j
#SBATCH -t 24:00:00
#SBATCH -A tagmap
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
angsd -b bams.traits -GL 1 $FILTERS $TODO8 -P 12 -out zz805' > zz805
zz805job=$(sbatch zz805 | grep "Submitted batch job" | perl -pe 's/\D//g')

# for -minMaf 0.01
export FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 152 -snp_pval 1e-5 -minMaf 0.01 "
export TODO8="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1"
echo '#!/bin/bash
#SBATCH -J zz8
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p largemem512GB
#SBATCH -o zz8.o%j
#SBATCH -e zz8.e%j
#SBATCH -t 24:00:00
#SBATCH -A tagmap
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
angsd -b bams.traits -GL 1 $FILTERS $TODO8 -P 12 -out zz801' > zz801
zz801job=$(sbatch zz801 | grep "Submitted batch job" | perl -pe 's/\D//g')

#zcat zz8.mafs.gz | wc -l
# 7218426

#---------------- skip
# rerunning for -doGeno 11 
export TODO11="-doMajorMinor 1 -doMaf 1 -doPost 2 -doGeno 11"
echo '#!/bin/bash
#SBATCH -J zz11
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p largemem512GB
#SBATCH -o zz11.o%j
#SBATCH -e zz11.e%j
#SBATCH -t 24:00:00
#SBATCH -A tagmap
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
angsd -b bams.qc -GL 1 $FILTERS $TODO11 -P 12 -out zz11 -checkBamHeaders 0' > zz11
zz11job=$(sbatch zz11 | grep "Submitted batch job" | perl -pe 's/\D//g')

# selecting good sites (not lumped paralogs) 
module load python2
echo "zcat zz11.geno.gz | python ~/bin/HetMajorityProb.py | awk '\$6 < 0.75 {print \$1\"\t\"\$2}' > allSites" >hp
ls5_launcher_creator.py -j hp -n hp -t 6:00:00 -e matz@utexas.edu -w 1 -a tagmap -q normal
hpjob=$(sbatch hp.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
#-------------------

zcat zz8.mafs.gz | cut -f 1,2 | tail -n +2 >allSites

>thi25
for R in `seq 1 96`; do
echo "thinner.pl infile=allSites interval=25000 criterion=random >z25k_sites_$R">>thi25;
done
ls5_launcher_creator.py -j thi25 -n thi25 -t 2:00:00 -e matz@utexas.edu -w 96 -a tagmap -q normal
thijob=$(sbatch thi25.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

# unzipping geno file
echo "zcat zz8.geno.gz > zz8.geno" >zz
ls5_launcher_creator.py -j zz -n zz -t 2:00:00 -e matz@utexas.edu -w 1 -a tagmap -q normal
zzjob=$(sbatch zz.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

module load gsl
export NB=`cat bams.qc | wc -l`
>ld
for R in `seq 1 96`; do
echo "awk 'NR==FNR{a[\$1\$2]; next} \$1\$2 in a{print}' z25k_sites_$R zz8.geno > zz_$R.geno && NS=\`wc -l zz_$R.geno\` && cut -f 1,2 zz_$R.geno > zz_$R.geno.sites && gzip zz_$R.geno -f && ngsLD --geno zz_$R.geno.gz --probs 1 --n_ind $NB --n_sites \$NS --max_kb_dist 0 --pos zz_$R.geno.sites --out z25k_$R.LD --n_threads 12 --extend_out 1 ">> ld;
done
ls5_launcher_creator.py -j ld -n ld -t 12:00:00 -e matz@utexas.edu -w 8 -a tagmap -q normal
LDjob=$(sbatch --dependency=afterok:$zzjob ld.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

# geno files for wgcna
#for R in `seq 1 96`; do
#awk 'NR==FNR{a[$1$2]; next} $1$2 in a{print}' z25k_sites_$R zz8.geno > zz8_$R.geno;
#done


#------ collecting sites that correlate better than 0.1:

>si
for R in `seq 1 96`; do
echo "awk '\$7>0.1' z25k_$R.LD | awk '{a=\$1 \"\\n\" \$2; print a}' >z25k_$R.LD.sites">>si;
done
ls5_launcher_creator.py -j si -n si -t 1:00:00 -e matz@utexas.edu -w 24 -a tagmap 
sijob=$(sbatch si.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
#sijob=$(sbatch --dependency=afterok:$LDjob si.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

cat *LD.sites | sort | uniq | perl -pe 's/:/\t/' | sort -k 1,1 -k 2,2n > LDsites2
wc -l LDsites2
#452342
#437885
#822482 rerun with 96

>thiL
for R in `seq 1 96`; do
echo "thinner.pl infile=LDsites2 interval=20000 criterion=random >LD_sites20_$R">>thiL;
done
ls5_launcher_creator.py -j thiL -n thiL -t 0:30:00 -e matz@utexas.edu -w 24 -a tagmap 
#thiLjob=$(sbatch --dependency=afterok:$sijob thiL.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
thiLjob=$(sbatch thiL.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
wc -l LD_sites20_1

# 12358 LD_sites20_1
# 12266
# 12279 rerun 96

# reanalyzing LD
module load gsl
export NB=`cat bams.qc | wc -l`
>LD
for R in `seq 1 96`; do
echo "awk 'NR==FNR{a[\$1\$2]; next} \$1\$2 in a{print}' LD_sites20_$R zz8.geno > LD_$R.geno && NS=\`wc -l LD_$R.geno\` && cut -f 1,2 LD_$R.geno > LD_$R.geno.sites && gzip LD_$R.geno -f && ngsLD --geno LD_$R.geno.gz --probs 1 --n_ind $NB --n_sites \$NS --max_kb_dist 0 --pos LD_$R.geno.sites --out LD_$R.LD --n_threads 12 --extend_out 1 ">> LD;
done
ls5_launcher_creator.py -j LD -n LD -t 2:00:00 -e matz@utexas.edu -w 4 -a tagmap -q normal
ld2job=$(sbatch --dependency=afterok:$thiLjob LD.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
#ld2job=$(sbatch LD.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

>sl
for R in `seq 1 96`; do
echo "cut -f 1,2,7 LD_$R.LD > LD_$R.LD.slim">>sl;
done
ls5_launcher_creator.py -j sl -n sl -t 0:30:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sljob=$(sbatch --dependency=afterok:$ld2job sl.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
#sljob=$(sbatch sl.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

>lm5
for R in `seq 1 96`; do
echo "Rscript ~/bin/ld2matrix.R LD_$R.LD.slim">>lm5;
done
ls5_launcher_creator.py -j lm5 -n lm5 -t 48:00:00 -e matz@utexas.edu -w 4 -a tagmap -q normal
ld2mjob=$(sbatch --dependency=afterok:$sljob lm5.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')


# serial WGCNA
>ldss
# for -doGeno 8 (posterior number of derived alleles):
for I in `seq 1 96`;do echo "Rscript ~/bin/LD_WGCNA_auto.R LD_$I.LD.slim_matrix.RData bams.qc LD_$I.geno.gz">>ldss;done
ls5_launcher_creator.py -j ldss -n ldss -t 0:30:00 -e matz@utexas.edu -w 6 -a tagmap -q normal 
#sbatch ldss.slurm
ldssjob=$(sbatch --dependency=afterok:$ld2mjob ldss.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

echo "Rscript ~/bin/collectRawData.R 96" >coll
ls5_launcher_creator.py -j coll -n coll -t 2:00:00 -e matz@utexas.edu -w 1 -a tagmap -q normal 
# sbatch coll.slurm
sbatch --dependency=afterok:$ldssmjob coll.slurm
