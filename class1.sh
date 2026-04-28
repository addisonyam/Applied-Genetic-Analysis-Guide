#set a shell variable  to path to class data.  Now we can use $DATADIR to refer to the longer directory name: 

export DATADIR=/projectnb/bs859/data/plink_tutorial  

#load plink software module:
module load plink/1.90b6.27 

#what's in the plink_tutorial directory we refer to as $DATADIR ?
ls -lat $DATADIR

#Look at the first 10 lines of the file wgas1_bin.fam 
head $DATADIR/wgas1_bin.fam     

#Look at the first 10 lines of the file wgas1_bin.bim 
head $DATADIR/wgas1_bin.bim     

# How many lines are there in the fam file?   "wc" == wordcount 
# The first number returned is the number of lines in the file
# *.fam is the pedigree file, which does not have a header, 
# so the number of lines==the number of individuals in the file
wc $DATADIR/wgas1_bin.fam  

# How many lines are there in the bim file? 
# *.bim is the map file, so the number of lines==the number of genetic variants in the data set

wc $DATADIR/wgas1_bin.bim 

##run some plink analyses. 

##compute the allele frequencies, and save them to a file freq1.frq (PLINK automatically puts the .frq extension on the filename.  The --out freq1 is providing the base of the name)
plink --bfile $DATADIR/wgas1_bin --freq --out freq1   

##sometimes we want genotype counts for all 3 genotypes instead of 
## just allele frequencies.  We can use --freqx to get that
plink --bfile $DATADIR/wgas1_bin --freqx --out freqx      
## The genotype counts are also found in the HWE output, and in the output of "--freq counts"

##compute the missing rates for individuals and variants, and save them to files: miss1.imiss miss1.lmiss 

plink --bfile $DATADIR/wgas1_bin --missing --out miss1    

##compute a hardy-weinberg test for each variant, and save the results to a file  hwe1.hwe

plink --bfile $DATADIR/wgas1_bin --hardy --out hwe1   

##how big are the files that are produced by the plink analyses we've performed?
wc freq1.frq
wc miss1.imiss   #individual missing rates
wc miss1.lmiss   #variant missing rates
wc hwe1.hwe

#what do the results in these files look like?  (look at firset 10 lines)
head freq1.frq
head miss1.imiss
head miss1.lmiss
head hwe1.hwe


awk '($3=="UNAFF"&&$9<0.0001)||NR==1{print $0}'  hwe1.hwe>hwe.unaff.0001.txt
awk 'NR==1||$6>0.03{print $0}' miss1.imiss>highmissingindiv.txt

cat hwe.unaff.0001.txt
cat highmissingindiv.txt

##filtering the file and creating a new, filtered data set "wgas2" using plink's standard order of operations:
plink --bfile $DATADIR/wgas1_bin --maf 0.01 --geno 0.05 --mind 0.05 --hwe 1e-4 --make-bed --out wgas2 
##

##filtering the file and creating a new, filtered data set "wgas2b" forcing plink to first omit variants with MAF<0.01
## and genotype missing rate>0.05, and THEN omit based on individual missing rate >0.05 and hwe p-value <0.0001 in controls:
plink --bfile $DATADIR/wgas1_bin --maf 0.01 --geno 0.05 --make-bed --out wgas2a  
plink --bfile wgas2a --mind 0.05 --hwe 1e-4 --make-bed --out wgas2b  


##sex check -- compare sex in the input file to the genotypes on chromosome X
plink --bfile wgas2 --check-sex --out sex 


##ld prune
plink --bfile wgas2 --indep-pairwise 10000 kb  1  0.2 --out try1
##Produces try1.prune.in and try1.prune.out (and .log).
## we want to keep the snps listed in try1.prune.in,
## or remove the snps listed in try1.prune.out
##its also a good idea to exclude X chromosome from this
## analysis --chr 1-22 keeps only the chromosomes 1 through 22

## compute the heterozyogote deficit F statistics:
plink --bfile wgas2 --chr 1-22 --extract try1.prune.in --het --out Fstat

module load R   ##need to load R module before using R

## run my Fstat.R program through R and put the output in Fstat.log (no need to load R or Rstudio to run a quick script):
Rscript Fstat.R>Fstat.log


##IBD estimation with the prune-in SNP subset:
plink --bfile wgas2 --chr 1-22 --extract try1.prune.in --genome --out ibd 
wc ibd.genome
head ibd.genome

awk '$10>0.05{print $0}' ibd.genome > ibd.gt05.genome
##using the --min 0.05 flag will keep only the IBD pairs that have pihat>=0.05
## so the plink command below accomplishes the same thing as the plink 
##command above + the awk command
plink --bfile wgas2 --chr 1-22 --extract try1.prune.in --genome --min 0.05  --out ibd05  


