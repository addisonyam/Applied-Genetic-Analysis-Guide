##load plink and gcta modules
module load gcta/1.94.1
module load plink/1.90b6.21

DATA2='/projectnb/bs859/data/meta/downloads'
cat $DATA2/UKB.v2.albuminuria.n382500.README.txt
zcat $DATA2/UKB.v2.albuminuria.n382500.sorted.tsv.gz|head


zcat $DATA2/UKB.v2.albuminuria.n382500.sorted.tsv.gz|wc


zcat $DATA2/UKB.v2.albuminuria.n382500.sorted.tsv.gz|awk 'NR==1||$10<5e-8 {print $0}'  >tophits.txt

wc tophits.txt
cat tophits.txt

##the top locus is on chromosome 10, has multiple variants with p<5e-50
zcat $DATA2/UKB.v2.albuminuria.n382500.sorted.tsv.gz|awk 'NR==1||$10<5e-50 {print $0}'  > p5e-50.txt


##remind me of the format of my GWAS summary statistics:
zcat $DATA2/UKB.v2.albuminuria.n382500.sorted.tsv.gz|head


###create reformatted file for GCTA conditional analysis 

zcat $DATA2/UKB.v2.albuminuria.n382500.sorted.tsv.gz|awk '{print $2,$5,$6,$7,$8,$9,$10,$11}' > togcta.txt

head togcta.txt

##conditional analysis, using 1000G europeans as the reference group from 
##which LD is estimated

##directory has 1000G plink files for all 5 1000G super populations and combined
ls /projectnb/bs859/data/1000G/plinkformat/

## 503 individuals in the EUR 1000G  
##77.8 million variants in the files
wc /projectnb/bs859/data/1000G/plinkformat/1000G_EUR.fam
wc /projectnb/bs859/data/1000G/plinkformat/1000G_EUR.bim

##is our top variant in the 1000G EUR data?
grep rs141640975 /projectnb/bs859/data/1000G/plinkformat/1000G_EUR.bim
##what is its frequency?
plink --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR --chr 10 --freq --out chr10freq --allow-extra-chr
grep rs141640975 chr10freq.frq

##run the conditional analysis for chromosome 10 only.  This will take a while!
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR --cojo-file togcta.txt --cojo-slct --chr 10 --out chr10 > chr10.log 

wc chr10.badsnps
wc chr10.freq.badsnps

wc chr10.jma.cojo

##re-run conditional analysis, excluding variants with MAF<0.01:
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR  --cojo-file togcta.txt  --cojo-slct --chr 10 --maf 0.01 --out chr10.maf01 > chr10.maf01.log


##re-run conditional analysis, with the wrong LD reference panel for these data:
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_AFR  --cojo-file togcta.txt  --cojo-slct --chr 10  --out chr10.AFRld > chr10.AFRld.log


##  extract the rsids from the variants with p<5e-50 
## we are going to use VEP (Variant Effect Predictor) to get
## annotations for these 18 variants.
cut -f2 p5e-50.txt > toVEP.txt 

