

##load R, plink, and eigensoft 
##smartpca is a program in the eigensoft module
module load R
module load plink/1.90b6.27
module load eigensoft  

#look at example parameter file
cat test.par

# plink fileset /projectnb/bs859/data/plink_tutorial/wgas2/wgas2_pruned.* 
# is the cleaned and pruned file set we used last week in class for estimating
# individual heterozygosity and pairwise IBD.
# Here, we use the same pruned data to do the PCA.
# The output files will be "test.evec" (the eigenvectors, or as we call them, 
# the PCs), and in test.out, there will be a lot of useful information that 
# would otherwise be put on the screen while smartpca is running.
# smartpca tests the association of the PCs with case status (or whatever
# phenotype is in the input *.fam file in the 6th column)

#run smartpca  using the parameter file test.par
smartpca -p test.par > test.out

more test.out

#Plot PCs by case status  - 
#I wrote an R script to do this that is
# somewhat flexible.  The command takes 4 arguements:
# 1) The name of the output file from smartpca
# 2 and 3) the two PCs you want to plot on x and y axes, respectively
# 4) Number of PCs in the file 
# this script assumes the output is from smartpca, so the 
# first column is the individual ID and the last column is 
# case status (from the plink fam file used to run smartpca)
# produces a simple scatterplot with name arg1.PC.arg2.arg3.jpeg 
Rscript --vanilla plotPCs.R  test.evec 1 2 10 
Rscript --vanilla plotPCs.R  test.evec 1 3 10 
Rscript --vanilla plotPCs.R test.evec 2 3 10 

## the  --vanilla argument when invoking Rscript makes sure that 
## user-specific R settings are ignored, and that there is no saving 
## or restoring of workspaces. This makes the script more portable to 
## other systems


# Print the first column of the file test.evec:
awk '{print $1}' test.evec

#We will add population information to evec file so we can plot the 
#PCs by origin of the sample instead of case status.  
#Note that for this sample, the individual IDs all start with 
#either JA or CH.  JA==Japan, and CH==China.  We will take advantage of
#that!
#Since I want to use my R script to do the plot, I have to replace case
#status with JA and CH but leave the other columns where my script
#expects them to be.

awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11, substr($1,1,2)}' test.evec>test.evec2

#last column of test.evec2 should now be CH or JA to indicate 
#the sample's population of origin, and columns 1-11 
#are identical to test.evec:
head test.evec2
tail test.evec2 
##the only other difference is that now the file has a single space delimiter.

#replot PCs, this time individuals are colored by 
# population instead of case status-- 
Rscript --vanilla plotPCs.R  test.evec2 1 2 10
Rscript --vanilla plotPCs.R test.evec2 1 3 10 
Rscript --vanilla plotPCs.R test.evec2 2 3 10 

########################################### 
##Analyses with TGEN cleaned data set below

## get short names for the directories we will be using in the next analyses:
DATADIR=/projectnb/bs859/data/tgen/cleaned
HAPMAPDIR=/projectnb/bs859/data/tgen/hapmap

##Use awk get a list of the SNPs in the TGEN_cleaned map (.bim) file where the two alleles are A and T,  or C and G:
awk '($5=="A"&&$6=="T")||($5=="T"&&$6=="A")||($5=="G"&&$6=="C")||($5=="C"&&$6=="G"){print $2}' $DATADIR/TGEN_cleaned.bim > atcgSNPs_omit.txt
head atcgSNPs_omit.txt
#How many SNPs are A/T or G/C?  
wc atcgSNPs_omit.txt




##Remove this list of SNPs from the TGEN_cleaned file set, and create a new file set called TGEN_cleaned1 
##(this file set will be in your current working directory)
plink   --bfile $DATADIR/TGEN_cleaned --exclude atcgSNPs_omit.txt --make-bed --out TGEN_cleaned1 --allow-no-sex

#how many SNPs are in the updated fileset?  Can look at the plink output, or use wc to check the bim file:
wc TGEN_cleaned1.bim

##Now we will attempt to merge the hapmap CEU data wtih our TGEN data. 
##This merge will keep all the variants that are just in CEU or just in TGEN
## by adding the --geno 0.03, we are omitting variants with  >3% missing, 
## which will exclude the variants that are only in one of the 2 data sets.  

plink   --bfile TGEN_cleaned1 --bmerge $HAPMAPDIR/ceu_tgen --geno 0.03 --make-bed --out merge1

##Since there were mismatched alleles, we need to "flip" our data for some variants,  before we can do the merge.  The *.missnp file gives the
##list of SNPs plink found that did not have matching alleles:
wc merge1-merge.missnp
head merge1-merge.missnp

# Flip the alleles of the list of snps, and make a new file set called TGEN_cleaned2:
plink --bfile TGEN_cleaned1 --flip merge1-merge.missnp --make-bed --out TGEN_cleaned2
#try to merge with the hapmap CEU again:
plink --bfile TGEN_cleaned2 --bmerge $HAPMAPDIR/ceu_tgen  --geno 0.03 --make-bed  --out merge2
##remove the 7 remaining problematic variants. We don't know exactly what the issue is, but need all of these for the PCA:
plink --bfile TGEN_cleaned2 --exclude merge2-merge.missnp --make-bed --out TGEN_cleaned3 --allow-no-sex

##merge the TGEN data with each of the 3 hapmap files and include only the variants in the TGEN file 
##create the list of variants we want to keep after merging:
cut -f2 TGEN_cleaned3.bim>keepsnps.txt  #list of variants in the TGEN_cleaned3 file
# Merge TGEN cleaned data with Hapmap CEU
plink --bfile TGEN_cleaned3 --allow-no-sex --bmerge $HAPMAPDIR/ceu_tgen --extract keepsnps.txt --make-bed --out merge3
#Merge the TGEN+CEU data with the hapmap CHBJPT samples:
plink --bfile merge3 --allow-no-sex --bmerge $HAPMAPDIR/chbjpt_tgen --extract keepsnps.txt --make-bed --geno 0.03 --out merge4
#THere were a few more  problem SNPs in this merge.  We will exclude them from the TGEN+CEU merged file and then re-merge with CHBJPT:
plink --bfile merge3 --allow-no-sex --exclude merge4-merge.missnp --make-bed --out merge3a
plink --bfile merge3a --allow-no-sex --bmerge $HAPMAPDIR/chbjpt_tgen --extract keepsnps.txt  --make-bed --geno 0.03 --out merge4

##now, merge the TGEN+CEU+CHBJPT with hapmap YRI:
plink --bfile merge4 --allow-no-sex --bmerge $HAPMAPDIR/yri_tgen  --extract keepsnps.txt --make-bed --geno 0.03 --out merge5

#Now, we prune prior to PCA
plink --bfile merge5 --indep-pairwise 10000kb 1 0.2 --out prune0.2 --allow-no-sex
##how many variants do we keep?
wc prune0.2.prune.in
## Keep the pruned-in variants, and keep only autosomal variants (chr 1-22)
## and create a new file set called TGEN_cleaned1a:
plink --bfile merge5 --extract prune0.2.prune.in --chr 1-22 --make-bed --out TGEN_hapmap_pruned --allow-no-sex

wc TGEN_hapmap_pruned.bim

#how many variants, and how many individuals, do we have in our final merged fileset? 
wc TGEN_hapmap_pruned.bim
wc TGEN_hapmap_pruned.fam

##note that the phenotype column in merge5.fam has 1 and 2 for controls and cases from tgen;  in the hapmap files, I gave
## the CEU phenotype= 3, the CHB/JPT phenotype =4, and the YRI phenotype = 5.  THis will make it easier to plot the results of the PCA
head  TGEN_hapmap_pruned.fam

##pca on hapmap and tgen data

## look at the smartpca hapmap parameter file that I provided.  Confirm that it specifies the analyses we want to do:
cat hapmap.par
## run smartpca using the hapmap-tgen merged file set:
smartpca -p hapmap.par > hapmaptgen.out

##look at the output (we put it in the file hapmaptgen.out)
more hapmaptgen.out

## take a quick look at the file that holds the PCs computed by smartpca:
## Note that the last column of the file is the phenotype from the *.fam file we used to run the PCA
## smartpca translates the "1" into "Control" and the "2" into "Case" in the evec file:
head TGEN_hapmap_pruned.evec

## Plot the first first three PCs by case/control/Hapmap population 

Rscript --vanilla plotPCs.R TGEN_hapmap_pruned.evec 1 2 10 
Rscript --vanilla plotPCs.R TGEN_hapmap_pruned.evec 1 3 10 
Rscript --vanilla plotPCs.R TGEN_hapmap_pruned.evec 2 3 10 

# identify the tgen data set outliers (selected PC values based on looking at the previous plot)
# note that first column is familyid:id, and 2nd--11th columns are the 10 PCs, and column 12 is the phenotype/hapmap population  
awk '$2<-0.01&&($12=="Case"||$12=="Control"){print $0}' TGEN_hapmap_pruned.evec> outliers.txt
## Who are the outliers?
awk '{print $1,$2,$3,$4}' outliers.txt

cut -d \t -f 1-4 outliers.txt
