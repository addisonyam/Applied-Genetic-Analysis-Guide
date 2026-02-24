## Homework 5: Association analysis of Rare Variants and Sequence Data
### BS859 Applied Genetic Analysis
### Addison Yam
### February 25, 2026

Use the annotated imputed genotype gds files in /projectnb/bs859/data/tgen/annotated_imputed_vcfs (created from the gene annotated imputed genotype vcfs in the same directory) to complete this homework assignment.
The TGEN Alzheimer disease phenotype data are in tgen.psam, and the genetic relationship data are in grm.rel and grm.rel.id in the same directory. Use the principal component data in /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt

```bash
# load the necessary modules
module load R 
module load plink2/alpha3.7
#make an alias for the directory with the data
export DATADIR=/projectnb/bs859/data/tgen/cleaned/
```


1. Perform gene-based tests on exonic variants on chromosome 19 twice. Once restricting to MAF<=0.05 and once restricting to MAF<=0.01 for chromosome 19 (we did the MAF<=0.05 analysis in class). Use the same other parameters that we used in class. Fill in the table below with the results for b. through g. Don’t forget to also answer h!
a. Show your SMMAT call for the MAF<=0.01 version and explain what each of the parameters in the call is
doing. Consult the SMMAT documentation!
b. How many individuals were included in the analysis?
c. How many of the gene-based tests included only 1 variant?
d. How many of the gene-based tests that include 2 or more variants have a cumulative minor allele count (cMAC)<10?
e. To achieve a family wise error rate of 0.05, what Bonferroni-adjusted significance level should we use to determine significant gene associations, if we use just the SMMAT-E test, and only the genes with cMAC>=10 and at least 2 variants?
f. How many genes are significant by the SMMAT-E test? List them
g. Identify the most significant gene from SMMAT-E with MAF<0.05. For BOTH the MAF<0.05 and MAF<0.01 analysis, report the SMMAT-E p-value for that gene and gene cMAC in the table.
h. What do your observations in g. tell you about the variants that are driving that gene association?

|  | MAF<0.05 SNPs | MAF<0.01 SNPs |
|  ---  |   ----  | --- |
| b. Number of individuals in the analysis |    |   |
| c. Number of gene-based tests with only 1 variant included |    |   |
| d. Gene-based tests with >= 2 variants and cMAC>=10 |    |   |
| e. Bonferroni-adjusted significance level using only SMMAT-E, and genes with cMAC>=10 and >=2 genetic variants |    |   |
| f. Number of significant genes |    |   |
| g. What is the SMMAT-E p-value for the most significant gene from the MAF<0.05 SMMAT-E analysis? |    |   |

2. Perform gene-based testing for the TGEN imputed data on chromosome 19, using a maximum allele frequency of 0.05 and flat weights (MAF.weights.beta=c(1,1)), so that all variants are equally weighted rather than up-weighting rarer variants.
a. Omit the gene-based tests with only a single variant and with cMAC<10. How many genes remain?
b. Create a scatter plot -log10(p) for the SMMAT-E p-value for the flat weights versus the default MAFweights.beta=c(1,25) weights.
c. Summarize and interpret your findings concerning the effect of the weighting scheme on the associations. What do your observations tell you about the most significant gene association?