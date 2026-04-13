## Homework 10: Mendelian Randomization
### BS859 Applied Genetic Analysis
### Addison Yam
### April 15, 2026

```bash

```

In class we investigated causal relationships between Body Mass Index (BMI) and Type 2 Diabetes (T2D), a dichotomous trait. We observed a causal relationship where for every SD increase in BMI, the odds of T2D increased by 3 times (OR 3.00). There are many other factors that have been observed to affect T2D risk. In this assignment, we will investigate causal relationships between High Density Lipoprotein cholesterol (HDL-c) and Fasting Glucose (FG) with T2D. 

Outcome Dataset (For Q1&2): Type 2 Diabetes GWAS Meta Analysis by Xue et al. 60,000 Cases and 600,000 controls of European ancestry that we used in class. (Xue, A., Wu, Y., Zhu, Z. et al. Genome-wide association analyses identify 143 risk variants and putative regulatory mechanisms for type 2 diabetes. Nat Commun 9, 2941 (2018). https://doi.org/10.1038/s41467-
018-04951-w)
Location:/projectnb/bs859/data/T2D/T2D_Xue_et_al_2018.txt

1. Fasting glucose (FG) is a blood biomarker used in the diagnosis of T2D, with higher FG indicating T2D. This trait has also been studied through genetic association studies. One example is a meta-analysis by Chen et al. ("The trans-ancestral genomic architecture of glycemic traits." Nature genetics 53.6 (2021): 840-860. http://nature.com/articles/s41588-
021-00852-9 ) Chen et al. (2021) analyzed 281,416 individuals of European ancestry. The file containing the
summary statistics for European ancestry individuals can be found in FG folder on scc: /projectnb/bs859/data/FG/MAGIC1000G_FG_EUR.tsv
There is a readme file that will tell you about the column names.
This file was downloaded from https://magicinvestigators.org/downloads/

a. Reformat the FG-GWAS file into GCTA format as we did in class for the T2D file and create the appropriate “exposure” and “outcome” files that list the file names for GSMR.
b. Use GSMR to test the causal hypothesis that FG results in a change in odds of T2D. You may use the default parameters for GSMR.

List the effect estimate, standard error P-value for the causal effect along with the number of SNP instruments

1c. Provide an Odds Ratio for the causal effect of FG and T2D, and interpret it with appropriate units. Is this estimate expected given what we know about FG & T2D?
1d. Provide a plot of the SNP instrument effects for FG & T2D. Comment on the features of the plot and if the plot agrees with the odds ratio from 1c.
1e. GSMR used the HEIDI outlier method to remove outlier SNPs. How many SNPs were removed from the FG exposure? Explain why these variants were removed.
1f. Report (and interpret) the result of causal analysis in the other direction (T2DFG), using the default GSMR parameters. 


Question 2

2a. Draw a causal diagram similar to the one seen on page 14 of the class 11 notes for the Mendelian Randomization (MR) framework addressed in Question 1. Namely, the causal relationship between FG & T2D (where FG is the exposure and T2D is the outcome).
2b. Explain briefly in your own words the following MR assumptions in the context of the causal diagram you drew in 2a.
Relevance Assumption:
Independence Assumption:
Exclusion Assumption:
2c. Given what we know about how T2D is diagnosed, which direction of the causal effects do you think is most likely correct?
2d. Thinking about the assumptions you explained in 2b., what possible explanations for the bidirectional result do you think are plausible?

Question 3
High Density Lipoprotein cholesterol (HDL-c) is a blood-lipid biomarker related to many health outcomes, but not directly involved in the diagnosis T2D. In general, higher levels of HDL-c is considered protective with regards to negative health outcomes such as heart disease. This trait has been extensively studied by GWAS, and meta-analysis summary statistics by Graham et al. (Graham, S.E., Clarke, S.L., Wu, KH.H. et al. The power of genetic diversity in genome-wide association studies of lipids. Nature 600, 675–679 (2021). https://doi.org/10.1038/s41586-021-04064-3) are provided for this assignment. The authors investigated genetic associations with HDL-c
among 1,320,016 individuals of European ancestry. The summary statistics can be found in /projectnb/bs859/data/lipids folder on scc: HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results

3a. Reformat the HDL-GWAS file into GCTA format as we did in class for the T2D file.
Make sure to include explicit NA values in the output.
Note: Alt is the effect allele
Create the appropriate “exposure” and “outcome” files for GSMR.
3b. Use GSMR to test the causal hypothesis that HDL-c results in a change in odds of T2D. You may use the default parameters for GSMR.
Please list the effect estimate, standard error P-value for the causal effect along with the number of SNP instruments.
3c. How many pleiotropic HDL SNPs were removed by GSMR?
3d. Provide an Odds Ratio for the causal effect of HDL-c and T2D, and interpret it with appropriate units. Is this estimate expected given what we know about HDL-c in general? 

3e. Provide a plot of the SNP instrument effects for HDL-c & T2D. Comment on the features of the plot and if the plot agrees with the odds ratio from 2c. Do you see any potential outliers remaining? (Remember that by default the HEIDI procedure removes outliers at threshold p<0.01.)
3f. Rerun GSMR to test for a causal relationship between HDL-c and T2D as in 2b, but use a more strict significance threshold to select SNPs of 1×10-10. Report the effect estimate and PValue and explain any similarities/differences from 2b.
3g. Rerun GSMR to test for a causal relationship between HDL-c and T2D as in 2b, but use a more strict significance threshold for the HEIDI outlier removal – use p=0.05. Plot the SNP effects, and report the effect estimate and P-Value and explain any similarities/differences from 2b. Have some of the more extreme instruments in the earlier plot been removed?