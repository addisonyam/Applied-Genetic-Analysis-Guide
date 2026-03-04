## Homework 7: Post-GWAS Analyses - Conditional Analyses and Find Mapping
### BS859 Applied Genetic Analysis
### Addison Yam
### March 18, 2026

```bash
# load the necessary modules
module load gcta/1.94.1
module load plink/1.90b6.21

#make an alias for the directory with the data
export DATADIR=/projectnb/bs859/data/tgen/annotated_imputed_vcfs/
```

Use the GWAS results from the Albuminuria GWAS located in `/projectnb/bs859/data/meta/downloads/UKB.v2.albuminuria.n382500.sorted.tsv.gz` for this assignment. This is the same GWAS we used for the in class computing example. The individuals included in this GWAS are primarily of European (and mostly British) ancestry.

We will focus on the region of chromosome 12 centered around rs2601006, the most significant SNP on that chromosome.

1)	Create a LocusZoom plot for rs2601006, with the most significant SNP as the LD reference SNP, between which all other LD is reported.   Include the LD legend in the plot, and display with European LD. Use a screenshot to record your plot, as the LocusZoom png download feature does not include LD reference or version of LocusZoom used.  Provide the plot here.
a. Interpret the plot:   do you think there may be more than one independent association in the region?  Why, or why not?
b. Do the variants that are highly associated with Albuminuria and in high LD with the top variant (in red) have similar effect sizes and allele frequencies to the top variant?
c. Do the variants  in modest LD with the top variant (in green) near 69.85MB have similar allele frequency and effect as the top variant?


2) a. Use GCTA to perform a conditional analysis on the full chromosome 12.  Use the 1000G Europeans as a reference sample.    Present the results, and interpret your findings:  do you have any concerns about the analysis? Is there evidence for multiple independent associations on this chromosome?  

b. Repeat the analysis, but use the 1000G Africans as a reference sample.  Present the results, and interpret your findings as if the UK Biobank data were an African ancestry sample.   do you have any concerns about the analysis? Is there evidence for multiple independent associations on this chromosome?  

c. Explain why it is important to match the LD panel ancestry with the GWAS participants ancestry, based on your analyses a) and b) above.

3)	If you were to do a full GWAS of this phenotype in a large African ancestry sample, would you expect to see more than one independent association on chromosome 12?  Why or why not?
- Answer: 

4)	Use VEP and/or or other tools to determine the functional annotations for rs2601006 and other SNPs that have p<5e-7 in within that region between and including positions 69,850,008 through 70,008,337  (this is the region of the Locus Zoom plot where there appear to be highly associated variants that are all in high LD).

a. How many variants were there in this region?
b. How many of the variants are intronic?
c. Are there any SNPs in the total set that are more likely to be functional or causal than others?  Present and explain your findings.
d. What does VEP tell you about the top SNP (lowest p-value) rs2601006?  Does this variant appear to have a known function?

