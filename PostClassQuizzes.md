## Post-Class 1: Genotype Data QC & PLINK
### BS859 Applied Genetic Analysis
### January 21, 2026
1. Which statement about call rate filtering is correct? Select up to 4 options
    - You should always remove low call-rate individuals first to maximize SNP retention.
    - **Low SNP call rate means many individuals lack a genotype for that marker.**
    - Call rate refers only to the number of SNPs on the chip, not to samples(individuals).
    - Call rate filtering is irrelevant for modern chips and can be skipped.

Feedback: SNP call rate = proportion of samples with a non-missing genotype at that marker. Call rate applies to both markers and samples.  While modern chips tend to have higher call rates, filtering out low call rate is still recommended.

2. In a case-control GWAS, which sample should you generally use when testing SNPs for departure from Hardy–Weinberg equilibrium (HWE) as a QC measure?
All individuals (cases + controls) always
    - Only affected individuals (cases)
    - **Only unaffected individuals (controls)**
    - Only related individuals

Feedback: Controls are used for HWE QC because cases may deviate due to true association; related individuals complicate HWE testing.

3. Because males have one X chromosome, X-chromosome SNPs outside the pseudo-autosomal region are expected to appear homozygous in males; therefore a high proportion of heterozygous calls on X in a male sample suggests a problem (e.g., assay error or sample mix-up).
    - **True**
    - False

Feedback: Males are hemizygous for non-PAR X SNPs, so heterozygous genotypes are unexpected and indicate possible issues.

4. Identical-by-state (IBS) is always less than identical-by-descent (IBD) for a pair of individuals.
    - True
    - **False**

Feedback: IBS counts alleles identical in state regardless of ancestry and is always greater than or equal to IBD (IBS ≥ IBD).

5. Why do we perform LD-pruning (creating a set of variants with low pairwise linkage disequilibrium) before some analyses?
    - To increase the number of variants for association testing
    - To remove monomorphic markers only
    - **To obtain an approximately independent set of markers for analyses like heterozygosity estimation, IBD estimation, and PCA**
    - To force all SNPs to have the same minor allele frequency

Feedback: LD-pruning reduces correlation among markers so methods that assume independence perform correctly.

6. A pair of samples has PI_HAT ≈ 0.5 in a PLINK IBD output. Which is the most likely relationship?
    - Unrelated
    - Monozygotic twins or duplicate sample
    - **Parent–child or full siblings (first-degree)**
    - Third-degree relatives (first cousins)

Feedback: PI_HAT ≈ 0.5 corresponds to sharing on average half their alleles IBD, consistent with first-degree relationships.

7. Which are indicators that a DNA sample might be contaminated (not a pure sample from one individual)?  Choose all that are correct.
    - **The individual has high IBD sharing with multiple other individuals in the - sample**
    - **The individual has an F statistic that is negative and large in absolute value**
    - The sample has a low call rate (high missing rate)
    - The sample has 0 IBD sharing with all other individuals in the sample
    - The individual has an F statistics that is positive and large in absolute value

Feedback: Contamination can induce false heterozygote calls (excess heterozygosity, negative F) and artificially inflate IBS/IBD with multiple samples.

## Post-Class Quiz 2: Population Structure Detection
### BS859 Applied Genetic Analysis
### January 28, 2026
1. Which definition best describes population structure in genetic studies?
    - A population where every individual has identical genotypes
    - **Deviation from random mating, producing subpopulations with different allele frequencies and correlations among variants**
    - Complete panmaxia across all geographic regions
    - A dataset with no missing genotypes

Feedback: Population structure occurs when mating is non-random (geographic separation, selection, inbreeding, assortative mating), producing subgroups with different genotype frequencies.

2. Population structure can cause spurious associations in GWAS when which of the following is true? Select up to 2 options
    - The study uses only unrelated individuals from a single homogenous population
    - Subpopulations differ in trait distributions but have identical genotype frequencies
    - **Subpopulations differ in both allele frequencies and trait distributions**
    - Subpopulations differ in allele frequencies but have identical trait distributions

Feedback: Spurious associations arise when ancestry (population structure) is associated with both genotype and phenotype. If both vary across subgroups, genotype-phenotype associations can be confounded.

3. Principal Components Analsysis (PCA) on the full set of genome-wide markers is a common and effective method to detect population structure in GWAS data. 
    - True
    - **False**

Feedback: PCA summarizes genetic variation into orthogonal axes (PCs) that often reflect ancestry; However, we do not use the full set of markers.  LD pruning is used to remove highly correlated markers and reduce correlation among markers before PCA.

4. When using PCA to infer broad ancestry (e.g., European vs. East Asian vs. African), which approach is recommended?
    - Run PCA on only you study damples without external references
    - Use only A/T and G/C SNPs for the PCA
    - Remove all SNPs with MAF > 0.05
    - **Include reference samples with known population labels (e.g., HapMap/1000G) in the PCA so study samples can be compared to known populations**

Feedback: Including reference samples lets you place your samples relative to known population clusters and identify ancestry or outliers

5. Why are A/T and G/C SNPs often removed before merging datasets from different sources for PCA? 
    - They have the highest allele frequencies
    - They are always multiallelic
    - **Their strand orientation is ambiguous (complements are the same), risking incorrect matching across datasets**
    - They are always on the sex chromosomes

Feedback: A/T and G/C SNPs look identical after strand flip (A<->T, G<->C), making it hard to determine correct allele orientation when merging datasets; excluding them avoids strand-mismatch errors.

## Post-Class Quiz 3: Common Varitant GWAS - Association Analysis & Mixed Models
### BS859 Applied Genetic Analysis
### February 4, 2026
1. In a standard GWAS using an additive genotype coding for a biallelic SNP (alleles A and C), with coded allele A and reference allele C, what does the regression coefficient β for the coded genotype represent in a logistic regression model?
    - The difference in mean phenotype between AA and CC
    - **The log of the odds ratio for disease for a genotype with one additional copy of the A allele compared to the reference genotype**
    - The odds ratio comparing heterozygotes to homozygous reference only
    - The variance in phenotype explained by the SNP

Feedback: In logistic regression with additive coding (0,1,2), the genotype coefficient β is the change in log odds per additional copy of the coded (effect) allele. exp(β) is the odds ratio for a one-allele increase.

2. When deciding whether to include principal components (PCs) of genetic data as covariates in a GWAS, which statement is correct?
    - You must always include the top 10 PCs to avoid any confounding
    - PCs are irrelevant if you use PLINK for association testing 
    - Including more PCs always increases power in logistic regression
    - **PCs need to be included if they are associated with the phenotype of interest.**

Feedback: If PCs are unrelated to the phenotype, they are unnecessary for avoiding confounding. One common practice is to adjust for PCs that are associated with the trait (or to include a set number of top PCs for conservative control). Including many PCs unnecessarily can reduce power in logistic regression.

3. The genomic control inflation factor λ<sub>GC</sub> is computed from the ratio of:
    - **Observed median test statistic to expected median test statistic under the null**
    - Observed mean p-value to 0.05
    - Number of significant SNPs to total SNPs tested
    - Observed genomic heritability to total phenotypic variance

Feedback: λ<sub>GC</sub>= (median observed χ2) / (median expected χ2 under null). It quantifies overall inflation of test statistics relative to null expectation.

4. In GWAS with millions of correlated tests, the conventional genome-wide significance threshold ~5×10^-8 is based on a false discovery rate of 0.05
    - True
    - **False**

Feedback: The conventional genome-wide threshold (~5×10^-8) comes from Bonferroni correction for roughly 1e6 independent tests (0.05/1e6). This is a pragmatic approximation given linkage disequilibrium between SNPs.

5. Which of the following best describes the purpose of QQ plots in GWAS?
    - To display genomic positions vs -log10(p) across the genome
    - To estimate per-SNP effect sizes
    - **To compare observed p-value quantiles to expected under the null to detect inflation and/or enrichment of associations**
    - To compare observed p-value quantiles to expected under the null to detect inflation or enrichment of associations

Feedback: QQ plots show observed vs expected p-value quantiles; deviations from the diagonal in the middle of the distribution can indicate inflation (population structure, bias) or true signal enrichment at the tail.

6. Which advantage of generalized linear mixed models (GLMMs) like those implemented in GMMAT or GENESIS is most relevant when analyzing binary traits in GWAS with population structure and distantly related individuals?
    - GLMMs avoid the need to compute principal components entirely.
    - **GLMMs model a random effect based on the genetic relationship matrix (GRM) to account for both cryptic relatedness and population structure while allowing a logistic link for binary outcomes.**
    - GLMMs always produce smaller p-values than standard logistic regression, increasing power in all settings.
    - GLMMs require no assumptions about the distribution of the outcome.

Feedback: GLMMs incorporate a random effect whose covariance is proportional to the GRM, which captures genetic similarity between individuals. This allows proper control for both relatedness and some forms of population structure in analyses of binary traits by combining the logistic (or other) link with random effects. They do not make PCs unnecessary, do not universally reduce p-values, and still rely on appropriate outcome distribution assumptions (e.g., Bernoulli for binary traits with a logistic link).


## Post-Class Quiz 4: Genotype imputation and analysis of imputed genotypes
### BS859 Applied Genetic Analysis
### February 11, 2026

1. Which reference panel attribute is generally most beneficial for imputation accuracy, especially for low-frequency variants?
- Smaller size but single-population ancestry match
- Using only HapMap sample
- Only including unrelated individuals with missing phenotype data
- **Larger size and diverse ancestries (including the study’s ancestry where possible)**

Feedback: Larger panels increase the chance of observing haplotypes that match study samples; including diverse ancestries helps, particularly for admixed samples and low-frequency variants. Small or older panels (HapMap) are less informative.

2. For very rare variants (very low MAF), the imputation R2 is a highly reliable indicator of the true imputation accuracy.
- True
- **False**

Feedback: At very low MAF, R2 is less reliable and can become noisy or biased; imputation accuracy is harder to estimate and often lower for rare variants.

3. Which of these statements about the imputation R2 (or INFO) metric is correct? (choose all that are correct)
- R2 is a per-individual measure of imputation certainty.
- **R2 estimates the squared correlation between imputed dosages and true genotypes and reflects the loss of information due to imputation; it becomes less reliable at very low MAF.**
- R2 is the exact squared correlation between imputed and true genotypes and is unbiased at all allele frequencies.
- R2 is only reported for genotyped (typed) variants, not imputed ones.

Feedback: R2 is an estimate of the squared correlation (imputed vs true) and indicates effective sample size (N * R2). It is less reliable at very low MAF. It is a per-variant metric, not per-individual, and is reported for imputed variants.

4. The imputed value commonly used as the predictor in regression (e.g., GWAS) is:
- The hard-called (best guess) genotype
- The imputation R2 metric
- The haplotype switch error rate
- **The allele dosage (expected allele count) derived from posterior genotype probabilities**

Feedback: Dosage (expected number of alternate alleles, e.g., 0–2) preserves imputation uncertainty and is typically used in regression. Best-guess GT loses uncertainty; R2 is a quality metric, not a predictor.

5. Suppose you imputed a cohort of N = 5,000 samples and observe an imputed SNP with R2 = 0.2 and MAF = 0.02. Approximately what is the effective sample size for association at this SNP, and what are the practical consequences?
- Effective N ≈ 5,000; include the SNP without concern. 
- Effective N ≈ 100; variant is essentially uninformative.
- **Effective N ≈ 1,000; power is reduced and the SNP may be excluded depending on study aims.**
- Effective N cannot be estimated from R2.

Feedback: Effective N ≈ N * R2 = 5,000 * 0.2 = 1,000. This reduced effective sample size means lower power to detect association; many analyses use R2 thresholds (e.g., 0.3) to exclude such variants, though for specific hypotheses one might keep them with caution.

6. Which of the following most directly improves the ability to meta-analyze studies that used different genotyping chips?
- Increasing per-sample sequencing depth
- Using principal components as covariates
- **Genotype imputation to a common dense reference panel**
- Applying stricter Hardy–Weinberg filters

Feedback: Different chips assay different variant sets. Imputing all studies to the same reference panel fills in variants and creates a common variant set for meta-analysis. PCs and QC are important but do not bridge differing chip content.

## Post-Class Quiz 5: Association analysis of Rare Variants and Sequence Data
### BS859 Applied Genetic Analysis
### February 18, 2026

1. Which of these is a recommended way to reduce multiple-testing burden when performing many variant-set tests per gene (different variant filters/weights)?
- Apply Bonferroni across every individual test without exception
- Only test single variants and ignore gene-based tests
- **Combine multiple test p-values per gene using ACAT to produce a single p-value**
- Use Ti/Tv ratio as a correction factor

2. Which of the following is the primary rationale for aggregating rare variants in gene- or region-based tests instead of testing each rare variant individually?

- Sequencing errors make single-variant tests impossible
- **Single-variant tests have low power for very low MAF variants**
- Aggregation removes the need for covariate adjustment
- Aggregation guarantees causal variants will be identified

3. A gene-set whose cumulative minor allele count (cMAC) is extremely low should still always be included in the multiple-testing correction because it might become significant after Bonferroni adjustment.
- True
- **False**

4. Which test is most appropriate when you expect many variants in a set to have effects in different directions (some increasing, some decreasing the trait)?
- Burden test
- CAST
- **SKAT (variance-component test)**
- Madsen-Browning burden

5. When designing gene-based rare-variant association tests, which variant inclusion strategy is most likely to increase statistical power by enriching the tested set for variants that are plausibly functional (i.e., more likely to affect protein function)?
- **Include all variants within a 1 Mb window around the gene, without additional filtering**
- Include only loss-of-function (LoF) variants and nonsynonymous variants predicted to be deleterious by functional annotation
- Include only intergenic variants with no predicted regulatory function
- Randomly select variants located within the gene region regardless of annotation

6. Which of the following is a proper reason to use a mixed-model framework (GLMM) for rare-variant set tests like SMMAT?
- **To account for relatedness or cryptic population structure via a genetic relationship matrix (GRM)**
- To remove the need to annotate variants by function
- Because mixed models give exact p-values without asymptotics
- To automatically increase the minor allele counts of rare variants

## Post-Class Quiz 6: Meta-Analysis for genetic studies
### BS859 Applied Genetic Analysis
### February 25, 2026

1. Which statement about the weighted (signed) Z-score meta-analysis is TRUE?

- It always gives the same result as inverse-variance pooling, regardless of study design
- It produces a pooled effect size estimate with standard error
- It ignores the direction of effect
**It does not require effect size estimates on the same scale across studies and combines signed Z-values weighted by sample size**

2. For meta-analysis of gene-based rare-variant SKAT tests, each study must provide:

- Gene-level p-values and sample sizes for each gene
- **Allele frequencies, single-variant score statistics, and the genotype covariance (Φ) matrix capturing LD structure**
- Single-variant effect estimates (β) and standard errors for each variant in the gene
Only burden test statistics for each gene from each study

3. Which of the following actions is required before combining effect estimates across studies?

- Ensure that all studies use identical genotyping platforms
- Standardize all effect sizes to have mean zero and variance one
- **Re-align effect estimates so they correspond to the same effect (counted) allele and correctly handle strand orientation (especially for A/T and C/G SNPs)**
- Exclude all SNPs showing heterogeneity before meta-analysis

4. By default, which meta-analysis method does METAL use when no scheme is specified in the control file?

- Fixed-effect inverse-variance weighting
- **Sample-size weighted signed Z-score meta-analysis**
- Random-effects REML meta-analysis
- Gene-based SKAT meta-analysis

5. Under a fixed-effects meta-analysis model, as study sample sizes increase, the study-specific effect estimates would converge to the same underlying true effect.

- **True**
- False

6. Which of the following statements about genomic-control (GC) correction in METAL is TRUE?

- **GC correction rescales test statistics using a genomic inflation factor (λ) and generally reduces inflated statistics, making p-values less significant**
- GC correction always increases the significance of strongly associated SNPs
- GC correction is unnecessary in meta-analysis because combining studies removes population stratification
- GC correction alters allele frequencies and effect estimates in the output

## Post-Class Quiz 7: Post-GWAS Analyses - Conditional Analysis and Fine Mapping
### BS859 Applied Genetic Analysis
### March 4, 2026

Question 1 - To approximate joint SNP effects from single-SNP marginal effects and a reference LD panel, which matrix is required or estimated?
- The environmental covariance matrix
- **X'X (the genotype variance-covariance matrix) estimated using allele frequencies and LD from a reference sample**
- The transcriptome co-expression matrix
- The inverse Fisher information matrix from unrelated traits

Question 2 - Which statement correctly contrasts the heuristic LD-threshold approach with penalized regression (e.g., LASSO) for fine mapping?
- Heuristic approaches estimate joint effects while penalized regression ignores multicollinearity
- Both methods are optimized for finding multiple causal variants in high-LD regions and always return identical sets
- **Heuristic uses pairwise LD thresholds around lead SNP but does not model joint effects; penalized regression models SNPs jointly and performs shrinkage/selection**
- Penalized regression never excludes causal variants

Question 3 - In Bayesian fine-mapping, what does PIP represent?
- **The posterior probability that a SNP is included in a causal model**
- The genome-wide adjusted p-value for a SNP
- The probability that a SNP is in high LD with the lead variant
- The estimated effect size after shrinkage

Question 4 - GCTA COJO conditional analysis using GWAS summary statistics does not require an external LD reference panel.
- True
- **False**

Question 5 - How can functional genomic annotation be incorporated into statistical fine-mapping?
- **By incorporating annotation into the prior model or using it to prioritize variants within a credible set**
- By directly estimating SNP–SNP correlations within the associated region
- By replacing the need for conditional or joint modeling
- By restricting analysis only to protein-coding variants

Question 6 - When annotating a list of SNPs with Ensembl VEP, which of the following is true?
- VEP reports a single summary consequence per variant and ignores transcript-specific differences
- VEP estimates linkage disequilibrium between nearby variants during annotation
- VEP performs Bayesian fine-mapping and calculates posterior inclusion probabilities
- **VEP can return multiple transcript-specific consequences per variant, including predicted functional impacts**