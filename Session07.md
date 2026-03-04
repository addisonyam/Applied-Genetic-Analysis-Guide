## Week 7: Post-GWAS Analyses - Conditional Analysis and Fine Mapping
### BS859 Applied Genetic Analysis
### March 4, 2026

- this exam is about how you interpret results and understanding of concepts

## PART 1: LECTURE CONTENT - CONCEPTUAL FOUNDATIONS

### 1. What is Fine Mapping?

#### 1.1 The Problem
GWAS identifies regions of the genome associated with a trait, but due to linkage disequilibrium (LD), there are often many variants with similar strength of association in an associated region.

Key insights:
- The variant with the smallest p-value is not typically the causal variant
- The causal variant may not even be genotyped or imputed, but may be in LD with the identified variants
- Multiple variants in a region could be independently associated with the trait

The goal of fine mapping: Determine which variants are most likely to be functional (causal) and quantify the strength of evidence.

#### 1.2 The Fine Mapping Pipeline
From initial GWAS to causal variant identification:

1. GWAS → Manhattan plot identifies genome-wide significant regions
2. Define regions of interest based on LD structure
3. Statistical fine mapping to identify independent signals and credible sets
4. Annotation of prioritized variants with genomic features
5. Functional follow-up studies

### 2. Understanding Associated Regions with LocusZoom
#### 2.1 What is LocusZoom?
LocusZoom is a visualization tool that displays:
- Association signals (-log10(p-values)) across a genomic region
- LD structure with a reference SNP (color-coded by r²)
- Gene annotations in the region
- Recombination rates

Web tool: https://my.locuszoom.org/

#### 2.2 Interpreting a LocusZoom Plot
Key elements:
- Y-axis: -log10(p-value) from your GWAS
- X-axis: Genomic position (chromosome coordinates)
- Color coding: LD (r²) with the reference SNP (usually the top SNP)
    - Red: High LD (r² > 0.8)
    - Orange/Yellow: Moderate LD (r² 0.5-0.8)
    - Green: Low LD (r² 0.2-0.5)
    - Blue: Very low LD (r² < 0.2)
- Bottom panel: Gene locations and orientations
- Recombination rate: Overlaid line showing local recombination hotspots

What to look for:
- Multiple distinct peaks suggest multiple independent signals
- A single broad peak suggests one signal with multiple correlated variants
- Variants in low LD but still significant suggest independent associations

#### 2.3 LocusZoom Workflow
1. Upload GWAS summary statistics (tab-delimited .txt or .gz)
2. Specify genome build (hg19/GRCh37 or hg38/GRCh38)
3. Define column mappings (chromosome, position, p-value, etc.)
4. Navigate to a region of interest
5. Select a reference SNP for LD coloring
6. Choose LD reference panel (e.g., 1000 Genomes EUR, AFR, etc.)

### 3. Identifying Independent Signals
#### 3.1 Why Identify Independent Signals?
- Multiple variants in a region may independently influence the trait
- Fine mapping should account for multiple causal variants
- Helps prioritize variants for functional studies
- Informs understanding of genetic architecture

#### 3.2 Forward Stepwise Regression (Individual-Level Data)
If you have individual-level data:

1. Test all SNPs for association (standard GWAS)
2. Take the most significant SNP, add it to model as a covariate
3. Test all remaining SNPs conditional on the first SNP
4. If a second SNP is significant, add it to model
5. Repeat until no more SNPs are significant at specified α-level

Limitation: Requires individual-level data from all studies.

#### 3.3 Stepwise Approach for Meta-Analysis
Without individual-level data:

1. Meta-analyst performs meta-analysis, identifies top SNP
2. Each study runs analysis conditional on top SNP, shares results
3. Meta-analyst performs conditional meta-analysis, identifies next top SNP
4. Repeat until no more significant signals

Limitation: Slow, tedious, requires multiple re-analyses by each study.

#### 3.4 Conditional Analysis with GWAS Summary Statistics (GCTA-COJO)
The key insight: Use LD from a reference sample (same population as GWAS) to approximate conditional analysis using only summary statistics.

How it works:

The joint model for n individuals and m SNPs is:
- y = Xb + e where y is phenotype, X is genotype matrix, b is joint effects

The marginal effects from single-SNP GWAS are:
- β̂ = D⁻¹X'y where D is diagonal of X'X

We can convert marginal to joint effects:
- X'y = Dβ̂
- b̂ = (X'X)⁻¹X'y = (X'X)⁻¹Dβ̂

Since we don't have X'X from meta-analysis, we estimate it from a reference sample with the same LD structure.

Requirements:
- GWAS summary statistics (SNP, effect allele, other allele, frequency, beta, se, p, N)
- LD reference panel from same population (PLINK format)
- Sufficient reference sample size (≥2000 recommended)

GCTA-COJO options:
| Option           | Purpose                                         |
|------------------|-------------------------------------------------|
| --cojo-file      | Input summary statistics file                   |
| --cojo-slct      | Stepwise selection to identify independent SNPs |
| --cojo-cond      | Conditional analysis on given SNP list          |
| --cojo-joint     | Estimate joint effects for all SNPs             |
| --cojo-p         | P-value threshold for selection (default 5e-8)  |
| --cojo-wind      | Window size for LD calculation (kb)             |
| --cojo-collinear | R² cutoff for collinearity (default 0.9)        |
| --maf            | Minor allele frequency threshold                |

#### 3.5 Limitations of Conditional Analysis
- Can't assume selected variants are causal
- Multiple signals could be due to untyped variant in LD with multiple SNPs
- Unlikely accurate when multiple causal SNPs are in high LD
- LD estimates for rare variants have wide confidence bounds
- Reference LD must match GWAS population

### 4. Fine Mapping Within a Region: Identifying a Causal Set
#### 4.1 Heuristic Approach
Select all SNPs within some LD threshold (e.g., r² > 0.8) of the lead SNP.

Limitations:
- Doesn't account for joint effects
- No objective measure of confidence
- May miss causal variants not in high LD with lead SNP

#### 4.2 Penalized Regression (Lasso, Elastic Net)
- Builds joint regression model for all SNPs in region
- Shrinks small effect estimates toward zero
- Performs variable selection

Best for: Prediction

Not ideal for: Identifying causal variants (will exclude one of two highly correlated causal SNPs)

#### 4.3 Bayesian Approaches
Core concept: Represent causal status with indicator vector c (1 = causal, 0 = non-causal). For m SNPs, there are 2ᵐ possible models.

Posterior probability of a model:

P(M_c | D) = P(D | M_c) × P(M_c) / Σ P(D | M) × P(M)

Where:
- P(D | M_c) is marginal likelihood (integrating over β)
- P(M_c) is prior probability of the model

Posterior Inclusion Probability (PIP):

PIPⱼ = Σ P(M | D) over all models where SNP j is causal

Advantages of PIP:
- Can be directly compared among SNPs
- Tends to select fewer SNPs than LD-based heuristics
- High PIP variants enriched for functional/pathogenic variants
- Better performance than stepwise or penalized regression in simulations

Factors affecting PIP:
- True effect size
- Sample size
- Number of SNPs in region
- LD structure (more correlation → poorer power)
- Prior model

#### 4.4 Bayesian Fine Mapping Software
From Schaid et al. 2018:
| Software       | Summary Stats? | Max Causal Variants | Search Method               |
|----------------|----------------|---------------------|-----------------------------|
| CAVIAR/eCAVIAR | Yes            | Fixed               | Exhaustive                  |
| FINEMAP        | Yes            | Fixed               | Shotgun stochastic search   |
| PAINTOR        | Yes            | Fixed/Computed      | Exhaustive/MCMC             |
| DAP            | Yes            | Computed            | DAP algorithm               |
| SuSiE          | Yes            | Computed            | Iterative Bayesian stepwise |
| JAM            | Yes            | Fixed/Computed      | Exhaustive/MCMC             |

#### 4.5 Population Diversity in Fine Mapping
Why multiple ancestries help:
- Different populations have different LD patterns
- Can separate effects of variants that are in LD in one population but not another
- Improves ability to pinpoint causal variants

### 5. Genomic Annotation
#### 5.1 What is Genomic Annotation?
Information about the biologic function of DNA sequence and genetic variants:

Protein coding:
- Loss-of-function > splice site > missense > synonymous > intronic
- Impact scores predicting deleteriousness (SIFT, PolyPhen, CADD)

Non-protein coding:
- Regulatory regions: promoters, enhancers
- Long non-coding RNA
- Transcription start sites
- Transcription factor binding sites
- Chromatin accessibility
- Histone modification patterns

Structural annotation:
- ORFs and gene structure
- Coding regions
- Location of regulatory motifs

#### 5.2 Annotation Resources
| Resource            | Description                                            |
|---------------------|--------------------------------------------------------|
| dbSNP               | Basic variant info, frequencies, clinical significance |
| Ensembl VEP         | Comprehensive variant annotation                       |
| UCSC Genome Browser | Multi-track visualization                              |
| GTEx                | eQTL expression data                                   |
| ENCODE              | Regulatory elements                                    |
| ClinVar             | Clinical significance                                  |
| GWAS Catalog        | Published associations                                 |

#### 5.3 Integrating Annotation into Fine Mapping
Ad hoc approach:
- Determine annotation for SNPs in credible set
- Prioritize for functional validation based on annotation

In models:
- Use annotation to define prior probabilities for variants
- May improve fine mapping when regions share enriched annotations

Caveat: Current understanding of genomic function may be too limited to accurately improve priors.

### 6. Important Considerations
#### 6.1 Multiple Testing in Conditional Analysis
- Stepwise selection identifies SNPs at p < threshold (usually 5e-8)
- Each conditional test is independent, so no additional multiple testing correction needed for selection
- But interpret with caution: selected SNPs are not necessarily causal

#### 6.2 LD Reference Panel Requirements
- Must come from same population as GWAS
- Sample size: ≥2000 minimum, 4000 optimal
- Clean reference sample: no cryptic relatedness
- Check allele frequency concordance between GWAS and reference

#### 6.3 Common Pitfalls
| Pitfall                              | Consequence               | Prevention           |
|--------------------------------------|---------------------------|----------------------|
| Using wrong LD reference             | False positives/negatives | Match ancestry       |
| Small reference sample               | Imprecise LD estimates    | Use ≥2000 samples    |
| Ignoring allele mismatches           | Incorrect strand/allele   | Check badsnps output |
| Interpreting selected SNPs as causal | Overconfidence            | Use Bayesian PIPs    |
| Rare variants in reference           | Unreliable LD             | Use MAF filter       |

#### 6.4 Final Thoughts
- Different software makes different assumptions
- Compare results across methods for robustness
- Consistent results increase confidence
- Statistical methods alone cannot determine causality
- Collaborate with lab scientists for functional validation

## PART 2: COMPUTING EXAMPLES - ALBUMINURIA GWAS
### 7. Albuminuria GWAS Dataset
#### 7.1 Study Information
- Phenotype: Albuminuria (excretion of albumin in urine)
- Sample size: 382,500 individuals (UK Biobank)
- Population: Primarily European (British) ancestry
- Source: Haas ME et al. (2018), AJHG 103, 461-473
- Data location: `/projectnb/bs859/data/meta/downloads/UKB.v2.albuminuria.n382500.sorted.tsv.gz`

#### 7.2 File Format
```text
SNP rs_marker_id Chr Position_hg19 Effect_allele Other_allele Effect_AF Beta SE Pvalue Sample_size
1:55326_T_C rs3107975 1 55326 C T 0.00246275 0.0138664 0.0164503 0.399268 382500
1:86028_T_C rs114608975 1 86028 C T 0.0423959 -0.000363436 0.00468298 0.93814 382500
```

Columns:
1. SNP: UK Biobank variant ID
2. rs_marker_id: rsID or reference panel ID
3. Chr: Chromosome
4. Position_hg19: Position in GRCh37 coordinates
5. Effect_allele: Effect allele (coded)
6. Other_allele: Non-effect allele
7. Effect_AF: Frequency of effect allele
8. Beta: Effect size (log(mg urine albumin/g urine creatinine))
9. SE: Standard error of beta
10. Pvalue: Association p-value
11. Sample_size: N for this variant

#### 7.3 Exploring the Data
```bash
# Count variants
zcat UKB.v2.albuminuria.n382500.sorted.tsv.gz | wc
# 11,709,858 variants

# Extract genome-wide significant hits (p < 5e-8)
zcat UKB.v2.albuminuria.n382500.sorted.tsv.gz | awk 'NR==1 || $10<5e-8 {print $0}' > tophits.txt
# 1,682 variants at 33 loci

# Extract highly significant hits (p < 5e-50) - CUBN region
zcat UKB.v2.albuminuria.n382500.sorted.tsv.gz | awk 'NR==1 || $10<5e-50 {print $0}' > p5e-50.txt
```

### 8. LocusZoom Visualization
#### 8.1 Accessing the Pre-Processed Results
The full GWAS was uploaded to LocusZoom and is available at:

https://my.locuszoom.org/gwas/359645/?

token=0d1adaeabd644ee0aba87c221cbfa3cb

#### 8.2 Key Observations in CUBN Region (Chromosome 10)
- Top SNP: rs141640975 (p = 1.75e-107)
- Effect allele frequency: 0.0026 (very rare)
- Many other variants also strongly associated
- Red variants have similar allele frequency and high LD in EUR 1000G
- Some variants far from top SNP, in low LD, but still genome-wide significant → possible independent signals

#### 8.3 Checking LD Reference
LocusZoom uses 1000 Genomes LD reference:
- For common variants: reliable
- For rare variants: may not be present or LD estimates imprecise due to small sample size (503 EUR individuals)

### 9. GCTA Conditional Analysis
#### 9.1 Prepare Summary Statistics File
GCTA requires input format:
```text
SNP A1 A2 freq b se p N
rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850
```

Reformat the albuminuria data:

```bash
zcat UKB.v2.albuminuria.n382500.sorted.tsv.gz | \
  awk '{print $2,$5,$6,$7,$8,$9,$10,$11}' > togcta.txt

# Check result
head togcta.txt
```
#### 9.2 LD Reference: 1000 Genomes EUR
```bash
# Location of 1000G EUR PLINK files
ls /projectnb/bs859/data/1000G/plinkformat/

# Check sample size
wc /projectnb/bs859/data/1000G/plinkformat/1000G_EUR.fam
# 503 individuals

# Check if top variant exists in reference
grep rs141640975 /projectnb/bs859/data/1000G/plinkformat/1000G_EUR.bim
# 10 rs141640975 0 16992011 A G

# Check allele frequency in reference
plink --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR \
  --chr 10 --freq --out chr10freq --allow-extra-chr
grep rs141640975 chr10freq.frq
# Frequency = 0.003976 (4 copies in 503 individuals)
```

#### 9.3 Run GCTA Conditional Analysis
```bash
# Load modules
module load gcta/1.94.1
module load plink/1.90b6.21

# Basic conditional analysis (all variants)
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR \
  --cojo-file togcta.txt \
  --cojo-slct \
  --chr 10 \
  --out chr10 > chr10.log

# With MAF filter (MAF ≥ 0.01)
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR \
  --cojo-file togcta.txt \
  --cojo-slct \
  --chr 10 \
  --maf 0.01 \
  --out chr10.maf01 > chr10.maf01.log

# Wrong LD reference (AFR instead of EUR)
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_AFR \
  --cojo-file togcta.txt \
  --cojo-slct \
  --chr 10 \
  --out chr10.AFRld > chr10.AFRld.log
```

#### 9.4 Understanding GCTA Output Files
| File               | Description                                            |
|--------------------|--------------------------------------------------------|
| chr10.jma.cojo     | Results of stepwise selection (joint model)            |
| chr10.ldr.cojo     | LD correlation matrix between selected SNPs            |
| chr10.badsnps      | SNPs with allele mismatches between GWAS and reference |
| chr10.freq.badsnps | SNPs with frequency mismatches                         |
| chr10.log          | Screen output with progress and warnings               |

#### 9.5 Interpreting .jma.cojo Output
Columns in joint model output:
| Col   | Name          | Description                          |
|-------|---------------|--------------------------------------|
| 1     | Chr           | Chromosome                           |
| 2     | SNP           | SNP identifier                       |
| 3     | bp            | Physical position                    |
| 4     | freq          | Frequency of effect allele in GWAS   |
| 5     | A1            | Effect allele                        |
| 6-8   | b, se, p      | Original GWAS effect, SE, p-value    |
| 9     | n             | Effective sample size                |
| 10    | freq_geno     | Frequency in reference sample        |
| 11-13 | bJ, bJ_se, pJ | Joint analysis effects, SE, p-value  |
| 14    | LD_r          | LD correlation with next SNP in list |

#### 9.6 Results Without MAF Filter
- 33 SNPs in final joint model
- Includes very rare variants (MAF < 0.01)
- Top SNP rs141640975 (MAF=0.0026) is selected

#### 9.7 Results With MAF ≥ 0.01 Filter
- 7 independent associations
- Top SNP rs141640975 excluded (MAF too low)
- Some variants not genome-wide significant in single-SNP GWAS become significant in joint model
- Last two variants are in other chromosome 10 regions (not CUBN)

#### 9.8 Results With Wrong LD Reference (AFR)
- Different LD structure leads to different selection
- Demonstrates importance of matching LD panel ancestry

#### 9.9 Checking for Problems
Always check log file and badsnps files:

```bash
# Count allele mismatches
wc chr10.badsnps
# 1116 SNPs with allele mismatches

# Check frequency mismatches
head chr10.freq.badsnps
```

Frequency mismatches may indicate:
- Strand issues
- Different allele frequencies between GWAS and reference
- Population differences
- Need for investigation if these SNPs are associated

### 10. Annotation with VEP
#### 10.1 What is VEP?
Variant Effect Predictor (VEP) from Ensembl:
- Annotates variants with genomic features
- Predicts functional consequences
- Provides gene, transcript, protein changes
- Includes population frequencies
- Links to known clinical significance

Web tool: http://www.ensembl.org/Homo_sapiens/Tools/VEP

#### 10.2 Preparing Input for VEP
```bash
# Extract rsIDs from top variants (p < 5e-50)
cut -f2 p5e-50.txt > toVEP.txt
# 18 variants in CUBN region
```
#### 10.3 VEP Input Options
- Paste data: Directly enter variant identifiers
- Upload file: Upload list of rsIDs or VCF
- Provide URL: Link to file

#### 10.4 VEP Output Information
For each variant, VEP provides:
- Gene: Ensembl gene ID and symbol
- Consequence: e.g., missense, intronic, downstream
- cDNA position: Position in transcript
- CDS position: Position in coding sequence
- Protein position: Amino acid position
- Amino acid change: e.g., R/W (Arg to Trp)
- SIFT score: Deleteriousness prediction
- PolyPhen score: Functional impact prediction
- Existing variation: Known rsIDs
- ClinVar: Clinical significance if available

#### 10.5 Using VEP Results
- Identify which variants are protein-coding vs. non-coding
- Check for missense, nonsense, splice site variants
- Look for high-impact consequences
- Cross-reference with ClinVar for known associations
- Prioritize variants for functional follow-up

### 11. Summary: Post-GWAS Analysis Workflow
##### Step 1: Identify Genome-Wide Significant Loci
- Extract variants with p < 5e-8
- Group into LD-based regions

##### Step 2: Visualize with LocusZoom
- Examine LD structure
- Identify potential independent signals
- Note variants with different LD patterns

##### Step 3: Conditional Analysis
- Use GCTA-COJO with appropriate LD reference
- Identify independent association signals
- Check for allele/frequency mismatches

##### Step 4: Fine Mapping
- Apply Bayesian methods (FINEMAP, SuSiE, CAVIAR)
- Calculate posterior inclusion probabilities
- Define credible sets of likely causal variants

##### Step 5: Annotation
- Annotate variants with VEP or similar tools
- Prioritize functional variants
- Check clinical significance databases

##### Step 6: Functional Follow-Up
- Collaborate with experimentalists
- Validate in relevant cell types/tissues
- Test in additional populations

### 12. Key Takeaways for Homework
##### For Question 1 (LocusZoom)
- Look at LD structure with top SNP
- Check for multiple distinct peaks
- Compare effect sizes and frequencies of variants in different LD bins
- Variants in high LD should have similar effect sizes and frequencies

##### For Question 2 (GCTA Conditional Analysis)
- Use 1000G EUR as reference (matches GWAS ancestry)
- Run with and without MAF filter
- Compare results with AFR reference
- Interpret: more independent signals with wrong reference?

##### For Question 3 (African Ancestry Expectations)
- Different LD structure in African populations
- Likely to see more independent signals (smaller LD blocks)
- May fine-map to smaller credible sets

##### For Question 4 (VEP Annotation)
- Count variants in region (69.85-70.01 Mb)
- Identify intronic vs. exonic variants
- Look for potentially functional variants (missense, splice site)
- Check top SNP rs2601006 annotation

### 13. Quick Reference: Key Commands for Homework
##### GCTA Conditional Analysis
```bash
# Load modules
module load gcta/1.94.1
module load plink/1.90b6.21

# Prepare summary stats (if needed)
# Your file should have: SNP A1 A2 freq b se p N

# Run with EUR reference
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR \
  --cojo-file your_summary_stats.txt \
  --cojo-slct \
  --chr 12 \
  --out chr12_EUR

# Run with AFR reference
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_AFR \
  --cojo-file your_summary_stats.txt \
  --cojo-slct \
  --chr 12 \
  --out chr12_AFR
```
##### VEP Annotation
- Go to https://www.ensembl.org/Homo_sapiens/Tools/VEP
- Upload list of rsIDs
- Select appropriate options
- Download results

**With these notes, you should be able to:**
- Understand the goals and methods of fine mapping
- Interpret LocusZoom plots
- Run and interpret GCTA conditional analyses
- Recognize the importance of matching LD reference ancestry
- Use VEP for variant annotation
- Complete Homework 7 without referring back to slides

