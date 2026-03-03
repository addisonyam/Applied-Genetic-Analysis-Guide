## Week 6: Meta-Analysis for Genetic Studies
### BS859 Applied Genetic Analysis
### Feburary 25, 2026

## PART 1: LECTURE CONTENT - CONCEPTUAL FOUNDATIONS

### 1. Why Meta-Analysis in Genetics?

#### 1.1 The Power Problem in GWAS
Genome-wide association studies test millions of SNPs, requiring a very stringent significance threshold to control false positives:
- Genome-wide significance: p < 5 × 10⁻⁸

This stringent threshold means:
- A single study may not have enough power to detect associated SNPs
- True associations with small-to-moderate effects can be missed
- Solution: Combine information across multiple studies to increase sample size and power

1.2 Two Approaches to Combine Data
Option 1: Pool individual-level data
- Combine raw data from all studies, analyze as one large sample
- Advantages: Maximum flexibility, consistent analysis methods
- Disadvantages:
    - Confidentiality restrictions may prevent data sharing
    - May not be desirable due to heterogeneity:
        - Different study designs requiring different statistical methods
        - Different population backgrounds (diet, environment, etc.)
        - Different trait measurements or definitions

Option 2: Meta-analysis ← Focus of this week
- Analyze each study separately, then combine summary statistics
- Advantages:
    - Can combine studies even without access to individual-level data
    - Preserves study-specific analyses that account for design differences
    - Widely used and accepted in genetics

### 2. The Meta-Analysis Framework
#### 2.1 Basic Setup
We have k studies (i = 1 to k), each testing the same null hypothesis:
- H₀: β = 0 (no association between SNP and trait)
Each study i provides:
- Effect estimate β̂ᵢ (e.g., regression coefficient from linear or logistic regression)
- Standard error SE(β̂ᵢ)
- P-value pᵢ
- Direction of effect (+ or -)
Goal: Obtain a single pooled estimate and test of H₀ using all k studies.

Key assumption: Independence of study results (each subject contributes to only one study)

2.2 First Step: Harmonize Effect Measures
All studies must report a consistent effect measure:
- For continuous traits: β coefficient from linear regression
- For binary traits: log odds ratio from logistic regression

These can be combined across studies because they share the same interpretation: change in trait per copy of the effect allele.

3. Fixed Effects Meta-Analysis
3.1 The Fixed Effects Assumption
Core assumption: All studies are estimating the same underlying population effect β

If any study had an infinitely large sample size, its effect estimate would exactly equal β

The only reason study estimates differ is random sampling error

Visual representation: All study estimates scatter around the same true value β, with spread determined by each study's precision.

3.2 Inverse Variance Weighting
The most common fixed effects approach weights each study by its precision:

Weight for study i: wᵢ ∝ 1 / SE(β̂ᵢ)²

Interpretation:

Larger studies → smaller SE → larger weight

Smaller studies → larger SE → smaller weight

Scaled weights (sum to 1):

wᵢ = [1/SE(β̂ᵢ)²] / [Σⱼ 1/SE(β̂ⱼ)²]

3.3 Pooled Effect Estimate and Standard Error
Combined effect estimate:

β̂ = Σ wᵢ β̂ᵢ

Standard error of combined estimate:

SE(β̂) = √[1 / Σ 1/SE(β̂ᵢ)²]

3.4 Test Statistic
Z_metaβ = β̂ / SE(β̂)

Under H₀: β = 0, Z_metaβ ~ N(0,1) (or Z_metaβ² ~ χ²₁)

