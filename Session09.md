## Week 9: Post GWAS analyses: Heritability 
### BS859 Applied Genetic Analysis
### March 25, 2026

Week 9: Post-GWAS Heritability Estimation
Comprehensive Study Notes – No Slides Needed
1. What is Heritability and Why Do We Estimate It?
After GWAS identifies associated SNPs and regions, we often want to know:

How much of the trait variation is due to genetics? (heritability)

Do two traits share genetic influences? (genetic correlation)

Which functional categories or cell types contribute to heritability? (partitioned heritability)

Heritability is the proportion of total phenotypic variance that is due to genetic factors.
It can be estimated from family data, individual-level GWAS data, or summary statistics.

2. Definitions: Variance Components
For a quantitative trait, we can decompose phenotypic variance:

text
σ²_T = σ²_G + σ²_E
where:

σ²_T = total phenotypic variance

σ²_G = genetic variance (additive + dominance)

σ²_E = environmental variance (including error)

Broad‑sense heritability (H²): proportion of variance explained by all genetic effects
H² = (σ²_A + σ²_D) / σ²_T

Narrow‑sense heritability (h²): proportion explained by additive genetic effects only
h² = σ²_A / σ²_T

For most GWAS, we focus on narrow‑sense heritability because additive effects are what SNPs capture.

3. Additive and Dominance Effects
For a single bi‑allelic locus (alleles B and b, with BB, Bb, bb):

a = half the difference between the two homozygote means (BB – bb)/2

d = deviation of the heterozygote mean from the midpoint of the two homozygotes

Genetic models:

Additive: d = 0

Dominant: d = +a

Recessive: d = –a

The degree of dominance = d/a.

4. Estimating Heritability from Relatives
Traditional methods use expected genetic resemblance among relatives.

Relationship	Additive coefficient (φ)	Dominance coefficient
MZ twins	1	1
Parent-offspring	1/2	0
Full sibs	1/2	1/4
DZ twins	1/2	1/4
Half-sibs	1/4	0
The covariance between relatives is a function of these coefficients and the variance components.
Example: Twin studies estimate h² ≈ 2(r_MZ – r_DZ) (r = intraclass correlation).

5. Estimating Heritability from Individual-Level GWAS Data (GCTA)
Instead of using known pedigree relationships, GCTA uses empirical genetic relationships estimated from all SNPs (the Genetic Relationship Matrix, GRM).

Mixed Model
For N individuals and M SNPs:

Standardize genotypes: X_{ij} = (g_{ij} – 2p_j) / √[2p_j(1–p_j)]

GRM: A = XXᵀ / M (each entry is the genetic relatedness between two individuals)

Model:

text
y = Xβ + ε
Treat β (SNP effects) as random: β ~ N(0, (h²/M) I)
Residuals: ε ~ N(0, (1–h²) I)

Then:

text
Var(y) = h² A + (1–h²) I
The heritability h² is estimated by restricted maximum likelihood (REML). This is implemented in GCTA (Yang et al. 2011).

For Dichotomous Traits
We can treat case–control status as a quantitative trait (0/1) but the heritability on this observed scale is not directly interpretable.
We convert to liability scale (assuming an underlying normal distribution with a threshold).
If K = population prevalence, p = proportion of cases in the sample, and Z = height of the standard normal density at the threshold t(K), then:

Observed scale h²_obs → liability scale h²_liab:

text
h²_liab = h²_obs × [K(1–K)]² / [p(1–p) Z²]
In GCTA, you can provide --pop-prev and --samp-prev to output liability-scale heritability.

6. LD Score Regression (LDSC)
LDSC estimates heritability from GWAS summary statistics without individual-level data.

Key Insight
Under a polygenic model (many small effects), the expected χ² statistic for a SNP j is:

text
E[χ²_j] = (N h² / M) ℓ_j + 1 + N α
where:

ℓ_j = LD score = Σ_k r²_{jk} (sum of squared LD with all other SNPs)

N = sample size

M = number of SNPs

α = inflation due to population stratification / cryptic relatedness

The intercept (1 + Nα) should be ≈ 1 if no stratification; if > 1, indicates inflation.

Heritability estimate:
From the slope b₁ = N h² / M → h² = (M/N) b₁.

Practical Steps (using LDSC software)
Format summary statistics using munge_sumstats.py:

Input: GWAS summary file (SNP, effect allele, other allele, beta/OR, SE, P, N, etc.)

Merge with a list of HapMap3 SNPs (to restrict analysis)

Output: a .sumstats.gz file with standardized columns.

Run heritability estimation:

text
ldsc.py --h2 <sumstats.gz> \
        --ref-ld-chr <ldscore_dir> \
        --w-ld-chr <ldscore_dir> \
        --out <output>
--ref-ld-chr: LD scores for the SNPs (e.g., from 1000 Genomes European population)

--w-ld-chr: LD scores used for regression weights (usually the same)

For binary traits, add --pop-prev K and --samp-prev p to obtain liability‑scale h².

Example: Alzheimer’s Disease (IGAP 2019)
Cases: 21,982, controls: 41,944 → p = 0.344

Population prevalence K assumed = 0.10

LDSC with 1000G EUR LD scores yields:

Observed scale h² = 0.0713 (SE 0.0114)

Liability scale h² = 0.0831 (SE 0.0133)

7. Genetic Correlation from Summary Statistics
The same principle extends to two traits. For two studies with sample sizes N₁ and N₂ and overlapping samples N_s, the expected product of Z‑scores is:

text
E[Z₁j Z₂j] = (√(N₁N₂) ρ_g / M) ℓ_j + (N_s ρ) / √(N₁N₂)
where:

ρ_g = genetic correlation between the two traits

ρ = phenotypic correlation among overlapping samples

The slope gives an estimate of ρ_g.

LDSC command:

text
ldsc.py --rg <sumstats1.gz>,<sumstats2.gz> \
        --ref-ld-chr <ldscore_dir> \
        --w-ld-chr <ldscore_dir> \
        --out <output>
Interpretation:

ρ_g = 0: no shared genetic architecture

ρ_g = 1: complete positive genetic correlation

ρ_g = –1: complete negative correlation

Example: AD vs. brain volume (height‑corrected) gave ρ_g = –0.063 (p = 0.53) – not significant.

8. Partitioning Heritability
LDSC can estimate heritability within functional annotations (e.g., coding, intronic, regulatory regions) by using annotated LD scores.

Each SNP is assigned a binary indicator for each annotation (e.g., 1 if in coding region, else 0).

LD scores are computed separately for each annotation, summing r² only over SNPs that also have the annotation.

Regression is done with multiple LD scores simultaneously, estimating the heritability contribution of each annotation.

Enrichment = (proportion of h² in annotation) / (proportion of SNPs in annotation)
If enrichment > 1, the annotation contributes more heritability than expected by chance.

Baseline model: 53 overlapping annotations (from Finucane et al. 2015).

Cell‑type specific enrichment:
For AD, SNPs in CNS‑relevant functional elements showed significant enrichment (proportion of h² much larger than proportion of SNPs), indicating that CNS biology is important for AD.

9. Important Considerations
LD score regression requires sufficient sample size (≥5000 individuals) to be reliable.

LD reference population must match the ancestry of the GWAS sample.

Imputation quality correlates with LD score; lower imputation quality can bias results.

For dichotomous traits, the liability‑scale transformation depends on an assumed population prevalence (which may be uncertain).

Intercept in LDSC can detect population stratification: if >1, it suggests inflation due to stratification or cryptic relatedness.

10. Potential Exam Questions (Based on First Exam)
Below are example questions in the style of the first exam. These test understanding of concepts, interpretation, and application.

1) Circle all correct statements about estimating heritability from GWAS summary statistics using LD score regression:
a) LD score regression can estimate heritability without individual-level data.
b) The intercept of the LD score regression reflects the average LD score across SNPs.
c) A positive intercept significantly greater than 1 suggests population stratification.
d) LD score regression requires that all SNPs be independent.

Correct: a, c

2) A GWAS of body mass index (BMI) in 300,000 Europeans yields LD score regression results:
Intercept = 1.02 (SE = 0.01)

Slope = 0.18 (SE = 0.02)

Number of SNPs used = 1,000,000

Sample size = 300,000

a) Calculate the SNP‑based heritability of BMI.
b) What does the intercept tell you about the study?

Answers:
a) h² = (M/N) × slope = (1,000,000 / 300,000) × 0.18 ≈ 0.6.
b) Intercept ≈ 1, so little to no inflation due to population stratification.

3) A researcher performs LD score regression for a binary trait with sample case proportion = 0.05 and assumed population prevalence = 0.01. The observed‑scale heritability is 0.04. Convert this to liability‑scale heritability. Show the formula and explain why the conversion is necessary.
Answer:
Formula: h²_liab = h²_obs × [K(1–K)]² / [p(1–p) Z²], where Z = dnorm(qnorm(K, lower.tail=F)).
The conversion is needed because the observed 0/1 scale depends on the sample's case proportion, whereas the liability scale is an underlying continuous trait that is more comparable across studies.

4) A genetic correlation analysis between two traits gives ρ_g = 0.75 (p = 2e-10). Interpret this result.
Answer:
The traits share 75% of their additive genetic variance, indicating strong shared genetic architecture. Variants that increase one trait tend to increase the other.

5) The table below shows LD score regression results for several functional annotations for a trait. For each annotation, decide whether it is enriched (enrichment >1) and what that implies.
Annotation	Prop. SNPs	Prop. h²	Enrichment
Coding	0.015	0.025	1.67
Intronic	0.388	0.32	0.82
Regulatory	0.063	0.21	3.33
Which annotation contributes the most heritability per SNP? What biological interpretation might you draw?

Answer:
Enrichment = (Prop. h²) / (Prop. SNPs).
Regulatory regions show the highest enrichment (3.33), meaning they contribute disproportionately to heritability. This suggests that variation in regulatory elements is important for the trait.

6) True or False (and explain):
a) GCTA and LDSC always give the same heritability estimate for the same trait.
b) LD score regression can only be used for continuous traits.

Answers:
a) False. They use different methods and data sources; GCTA uses individual‑level data and the empirical GRM, while LDSC uses summary statistics and a reference LD panel. They can differ, especially if the LD reference does not perfectly match the study population.
b) False. LDSC can be used for binary traits by converting to liability scale.

7) A student runs LD score regression on a meta‑analysis of 10 European studies and obtains an intercept of 1.12. What are two possible causes? How could the researcher address this?
Answer:
Possible causes: population stratification within or across studies, cryptic relatedness, or sample overlap.
To address, the researcher could compute ancestry‑principal components and include them as covariates, or use a more homogeneous LD reference (e.g., ancestry‑matched). They could also use the intercept itself to adjust the heritability estimate (the intercept captures inflation, and LDSC automatically accounts for it in the slope, but a high intercept may indicate poor quality control).

11. Summary of Key Commands (from class)
bash
# 1. Munge summary statistics
munge_sumstats.py --sumstats file.gz --snp MarkerName \
  --N-cas 21982 --N-con 41944 --a1 Effect_allele --a2 Non_Effect_allele \
  --signed-sumstats Beta,0 --merge-alleles w_hm3.snplist --out AD

# 2. LD score regression for heritability
ldsc.py --h2 AD.sumstats.gz --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ --pop-prev 0.10 --samp-prev 0.344 --out AD_h2

# 3. Genetic correlation
ldsc.py --rg AD.sumstats.gz,BV.HGT.sumstats.gz --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ --out AD_BV_rg

# 4. Partitioned heritability (baseline annotations)
ldsc.py --h2 AD.sumstats.gz --ref-ld-chr baseline. --w-ld-chr weights.hm3_noMHC. \
  --overlap-annot --frqfile-chr 1000G.EUR.QC. --pop-prev 0.10 --samp-prev 0.344 \
  --print-coefficients --out AD_part

# 5. Cell-type enrichment (e.g., CNS)
ldsc.py --h2 AD.sumstats.gz --ref-ld-chr cell_type_group.3.,baseline. \
  --w-ld-chr weights.hm3_noMHC. --overlap-annot --frqfile-chr 1000G.EUR.QC. \
  --pop-prev 0.10 --samp-prev 0.344 --out AD_CNS --print-coefficients
With these notes, you should be able to:

Understand the concepts of heritability and genetic correlation

Explain how GCTA and LDSC work

Interpret LDSC output (intercept, slope, h², enrichment)

Recognize when to use liability‑scale transformation

Answer exam‑style questions about these topics.