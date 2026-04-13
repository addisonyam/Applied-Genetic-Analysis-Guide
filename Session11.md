## Week 11: Mendelian Randomization
### BS859 Applied Genetic Analysis
### April 8, 2026

These notes are designed to help you understand the concepts, assumptions, and practical implementation of Mendelian Randomization, especially using GCTA‑GSMR.
We will cover:
1. Why MR? – Problems with observational studies.
2. Core principles – Instrumental variables, key assumptions.
3. How MR works – Single variant, multiple variants, inverse‑variance weighted (IVW).
4. Pleiotropy and how to handle it – Median, MR‑Egger, HEIDI‑outlier, multivariable MR.
5. Two‑sample MR – Advantages, requirements, biases.
6. GSMR method in GCTA – Step‑by‑step workflow, file formats, commands.
7. Interpretation of results – Log odds ratios, ORs, plots, reverse MR.
8. Practical tips for the homework – Reformatting, running GSMR, adjusting parameters.

### 1. Why Mendelian Randomization?
##### The problem with observational studies
- Observational studies can show associations (e.g., high LDL → heart disease), but association does not imply causation.
- Confounding – a third variable (e.g., smoking, diet) affects both exposure and outcome.
- Reverse causation – the outcome may influence the exposure (e.g., early disease changes behaviour).
#### The gold standard: Randomised Controlled Trial (RCT)
- Random assignment eliminates confounding.
- But RCTs are expensive, often unethical, or impractical (e.g., lifelong exposure).
#### The MR solution
- Use genetic variants as instrumental variables (IVs).
- Genotypes are randomly assigned at conception → naturally randomised.
- If a genetic variant is associated with an exposure (e.g., LDL cholesterol), and the exposure truly causes the outcome (e.g., heart disease), then the variant should be associated with the outcome.
- MR mimics an RCT, but with lifelong exposure and without needing to intervene.

### 2. Instrumental Variables (IVs) and the Three Core Assumptions
In MR, a genetic variant Z is used as an instrument for exposure X to infer its causal effect on outcome Y.

Diagram (from slide 14)
```text
     U (confounders)
       ↓      ↓
Z  →  X  →  Y
```
- Z = genetic variant (e.g., SNP).
- X = exposure (e.g., BMI, HDL cholesterol).
- Y = outcome (e.g., type 2 diabetes).
- U = unmeasured confounders.

| Assumption            | Meaning                                                       | How to check (roughly)                                         |
|-----------------------|---------------------------------------------------------------|----------------------------------------------------------------|
| Relevance             | Z is associated with X.                                       | F-statistic > 10, or p < 5×10⁻⁸ in GWAS.                       |
| Independence          | Z is independent of confounders U.                            | Check for population stratification; use principal components. |
| Exclusion restriction | Z affects Y only through X (no direct effect, no pleiotropy). | Biological knowledge; sensitivity analyses (MR‑Egger, HEIDI).  |

If these hold, then `β_XY = β_ZY / β_ZX` gives an unbiased causal estimate.

### 3. Basic MR Estimation
##### Single variant (slide 15–16)
- Regress X on Z: get β_ZX (effect of SNP on exposure).
- Regress Y on Z: get β_ZY (effect of SNP on outcome).
- Causal estimate: β_XY = β_ZY / β_ZX.
- Test for causality: H0: β_ZY = 0 (i.e., no association between Z and Y).
##### Multiple variants (slide 33–34)
- Using many SNPs increases power.
- Inverse‑variance weighted (IVW) – most common:
    - For each SNP i, compute β_XYi = β_ZYi / β_ZXi.
    - Combine as weighted average: β_XY = Σ (β_XYi * w_i) / Σ w_i, where w_i = 1 / var(β_XYi).
    - Equivalent to a linear regression through the origin of β_ZY on β_ZX (slide 34).
##### Two‑sample MR (slide 30)
- Use one GWAS for Z→X and a different GWAS for Z→Y.
- Samples must be from the same population.
- Overlap between samples can bias towards the observational association.
- Two‑sample MR is powerful because you can use large, independent GWAS.

### 4. Pleiotropy – The Main Challenge
##### What is pleiotropy?
- A genetic variant affects multiple traits.
- Horizontal pleiotropy (bad for MR): Z affects Y through a pathway other than X.
    - Violates exclusion restriction.
    - Biases causal estimate.
- Vertical pleiotropy (okay): Z affects X, which affects a mediator, which affects Y – still only through X.

##### How to detect and correct for pleiotropy
| Method           | How it works                                                           | When it works                                                                                      |
|------------------|------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------|
| Weighted median  | Take the median of all β_XYi (weighted by precision).                  | If ≥50% of the weight comes from valid instruments.                                                |
| MR‑Egger         | Regression with an intercept: β_ZYi = β_0 + β_XY * β_ZXi + ε_i.        | If pleiotropic effects are independent of β_ZX (INstrument Strength Independent of Direct Effect). |
| HEIDI‑outlier    | Compare each SNP’s ratio estimate to a target SNP; remove outliers.    | Removes SNPs with large pleiotropic effects.                                                       |
| Multivariable MR | Include multiple exposures (e.g., LDL and triglycerides) in the model. | If all pathways are included, pleiotropy disappears.                                               |

In GSMR (GCTA), the default is to use HEIDI‑outlier removal (p < 0.01) plus LD clumping.

### 5. GSMR – Generalized Summary Mendelian Randomization
GSMR is an implementation in GCTA that:
- Uses multiple SNPs as instruments.
- Accounts for LD among SNPs (by pruning).
- Uses a generalised least squares (GLS) approach that accounts for sampling variance in both β_ZX and β_ZY (more accurate than first‑order delta).
- Performs HEIDI‑outlier test to remove pleiotropic SNPs.

##### Input requirements (slide 70–71)
1. Exposure GWAS – summary statistics with columns: SNP, A1, A2, freq, b, se, p, N.
2. Outcome GWAS – same format, non‑overlapping sample (preferably).
3. LD reference panel – individual‑level genotypes (e.g., 1000 Genomes) to estimate LD for pruning.

##### File format (slide 82)
```text
SNP  A1  A2  freq  b  se  p  N
rs123  A  G  0.25  0.05  0.01  1e-8  50000
...
```
- A1 = effect allele (the one whose effect is given by b).
- freq = frequency of A1.

##### Running GSMR (slide 87)
```bash
gcta64 --gsmr-file exposure.txt outcome.txt \
       --bfile /path/to/ld_reference \
       --gsmr-direction 0 \
       --out my_output \
       --effect-plot
```
- exposure.txt and outcome.txt are simple text files listing names and file paths.
Example exposure.txt:

```text
BMI bmi.ss.txt
```

Example outcome.txt:
```text
T2D T2D.ss.txt
T2De T2De.ss.txt
```
- --gsmr-direction 0 : exposure → outcome.
- --gsmr-direction 1 : outcome → exposure (reverse).
- --gsmr-direction 2 : both directions.
- --effect-plot produces a file *.eff_plot.gz for plotting.

##### Output files (slide 88)
- .log – always check first for warnings (allele frequency differences, MAF < 0.01, etc.)
- .gsmr – the main results table: Exposure, Outcome, bxy (causal effect), se, p, number of SNPs.
- .pleio_snps – list of SNPs removed by HEIDI‑outlier.
- .freq.badsnps – SNPs with large allele frequency mismatch between GWAS and reference.

##### Plotting with R (slides 96–97)
GCTA provides gsmr_plot.r. You run a wrapper script plot.R:

```r
source("gsmr_plot.r")
gsmr_data <- read_gsmr_data("my_output.eff_plot.gz")
jpeg("plot.jpeg", width=960, height=480)
par(mfrow=c(1,2))
plot_gsmr_effect(gsmr_data, "BMI", "T2D", colors()[75])
plot_gsmr_effect(gsmr_data, "BMI", "T2De", colors()[75])
dev.off()
```
The plot shows each SNP’s β_ZX (x‑axis) vs β_ZY (y‑axis).
- The dashed line is the estimated causal slope β_XY.
- If all instruments are valid, points should scatter around that line.
- Outliers (far from line) are potential pleiotropic SNPs.

### 6. Interpreting Results
##### Causal effect bxy
- For a continuous exposure (e.g., BMI in SD units) and a binary outcome (e.g., T2D), bxy is the log odds ratio for a 1‑SD increase in exposure.
- Odds ratio (OR) = exp(bxy).
- Example from slides (BMI → T2D): bxy = 1.097 → OR = exp(1.097) = 3.00.
Interpretation: each 1‑SD increase in BMI multiplies the odds of T2D by 3.
##### Reverse MR
- If you also run --gsmr-direction 1 (outcome → exposure), you might get a significant estimate in both directions.
- This could be:
    - True bidirectional causality (rare).
    - Horizontal pleiotropy (SNPs affect both traits independently).
    - Genetic correlation (shared biology, not causation).
    - Weak instruments in one direction.
- In the BMI‑T2D example, reverse MR gave bxy = -0.0258 (T2D → BMI). That is very small but significant (p tiny).
Likely not a meaningful causal effect; possibly due to pleiotropy or sample overlap.

### 7. Practical Steps for Your Homework (FG & HDL)
You will apply GSMR to two exposures: Fasting Glucose (FG) and HDL cholesterol, with T2D as outcome.

##### Step 1 – Reformat the GWAS summary statistics
For each file, you need columns: SNP A1 A2 freq b se p N.
- FG file: /projectnb/bs859/data/FG/MAGIC1000G_FG_EUR.tsv
Read the README to identify which columns correspond to:
    - SNP ID
    - Effect allele (A1)
    - Other allele (A2)
    - Frequency of A1 (freq)
    - Beta (b)
    - Standard error (se)
    - P‑value (p)
    - Sample size (N) – usually reported in the paper, may be constant.

Use awk to extract and reorder, similar to the BMI example.
- HDL file:` /projectnb/bs859/data/lipids/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results`
Note: "Alt" is the effect allele.
Again, identify columns and use awk to produce GCTA format.
Include explicit NA values if missing (e.g., NA).

##### Step 2 – Remove duplicate SNPs (optional but recommended)
```bash
sort -t ' ' -k 1,1 -u myfile.ss.txt > myfile.ss.dedup.txt
mv myfile.ss.dedup.txt myfile.ss.txt
```
##### Step 3 – Create exposure and outcome list files
For FG → T2D:

```bash
echo "FG FG.ss.txt" > exposure.txt
echo "T2D T2D.ss.txt" > outcome.txt   # use Xue et al. file
```
For HDL → T2D, similar.

##### Step 4 – Run GSMR with default parameters
```bash
gcta64 --gsmr-file exposure.txt outcome.txt \
       --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR \
       --gsmr-direction 0 \
       --out FG_T2D_gsmr \
       --effect-plot
```
##### Step 5 – Check the log file
```bash
less FG_T2D_gsmr.log
```
Look for:
- Number of genome‑wide significant SNPs in common.
- Warnings about allele frequency differences.
- Number of SNPs after LD clumping (default r² < 0.05).
- Number of SNPs removed by HEIDI‑outlier (these are pleiotropic).

##### Step 6 – Extract results
```bash
cat FG_T2D_gsmr.gsmr
```
This gives bxy, se, p, and number of SNPs.

Compute OR = exp(bxy), and 95% CI = exp(bxy ± 1.96*se).

##### Step 7 – Plot
Use plot.R adapted for your exposure and outcome names.
Comment on the plot: do most SNPs lie near the dashed line? Are there outliers?

##### Step 8 – Reverse MR
Change --gsmr-direction 1 (or use 2) and rerun. Interpret.

##### Step 9 – Sensitivity analyses (for HDL)
- Change SNP selection threshold – use --gsmr-snp-min 1e-10 (instead of default 5e‑8) to select only very strong instruments.
- Change HEIDI threshold – use --heidi-thresh 0.05 (instead of default 0.01) to remove more SNPs.

Then compare results: does the causal estimate change? Are outliers removed?

### 8. Common Pitfalls & How to Avoid Them
| Pitfall                            | Consequence                                   | Solution                                                            |
|------------------------------------|-----------------------------------------------|---------------------------------------------------------------------|
| Sample overlap                     | Bias towards observational estimate           | Use independent GWAS for exposure and outcome.                      |
| Population stratification          | Violates independence assumption              | Use GWAS that adjusted for PCs; check allele frequency consistency. |
| Weak instruments (low F‑statistic) | Bias, especially in MR‑Egger                  | Use only genome‑wide significant SNPs (p < 5e‑8).                   |
| Horizontal pleiotropy              | Biased causal estimate                        | Use sensitivity methods (HEIDI, MR‑Egger, weighted median).         |
| LD between instruments             | Overlap in information, can cause overfitting | Prune SNPs (GSMR does this by default: r² < 0.05).                  |
| Different allele coding            | Wrong sign of β_ZX or β_ZY                    | Ensure A1 (effect allele) is consistent across files.               |

### 9. Key Terminology – Quick Reference
- Exposure – The risk factor you think causes the outcome (e.g., BMI, FG, HDL).
- Outcome – The disease or trait of interest (e.g., T2D).
- Instrument – Genetic variant (SNP) used to proxy the exposure.
- Relevance – SNP is associated with exposure.
- Independence – SNP is not confounded with outcome.
- Exclusion restriction – SNP only affects outcome through exposure.
- Pleiotropy – SNP affects multiple traits.
- Horizontal pleiotropy – SNP affects outcome via a pathway separate from the exposure (bad).
- Vertical pleiotropy – SNP affects exposure, which then affects a mediator (still okay).
- HEIDI‑outlier – Test to remove SNPs with heterogeneous ratio estimates.
- IVW – Inverse‑variance weighted estimator (standard fixed‑effect meta‑analysis of ratio estimates).
- GSMR – Generalised Summary Mendelian Randomisation (implements IVW with GLS and HEIDI).
- bxy – Log odds ratio for a 1‑SD increase in exposure (if outcome is binary).
- OR – Odds ratio = exp(bxy).

### 10. Preparing for the Exam
You should be able to:
- Explain why MR is useful and how it overcomes confounding.
- State the three core assumptions and give an example of violation for each.
- Calculate a simple MR estimate given β_ZX and β_ZY.
- Describe the IVW method and when it is consistent.
- Distinguish between horizontal and vertical pleiotropy.
- Name at least three methods to handle pleiotropy and their key conditions (e.g., weighted median requires ≥50% valid instruments).
- Interpret GSMR output – what bxy, se, p, and the plot tell you.
- Explain the purpose of reverse MR and why bidirectional results can be misleading.
- Identify potential biases in two‑sample MR (sample overlap, different populations, weak instruments).