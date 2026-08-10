# ANM_MDD_Anxiety_MitoBrain

This repository contains analysis pipelines developed for investigating the genetic architecture of **age at natural menopause (ANM)** and its relationship with **major depressive disorder (MDD)**, **anxiety**, **brain cell types**, and **mitochondrial biology**.

## Research Overview

The project integrates genome-wide association studies (GWAS), statistical genetics, fine-mapping, functional genomics, Mendelian randomization, and single-cell analyses to identify and prioritize genetic variants and genes underlying reproductive aging and neuropsychiatric biology.

## Current Analysis Pipelines

### ABF Fine-Mapping of Age at Natural Menopause

The script `ANM_ABF_finemapping.R` performs Approximate Bayes Factor (ABF) fine-mapping using the **coloc** R package.

Main features include:

- Reading GWAS summary statistics
- Defining ±500 kb regions around lead GWAS variants
- Running Bayesian ABF fine-mapping for each locus
- Computing posterior inclusion probabilities (PIP) for individual SNPs
- Constructing 90% credible sets
- Generating locus-level summaries
- Exporting publication-ready CSV and Excel output files

### Fine-Mapping Software

- R (≥ 4.2)
- data.table
- coloc
- openxlsx


## GTEx v8 eQTL-Based Mendelian Randomization

Gene-expression Mendelian randomization (MR) was performed to test whether genetically predicted expression of ANM-prioritized genes was associated with age at natural menopause.

Analyses were restricted to 1,294 genes prioritized from the ANM GWAS using FUMA.

Conditionally independent cis-eQTLs were obtained from GTEx v8 for 15 tissues:

- 13 brain tissues: Amygdala, Anterior cingulate cortex BA24, Caudate basal ganglia, Cerebellar Hemisphere, Cerebellum, Cortex, Frontal Cortex BA9, Hippocampus, Hypothalamus, Nucleus accumbens basal ganglia, Putamen basal ganglia, Spinal cord cervical c-1, and Substantia nigra
- Ovary
- Pituitary

Because GTEx v8 uses GRCh38 whereas the ANM GWAS uses GRCh37/hg19, GTEx variant coordinates were converted from GRCh38 to GRCh37 using UCSC liftOver before matching and allele harmonization with the ANM GWAS.


### Input Data

The primary MR analysis requires:

- `ANM_GTEx_v8_MR_final_instruments.txt`  
  Harmonized GTEx v8 cis-eQTL–ANM GWAS instrument-level data.

- `ANM_gene_prioritization_FINAL_INTEGRATED_CLEAN.xlsx`  
  FUMA-prioritized ANM genes and corresponding gene symbols.

- `ANM_1000G_EUR_MR_PLINK/`  
  Chromosome-specific PLINK2 reference files generated from 1000 Genomes Project Phase 3 European-ancestry samples (GRCh37), used for LD estimation.

- PLINK2  
  Used to calculate signed LD correlation matrices among multiple cis-eQTL instruments.


### Primary MR

The primary MR analysis used GTEx v8 conditionally independent cis-eQTL instruments satisfying:

- `pval_beta < 0.05`
- F-statistic > 10
- MAF ≥ 0.01

For gene–tissue pairs represented by one genetic instrument, MR effects were estimated using the **Wald ratio**.

For gene–tissue pairs represented by two or more instruments, signed linkage disequilibrium (LD) was estimated from the **1000 Genomes Project Phase 3 European reference panel** and incorporated into **LD-aware generalized inverse-variance weighted (IVW) MR**.

The LD correlation matrix was oriented to the GTEx exposure allele before calculation of the generalized IVW estimate.

Benjamini–Hochberg false discovery rate correction was applied across gene–tissue MR tests, with `MR_FDR < 0.05` considered statistically significant.


### Sensitivity Analysis

A sensitivity MR analysis was performed using the same instrument-quality criteria as the primary analysis, with the additional requirement:

- GTEx nominal cis-eQTL `P < 5 × 10^-8`

This stricter analysis assessed whether the MR findings were robust to restriction to cis-eQTL instruments showing very strong associations with gene expression.

Genes significant at `MR_FDR < 0.05` in both the primary and sensitivity analyses were retained as the **robust MR-supported gene set**.


### Main Results

The final robust MR analysis identified **358 genes** that were significant at `MR_FDR < 0.05` in both the primary and sensitivity analyses.

These 358 genes were retained as the final robust MR-supported gene set for downstream gene prioritization.


### Primary MR Script

#### `ANM_GTEx_v8_MR_FINAL_LDaware.py`

This script performs the primary GTEx v8 gene-expression MR analysis. It:

1. reads harmonized GTEx cis-eQTL instruments;
2. adds gene symbols from the FUMA-prioritized gene table;
3. matches instruments to chromosome-specific 1000 Genomes EUR PLINK reference files;
4. verifies allele correspondence with the LD reference;
5. performs Wald-ratio MR for single-instrument gene–tissue pairs;
6. performs LD-aware generalized IVW MR for multi-instrument gene–tissue pairs;
7. calculates MR beta, standard error, Z-score, P-value, and 95% confidence interval;
8. applies Benjamini–Hochberg FDR correction; and
9. saves the complete MR results and QC tables.


### Sensitivity MR Script

#### `ANM_GTEx_v8_MR_SENSITIVITY_P5e8.py`

This script repeats the MR analysis after additionally restricting instruments to GTEx cis-eQTLs with `pval_nominal < 5 × 10^-8`.

The resulting significant genes are compared with those from the primary MR analysis to identify MR-supported genes robust to stricter instrument selection.


### MR Software Requirements

- Python 3
- pandas
- NumPy
- SciPy
- openpyxl
- PLINK2
- UCSC liftOver
- 1000 Genomes Project Phase 3 European-ancestry reference data


## Citation

If you use this code, please cite the associated publication when available.

---

**Author**

Suzan Farhang-Sardroodi  
University of Toronto
