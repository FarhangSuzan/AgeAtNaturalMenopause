# ANM_MDD_Anxiety_MitoBrain

This repository contains analysis pipelines developed for investigating the genetic architecture of **age at natural menopause (ANM)** and its relationship with **major depressive disorder (MDD)**, **anxiety**, **brain cell types**, and **mitochondrial biology**.

## Research Overview

The project integrates genome-wide association studies (GWAS), statistical genetics, fine-mapping, functional genomics, and single-cell analyses to identify candidate causal variants and genes underlying reproductive aging and neuropsychiatric disorders.

## Current Analysis Pipelines

### ABF Fine-Mapping of Age at Natural Menopause

The script `ANM_ABF_finemapping.R` performs Approximate Bayes Factor (ABF) fine-mapping using the **coloc** R package.

Main features include:

- Reading GWAS summary statistics
- Defining ±500 kb regions around lead GWAS variants
- Running Bayesian ABF fine-mapping for each locus
- Computing posterior probabilities (PIP) for every SNP
- Constructing 90% credible sets
- Generating locus-level summaries
- Exporting publication-ready CSV and Excel output files

## Software

- R (≥ 4.2)
- data.table
- coloc
- openxlsx
## GTEx v8 eQTL-based Mendelian randomization

Gene-expression Mendelian randomization (MR) was performed to test whether genetically predicted expression of ANM-prioritized genes was associated with age at natural menopause (ANM).

### Input data

The analysis requires:

- `ANM_GTEx_v8_MR_final_instruments.txt`  
  Harmonized GTEx v8 cis-eQTL–ANM GWAS instrument-level data.

- `ANM_gene_prioritization_FINAL_INTEGRATED_CLEAN.xlsx`  
  FUMA-prioritized ANM genes and gene symbols.

- `ANM_1000G_EUR_MR_PLINK/`  
  Chromosome-specific PLINK2 reference files generated from 1000 Genomes Project Phase 3 European-ancestry samples (GRCh37), used for LD estimation.

- PLINK2  
  Used to calculate signed LD correlation matrices among multiple cis-eQTL instruments.

### Primary MR

The primary MR analysis used GTEx v8 conditionally independent cis-eQTL instruments satisfying:

- `pval_beta < 0.05`
- F-statistic > 10
- MAF >= 0.01

For gene–tissue pairs with one genetic instrument, MR effects were estimated using the Wald ratio.

For gene–tissue pairs with two or more instruments, signed linkage disequilibrium was estimated from the 1000 Genomes Phase 3 European reference panel and incorporated into generalized inverse-variance weighted (IVW) MR.

The LD matrix was oriented to the GTEx exposure allele before calculation of the generalized IVW estimate.

Benjamini–Hochberg false discovery rate correction was applied across all gene–tissue MR tests, with `MR_FDR < 0.05` considered statistically significant.

### Sensitivity analysis

A sensitivity MR analysis was performed using the same pipeline and instrument-quality criteria, with the additional requirement:

- GTEx nominal cis-eQTL `P < 5 × 10^-8`

This analysis tested whether MR findings were robust to restriction to very strongly associated cis-eQTL instruments.

Genes significant at `MR_FDR < 0.05` in both the primary and sensitivity analyses were retained as the robust MR-supported gene set.

### Main results

Primary MR:
- 894 genes tested
- 4,175 gene–tissue MR tests
- 2,477 significant gene–tissue associations
- 641 significant genes

Sensitivity MR:
- 472 genes tested
- 2,226 gene–tissue MR tests
- 1,477 significant gene–tissue associations
- 359 significant genes

Robust overlap:
- 358 genes significant in both analyses
- 352 genes with at least one identical significant tissue in both analyses
- ### `ANM_GTEx_v8_MR_FINAL_LDaware.py`

This script performs the primary GTEx v8 gene-expression MR analysis. It:

1. reads harmonized GTEx cis-eQTL instruments;
2. adds gene symbols from the FUMA-prioritized gene table;
3. matches instruments to chromosome-specific 1000 Genomes EUR PLINK reference files;
4. verifies allele correspondence with the LD reference;
5. performs Wald-ratio MR for single-instrument gene–tissue pairs;
6. performs LD-aware generalized IVW MR for multi-instrument pairs;
7. calculates MR beta, SE, Z-score, P-value, and 95% confidence interval;
8. applies Benjamini–Hochberg FDR correction; and
9. saves the complete MR results and QC tables.

## Citation

If you use this code, please cite the associated publication when available.

---

**Author**

Suzan Farhang-Sardroodi  
University of Toronto
