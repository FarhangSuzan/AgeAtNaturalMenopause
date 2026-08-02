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

## Citation

If you use this code, please cite the associated publication when available.

---

**Author**

Suzan Farhang-Sardroodi  
University of Toronto
