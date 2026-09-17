# ANM_MDD_Anxiety_MitoBrain

This repository contains analysis pipelines developed for investigating the
genetic architecture of **age at natural menopause (ANM)** and its relationship
with **major depressive disorder (MDD)**, **anxiety**, **brain cell types**,
**gene expression**, and **mitochondrial biology**.

## Research Overview

The project integrates genome-wide association studies (GWAS), statistical
genetics, fine-mapping, functional genomics, Mendelian randomization, Bayesian
colocalization, and single-cell analyses to identify and prioritize genetic
variants, genes, and regulatory mechanisms associated with reproductive aging
and neuropsychiatric traits.

## Current Analysis Pipelines

### ABF Fine-Mapping of Age at Natural Menopause

The script `ANM_ABF_finemapping.R` performs Approximate Bayes Factor (ABF)
fine-mapping using the **coloc** R package.

Main features include:

- Reading ANM GWAS summary statistics
- Using FUMA-defined genomic risk-locus boundaries for locus-specific fine-mapping
- Running Bayesian ABF fine-mapping for each locus
- Computing posterior inclusion probabilities (PIP) for individual SNPs
- Constructing 90% credible sets
- Identifying the highest-PIP representative variant within each fine-mapped locus
- Generating locus-level summaries
- Exporting publication-ready CSV and Excel output files

### Fine-Mapping Software

- R (≥ 4.2)
- data.table
- coloc
- openxlsx
## FUSION Transcriptome-Wide Association Analysis

Transcriptome-wide association analysis (TWAS) was performed using the
**FUSION** framework to identify genes whose genetically predicted expression
was associated with age at natural menopause (ANM).

Analyses were performed separately in three GTEx v8 tissues selected for their
relevance to reproductive and neuroendocrine biology:

- Ovary
- Pituitary
- Brain - Hypothalamus

Precomputed **GTEx v8 European-ancestry (EUR) gene-expression prediction
models** distributed through FUSION were used. Linkage disequilibrium (LD) was
estimated using the **1000 Genomes Project European reference panel** supplied
for FUSION analyses.

### Input Data

The marginal TWAS analysis requires:

- `ANM_FUSION.sumstats`  
  ANM GWAS summary statistics formatted for FUSION with the required columns
  `SNP`, `A1`, `A2`, and `Z`, where the GWAS Z-score was calculated as
  `beta / SE`.

- GTEx v8 EUR FUSION prediction models  
  Precomputed gene-expression prediction weights for ovary, pituitary, and
  brain hypothalamus.

- `LDREF/1000G.EUR.*`  
  Chromosome-specific 1000 Genomes European LD reference files used by FUSION.

The ANM GWAS used for this analysis was in **GRCh37/hg19** coordinates.
Precomputed GTEx v8 EUR prediction models distributed through FUSION were used
directly; no user-performed GTEx GRCh38-to-GRCh37 liftOver step was included in
the TWAS pipeline.

### Marginal TWAS

Marginal TWAS was performed using `FUSION.assoc_test.R`.

A representative command for ovary is:

```bash
Rscript FUSION/FUSION.assoc_test.R \
  --sumstats ANM_FUSION.sumstats \
  --weights WEIGHTS/Ovary/GTExv8.EUR.Ovary.pos \
  --weights_dir WEIGHTS/Ovary/ \
  --ref_ld_chr LDREF/1000G.EUR. \
  --chr 1 \
  --out RESULTS/Ovary/ANM_Ovary_chr1.dat
```

The analysis was repeated for chromosomes **1–22** and separately for
**ovary, pituitary, and brain hypothalamus**, using the corresponding GTEx v8
EUR prediction-model files.

For each tested gene-expression prediction model, FUSION reports a marginal
TWAS association statistic, including the TWAS Z-score (`TWAS.Z`) and
TWAS P-value (`TWAS.P`).

### Transcriptome-Wide Significance

Transcriptome-wide significance was determined separately within each tissue
using a Bonferroni correction:

$$
P_{\mathrm{TWAS}} < \frac{0.05}{N_{\mathrm{models}}}
$$

where \(N_{\mathrm{models}}\) is the total number of gene-expression prediction
models tested in the corresponding tissue.

Only TWAS associations passing the tissue-specific Bonferroni threshold were
carried forward to joint/conditional analysis.

### Joint/Conditional TWAS

Significant marginal TWAS associations were subsequently evaluated using
FUSION joint/conditional analysis with `FUSION.post_process.R`.

This analysis accounts for correlation among genetically predicted expression
signals within genomic loci and evaluates which TWAS associations remain
supported after conditioning on correlated prediction models.

A representative command is:

```bash
Rscript FUSION/FUSION.post_process.R \
  --sumstats ANM_FUSION.sumstats \
  --input RESULTS/JOINT/Ovary/ANM_Ovary_chr17_significant.dat \
  --out RESULTS/JOINT/Ovary/ANM_Ovary_chr17_joint \
  --ref_ld_chr LDREF/1000G.EUR. \
  --chr 17 \
  --plot \
  --plot_legend joint \
  --locus_win 500000
```

The procedure was repeated for chromosomes containing transcriptome-wide
significant TWAS associations and separately for each of the three tissues.

A **500-kb flanking locus buffer** was used for joint/conditional
post-processing and regional visualization.

### Joint/Conditional TWAS Outputs

FUSION produces two principal joint/conditional result files:

- `*.joint_included.dat`  
  TWAS predictors retained in the joint model. These files contain the marginal
  TWAS statistics together with the corresponding joint estimates and
  P-values.

- `*.joint_dropped.dat`  
  Marginally significant TWAS predictors that were not retained in the joint
  model after accounting for correlated predicted-expression signals. These
  files contain their corresponding conditional statistics.

Thus, the analysis distinguishes between:

1. **Marginal TWAS associations** — gene-expression prediction models
   individually associated with ANM before conditioning; and

2. **Joint/conditional TWAS associations** — results obtained after accounting
   for correlations among genetically predicted expression signals within
   loci.

Joint/conditional results were used to refine the set of ANM-associated
gene-expression signals while recognizing that TWAS association alone does
not establish a causal effect of gene expression on ANM.

### FUSION Software and Resources

The TWAS analysis used:

- R
- FUSION `FUSION.assoc_test.R`
- FUSION `FUSION.post_process.R`
- GTEx v8 EUR precomputed gene-expression prediction models
- 1000 Genomes Project European LD reference data

FUSION software, documentation, GTEx v8 prediction models, and reference
resources are available at:

- http://gusevlab.org/projects/fusion/

For the FUSION TWAS methodology, see:

Gusev A, et al. *Integrative approaches for large-scale transcriptome-wide
association studies.* Nature Genetics. 2016.


## GTEx v8 eQTL-Based Mendelian Randomization

Gene-expression Mendelian randomization (MR) was performed to test whether
genetically predicted expression of ANM-prioritized genes was associated with
age at natural menopause.

Analyses were restricted to **1,294 genes** prioritized from the ANM GWAS using
FUMA.

Conditionally independent cis-eQTLs were obtained from GTEx v8 for 15 tissues:

- 13 brain tissues: Amygdala, Anterior cingulate cortex BA24, Caudate basal
  ganglia, Cerebellar Hemisphere, Cerebellum, Cortex, Frontal Cortex BA9,
  Hippocampus, Hypothalamus, Nucleus accumbens basal ganglia, Putamen basal
  ganglia, Spinal cord cervical c-1, and Substantia nigra
- Ovary
- Pituitary

Because GTEx v8 uses GRCh38 whereas the ANM GWAS used for the MR workflow is in
GRCh37/hg19, GTEx variant coordinates were converted from GRCh38 to GRCh37
using UCSC liftOver before matching and allele harmonization with the ANM GWAS.


### Input Data

The primary MR analysis requires:

- `ANM_GTEx_v8_MR_final_instruments.txt`  
  Harmonized GTEx v8 cis-eQTL–ANM GWAS instrument-level data.

- `ANM_gene_prioritization_FINAL_INTEGRATED_CLEAN.xlsx`  
  FUMA-prioritized ANM genes and corresponding gene symbols.

- `ANM_1000G_EUR_MR_PLINK/`  
  Chromosome-specific PLINK2 reference files generated from 1000 Genomes
  Project Phase 3 European-ancestry samples (GRCh37), used for LD estimation.

- PLINK2  
  Used to calculate signed LD correlation matrices among multiple cis-eQTL
  instruments.


### Primary MR

The primary MR analysis used GTEx v8 conditionally independent cis-eQTL
instruments satisfying:

- `pval_beta < 0.05`
- F-statistic > 10
- MAF ≥ 0.01

Instrument strength was calculated as:

`F = (beta_eQTL / SE_eQTL)^2`

For gene–tissue pairs represented by one usable genetic instrument, MR effects
were estimated using the **Wald ratio**.

For gene–tissue pairs represented by two or more usable instruments, signed
linkage disequilibrium (LD) was estimated from the **1000 Genomes Project
Phase 3 European reference panel** and incorporated into **LD-aware generalized
inverse-variance weighted (IVW) MR**.

The LD correlation matrix was oriented to the GTEx exposure allele before
calculation of the generalized IVW estimate.

When fewer than two instruments could be evaluated using the LD reference or
the LD-based estimator could not be stably computed, the analysis reverted to
a single-instrument Wald ratio.

Benjamini–Hochberg false discovery rate correction was applied across
gene–tissue MR tests, with `MR_FDR < 0.05` considered statistically significant.


### Sensitivity Analysis

A sensitivity MR analysis was performed using the same instrument-quality
criteria as the primary analysis, with the additional requirement:

- GTEx nominal cis-eQTL `P < 5 × 10^-8`

Thus, the sensitivity analysis retained the primary conditionally independent
cis-eQTL, F-statistic, and MAF criteria while additionally requiring the
nominal SNP–gene expression association to reach genome-wide significance.

This stricter analysis assessed whether the MR findings remained significant
when restricted to more strongly supported cis-eQTL instruments.

Genes significant at `MR_FDR < 0.05` in both the primary and sensitivity
analyses were retained as the **robust MR-supported gene set**.


### Main Results

The final robust MR analysis identified **358 genes** that were significant at
`MR_FDR < 0.05` in both the primary and sensitivity analyses.

These 358 genes were retained as the final robust MR-supported gene set for
downstream gene prioritization and colocalization.


### Primary MR Script

#### `ANM_GTEx_v8_MR_FINAL_LDaware.py`

This script performs the primary GTEx v8 gene-expression MR analysis. It:

1. reads harmonized GTEx cis-eQTL instruments;
2. adds gene symbols from the FUMA-prioritized gene table;
3. matches instruments to chromosome-specific 1000 Genomes EUR PLINK reference files;
4. verifies allele correspondence with the LD reference;
5. performs Wald-ratio MR for single-instrument gene–tissue pairs;
6. performs LD-aware generalized IVW MR for multi-instrument gene–tissue pairs;
7. uses a single-instrument Wald-ratio fallback when LD-based estimation cannot
   be performed using at least two usable instruments;
8. calculates MR beta, standard error, Z-score, P-value, and 95% confidence interval;
9. applies Benjamini–Hochberg FDR correction; and
10. saves the complete MR results and QC tables.


### Sensitivity MR Script

#### `ANM_GTEx_v8_MR_SENSITIVITY_P5e8.py`

This script repeats the MR analysis while retaining the primary
instrument-quality criteria and additionally restricting instruments to GTEx
cis-eQTLs with `pval_nominal < 5 × 10^-8`.

The resulting significant genes are compared with those from the primary MR
analysis to identify MR-supported genes robust to stricter instrument selection.


### MR Software Requirements

- Python 3
- pandas
- NumPy
- SciPy
- openpyxl
- PLINK2
- UCSC liftOver
- 1000 Genomes Project Phase 3 European-ancestry reference data


## GTEx v8 Colocalization of MR-Prioritized Genes

Bayesian colocalization was performed to evaluate whether ANM GWAS associations
and GTEx v8 cis-eQTL associations for robust MR-supported genes were consistent
with a shared causal variant.

For each eligible gene, the GTEx v8 tissue corresponding to
`MR_Best_shared_tissue` was selected. A symmetric **±500 kb** region centered
on the representative ANM fine-mapped variant was defined in **GRCh38**
coordinates.

Regional GTEx v8 cis-eQTL summary statistics were accessed through the
tabix-indexed imported datasets distributed by the **eQTL Catalogue**. Unlike
the MR analysis, which uses selected conditionally independent cis-eQTL
instruments, the colocalization analysis uses all matched and harmonized
regional cis-eQTL and ANM GWAS variants that pass variant-level quality control.

The eQTL Catalogue GTEx v8 tabix file locations are documented here:

- `https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/tabix_ftp_paths_imported.tsv`


### Colocalization Workflow

The colocalization workflow consists of two scripts:

#### `04_prepare_GTEx_coloc_inputs.R`

This script:

1. reads the robust MR-supported gene manifest;
2. selects the MR-best shared tissue for each eligible gene;
3. locates the corresponding ANM representative fine-mapped variant in GRCh38;
4. defines a symmetric ±500 kb region;
5. queries the corresponding GTEx v8 regional cis-eQTL summary statistics from
   the eQTL Catalogue tabix resources;
6. retains rows corresponding to the target Ensembl gene; and
7. saves one local gene-by-tissue regional file for downstream colocalization.


#### `05_MR_prioritized_gene_colocalization.R`

This script:

1. reads the prepared regional GTEx v8 cis-eQTL files;
2. extracts the corresponding ANM GWAS variants;
3. matches variants by genomic position and rsID;
4. harmonizes alleles between GTEx and ANM;
5. removes allele mismatches and invalid MAF values;
6. estimates phenotype standard deviations required by `coloc`;
7. runs `coloc.abf`;
8. calculates regional posterior probabilities PP.H0–PP.H4;
9. calculates SNP-level `SNP.PP.H4` values;
10. identifies the top candidate shared variant based on `SNP.PP.H4`;
11. performs sensitivity analyses across alternative `p12` priors; and
12. exports SNP-level and regional colocalization summaries.


### Colocalization Parameters

Primary prior probabilities were:

- `p1 = 1 × 10^-4`
- `p2 = 1 × 10^-4`
- `p12 = 1 × 10^-5`

Sensitivity analyses evaluated:

- `p12 = 1 × 10^-6`
- `p12 = 5 × 10^-6`
- `p12 = 1 × 10^-5`
- `p12 = 5 × 10^-5`
- `p12 = 1 × 10^-4`

Regional colocalization was interpreted as:

- **Strong H4 support:** `PP.H4 ≥ 0.80`
- **Moderate H4 support:** `0.50 ≤ PP.H4 < 0.80`
- **H3 favored / distinct signals:** `PP.H3 > PP.H4`

`SNP.PP.H4` is interpreted separately from regional `PP.H4`: it represents the
posterior probability assigned to an individual candidate shared variant
conditional on H4, whereas regional `PP.H4` quantifies support for the
shared-causal-variant hypothesis itself.


### Colocalization Results

Among **358 robust MR-supported genes**, **352** had at least one gene–tissue
association significant at `MR_FDR < 0.05` in both the primary and sensitivity
MR analyses and were therefore eligible for the predefined best-tissue
colocalization workflow.

Of these:

- **344 genes** had usable GTEx v8 regional cis-eQTL data and completed colocalization;
- **8 genes** lacked target-gene GTEx cis-eQTL rows within the predefined
  best-tissue ±500 kb window;
- **44 genes** showed strong colocalization support with `PP.H4 ≥ 0.80`;
- **13 genes** showed moderate H4 support with `0.50 ≤ PP.H4 < 0.80`; and
- **278 genes** were H3-favored (`PP.H3 > PP.H4`), indicating greater posterior
  support for distinct rather than shared causal signals within the tested region.

Strong colocalization was incorporated as an additional gene-prioritization
evidence layer. The original `Evidence_count` was retained, and
`Evidence_count_with_coloc` adds one point for genes with
`Coloc_strong_H4 = TRUE`.


### Colocalization Software Requirements

- R
- data.table
- openxlsx
- coloc
- seqminer


## Data Sources

- **ANM GWAS:** GWAS Catalog study accession `GCST90292553`
- **GTEx v8:** GTEx Portal
- **GTEx v8 regional cis-eQTL summary statistics for colocalization:** eQTL Catalogue imported GTEx v8 tabix resources
- **1000 Genomes Project Phase 3 European reference panel:** International Genome Sample Resource
- **MDD and anxiety GWAS:** Psychiatric Genomics Consortium


## Citation

If you use this code, please cite the associated publication when available.

---

**Author**

Suzan Farhang-Sardroodi
