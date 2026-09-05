# ============================================================
# Age at Natural Menopause (ANM) Colocalization Analysis
# ============================================================
#
# Author: Suzan Farhang-Sardroodi
#
# This code was developed by Suzan Farhang-Sardroodi for the
# colocalization analyses of age at natural menopause (ANM)
# reported in this study.
#
# The analysis evaluates whether ANM GWAS and GTEx v8 cis-eQTL
# association signals for MR-prioritized genes are consistent
# with a shared regional causal variant.
#
# Please cite the associated publication when using or adapting
# this code.
# =============================================================================
# 05_MR_prioritized_gene_colocalization.R
#
# Purpose
# -------
# Systematically test whether ANM GWAS associations and GTEx v8 cis-eQTL
# associations for MR-prioritized genes are consistent with a shared
# regional causal variant.
#
# Analysis design
# ---------------
# For each robust MR-supported gene:
#
#   * Tissue:
#       MR_Best_shared_tissue
#
#   * Locus:
#       symmetric +/-500 kb around the independent ANM representative SNP
#
#   * Genome build:
#       GRCh38
#
#   * Colocalization method:
#       coloc.abf
#
#   * Primary priors:
#       p1  = 1e-4
#       p2  = 1e-4
#       p12 = 1e-5
#
#   * Strong colocalization:
#       PP.H4 >= 0.80
#
#   * Moderate H4 support:
#       0.50 <= PP.H4 < 0.80
#
#   * H3-favored:
#       PP.H3 > PP.H4
#
# Interpretation
# --------------
# PP.H4 is the regional posterior probability supporting a shared causal
# variant between the ANM GWAS and eQTL association.
#
# SNP.PP.H4 is a different quantity: it distributes posterior probability
# among individual variants conditional on H4. Therefore, a high SNP.PP.H4
# does not by itself establish regional colocalization when PP.H4 is low.
#
# Author
# ------
# Suzan Farhang-Sardroodi
# =============================================================================


suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(coloc)
})


# =============================================================================
# 1. PROJECT PATHS
# =============================================================================

PROJECT_DIR <- "."


MASTER_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "processed",
  "ANM_ROBUST_MR_GENES_FOR_COLOC.xlsx"
)


MASTER_SHEET <- "COLOC_ready"


PREP_LOG <- file.path(
  PROJECT_DIR,
  "data",
  "processed",
  "GTEx_coloc_inputs",
  "GTEx_coloc_input_preparation_log.xlsx"
)


ANM_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "gwas",
  "ANM_GCST90292553_GRCh38.txt"
)


GTEX_DIR <- file.path(
  PROJECT_DIR,
  "data",
  "processed",
  "GTEx_coloc_inputs"
)


OUT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "colocalization"
)


dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 2. COLOCALIZATION SETTINGS
# =============================================================================

ANM_N <- 173424L


P1 <- 1e-4

P2 <- 1e-4

P12 <- 1e-5


P12_SENSITIVITY <- c(
  1e-6,
  5e-6,
  1e-5,
  5e-5,
  1e-4
)


# =============================================================================
# 3. READ MR MASTER TABLE
# =============================================================================

message("Reading MR-prioritized gene table...")

master <- as.data.table(
  read.xlsx(
    MASTER_FILE,
    sheet = MASTER_SHEET
  )
)


# =============================================================================
# 4. READ GTEx PREPARATION LOG
# =============================================================================

prep <- as.data.table(
  read.xlsx(
    PREP_LOG,
    sheet = "Preparation_log"
  )
)


prep_ok <- prep[
  Status %in% c(
    "OK",
    "EXISTING_OK"
  )
]


message(
  "Prepared gene-by-tissue datasets available: ",
  nrow(prep_ok)
)


# =============================================================================
# 5. READ ANM GWAS ONCE
# =============================================================================

message("Reading ANM GRCh38 GWAS...")

anm <- fread(
  ANM_FILE
)


anm[, chromosome := as.character(chromosome)]

anm[, base_pair_location := as.integer(base_pair_location)]

anm[, effect_allele := toupper(effect_allele)]

anm[, other_allele := toupper(other_allele)]


# =============================================================================
# 6. OUTPUT HOLDERS
# =============================================================================

summary_list <- list()

status_list <- list()


# =============================================================================
# 7. MAIN GENE-BY-TISSUE LOOP
# =============================================================================

for (i in seq_len(nrow(prep_ok))) {


  gene <- as.character(
    prep_ok$Gene[i]
  )


  ensg <- as.character(
    prep_ok$Ensembl_ID[i]
  )


  tissue <- as.character(
    prep_ok$Tissue[i]
  )


  representative_snp <- as.character(
    prep_ok$Representative_SNP[i]
  )


  chr <- as.character(
    prep_ok$Chromosome[i]
  )


  center_pos <- as.integer(
    prep_ok$Center_GRCh38[i]
  )


  region_start <- as.integer(
    prep_ok$Region_start[i]
  )


  region_end <- as.integer(
    prep_ok$Region_end[i]
  )


  message("")
  message("============================================================")
  message("COLOC ", i, " of ", nrow(prep_ok))
  message("Gene: ", gene)
  message("Tissue: ", tissue)
  message("Representative SNP: ", representative_snp)
  message("Region: chr", chr, ":", region_start, "-", region_end)
  message("============================================================")


  # ===========================================================================
  # 7A. MATCH MR RESULT
  # ===========================================================================

  mr_row <- master[
    Gene == gene &
    Ensembl_ID == ensg
  ]


  if (nrow(mr_row) == 0) {

    status_list[[length(status_list) + 1]] <- data.table(
      Gene = gene,
      Ensembl_ID = ensg,
      Tissue = tissue,
      Status = "NO_MASTER_ROW"
    )

    next
  }


  mr_row <- mr_row[1]


  # ===========================================================================
  # 7B. READ REGIONAL GTEx cis-eQTL DATA
  # ===========================================================================

  gtex_file <- file.path(
    GTEX_DIR,
    paste0(
      gene,
      "_",
      tissue,
      "_GTEx_v8_GRCh38_symmetric.tsv"
    )
  )


  if (!file.exists(gtex_file)) {

    status_list[[length(status_list) + 1]] <- data.table(
      Gene = gene,
      Ensembl_ID = ensg,
      Tissue = tissue,
      Status = "GTEX_FILE_MISSING"
    )

    next
  }


  eq <- fread(
    gtex_file
  )


  eq <- eq[
    gene_id == ensg
  ]


  if (nrow(eq) == 0) {

    status_list[[length(status_list) + 1]] <- data.table(
      Gene = gene,
      Ensembl_ID = ensg,
      Tissue = tissue,
      Status = "NO_GTEX_GENE_ROWS"
    )

    next
  }


  eq[, chromosome := as.character(chromosome)]

  eq[, position := as.integer(position)]

  eq[, ref := toupper(ref)]

  eq[, alt := toupper(alt)]


  eq <- eq[
    chromosome == chr &
    position >= region_start &
    position <= region_end
  ]


  N_GTEX_REGION <- nrow(
    eq
  )


  if (N_GTEX_REGION == 0) {
    next
  }


  # ===========================================================================
  # 7C. DERIVE GTEx SAMPLE SIZE
  #
  # eQTL Catalogue provides allele number (AN).
  # For diploid autosomal variants:
  #
  #      N ~= AN / 2
  #
  # The maximum valid allele number in the regional dataset is used.
  # ===========================================================================

  valid_an <- eq$an[
    !is.na(eq$an) &
    is.finite(eq$an) &
    eq$an > 0
  ]


  if (length(valid_an) == 0) {
    next
  }


  GTEX_N <- as.integer(
    round(
      max(valid_an) / 2
    )
  )


  # ===========================================================================
  # 7D. ANM REGIONAL DATA
  # ===========================================================================

  anm_region <- anm[
    chromosome == chr &
    base_pair_location >= region_start &
    base_pair_location <= region_end
  ]


  N_ANM_REGION <- nrow(
    anm_region
  )


  if (N_ANM_REGION == 0) {
    next
  }


  # ===========================================================================
  # 7E. MATCH GTEx AND ANM BY CHROMOSOME + POSITION
  # ===========================================================================

  m <- merge(
    eq,
    anm_region,
    by.x = c(
      "chromosome",
      "position"
    ),
    by.y = c(
      "chromosome",
      "base_pair_location"
    )
  )


  N_POSITION_OVERLAP <- nrow(
    m
  )


  if (N_POSITION_OVERLAP == 0) {
    next
  }


  # ===========================================================================
  # 7F. EXACT rsID MATCH
  # ===========================================================================

  m[, rsid_match :=
      rsid == rs_id]


  N_RSID_MATCH <- sum(
    m$rsid_match,
    na.rm = TRUE
  )


  m <- m[
    rsid_match == TRUE
  ]


  if (nrow(m) == 0) {
    next
  }


  # ===========================================================================
  # 7G. ALLELE HARMONIZATION
  #
  # GTEx beta is reported for ALT relative to REF.
  #
  # If ANM effect allele = GTEx ALT:
  #       retain ANM beta.
  #
  # If ANM effect allele = GTEx REF:
  #       reverse the sign of the ANM beta.
  # ===========================================================================

  m[, allele_status := fifelse(

    effect_allele == alt &
    other_allele == ref,

    "SAME",

    fifelse(

      effect_allele == ref &
      other_allele == alt,

      "REVERSED",

      "MISMATCH"
    )
  )]


  N_SAME <- sum(
    m$allele_status == "SAME"
  )


  N_REVERSED <- sum(
    m$allele_status == "REVERSED"
  )


  N_MISMATCH <- sum(
    m$allele_status == "MISMATCH"
  )


  m <- m[
    allele_status != "MISMATCH"
  ]


  if (nrow(m) == 0) {
    next
  }


  m[, beta_ANM_harmonized :=
      fifelse(
        allele_status == "SAME",
        beta.y,
        -beta.y
      )]


  # ===========================================================================
  # 7H. MINOR-ALLELE FREQUENCY QC
  # ===========================================================================

  m[, ANM_MAF := pmin(
    effect_allele_frequency,
    1 - effect_allele_frequency
  )]


  m[, GTEx_MAF := maf]


  N_BEFORE_MAF <- nrow(
    m
  )


  m <- m[
    !is.na(ANM_MAF) &
    !is.na(GTEx_MAF) &
    ANM_MAF > 0 &
    ANM_MAF < 0.5 &
    GTEx_MAF > 0 &
    GTEx_MAF < 0.5
  ]


  N_REMOVED_MAF <-
    N_BEFORE_MAF -
    nrow(m)


  N_FINAL_SNPS <- nrow(
    m
  )


  if (N_FINAL_SNPS == 0) {
    next
  }


  # ===========================================================================
  # 7I. BUILD coloc DATASETS
  #
  # All matched regional variants surviving QC are supplied to coloc.abf.
  #
  # This differs from the MR analysis, which uses selected genetic
  # instruments. Coloc requires the regional association pattern.
  # ===========================================================================

  D1 <- list(

    beta =
      m$beta_ANM_harmonized,

    varbeta =
      m$standard_error^2,

    N =
      ANM_N,

    MAF =
      m$ANM_MAF,

    snp =
      m$rsid,

    position =
      m$position,

    type =
      "quant"
  )


  D2 <- list(

    beta =
      m$beta.x,

    varbeta =
      m$se^2,

    N =
      GTEX_N,

    MAF =
      m$GTEx_MAF,

    snp =
      m$rsid,

    position =
      m$position,

    type =
      "quant"
  )


  # ===========================================================================
  # 7J. DATASET VALIDATION
  # ===========================================================================

  check_dataset(
    D1
  )


  check_dataset(
    D2
  )


  # ===========================================================================
  # 7K. PHENOTYPE SD ESTIMATION
  #
  # coloc requires sdY for quantitative traits.
  #
  # Where direct sdY values were unavailable, sdY was estimated from:
  #
  #     sampling variance
  #     MAF
  #     sample size
  #
  # using coloc::sdY.est.
  # ===========================================================================

  sdY_ANM <- suppressWarnings(
    coloc:::sdY.est(
      vbeta = D1$varbeta,
      maf = D1$MAF,
      n = D1$N
    )
  )


  sdY_GTEx <- suppressWarnings(
    coloc:::sdY.est(
      vbeta = D2$varbeta,
      maf = D2$MAF,
      n = D2$N
    )
  )


  D1$sdY <- sdY_ANM

  D2$sdY <- sdY_GTEx


  # ===========================================================================
  # 7L. PRIMARY coloc.abf ANALYSIS
  # ===========================================================================

  res <- coloc.abf(
    dataset1 = D1,
    dataset2 = D2,
    p1 = P1,
    p2 = P2,
    p12 = P12
  )


  PP_H0 <- as.numeric(
    res$summary["PP.H0.abf"]
  )

  PP_H1 <- as.numeric(
    res$summary["PP.H1.abf"]
  )

  PP_H2 <- as.numeric(
    res$summary["PP.H2.abf"]
  )

  PP_H3 <- as.numeric(
    res$summary["PP.H3.abf"]
  )

  PP_H4 <- as.numeric(
    res$summary["PP.H4.abf"]
  )


  message(
    "PP.H4 = ",
    signif(PP_H4, 6),
    " | PP.H3 = ",
    signif(PP_H3, 6)
  )


  # ===========================================================================
  # 7M. SNP-LEVEL H4 POSTERIORS
  #
  # SNP.PP.H4 ranks candidate shared variants conditional on H4.
  #
  # It must not be interpreted as the regional probability of
  # colocalization.
  # ===========================================================================

  snp_results <- as.data.table(
    res$results
  )


  setorder(
    snp_results,
    -SNP.PP.H4
  )


  snp_results[, rank := .I]


  top_snp <- snp_results[1]


  # ===========================================================================
  # 7N. SHARED-PRIOR SENSITIVITY ANALYSIS
  # ===========================================================================

  prior_sensitivity <- data.table(
    p12 = P12_SENSITIVITY,
    PP_H3 = NA_real_,
    PP_H4 = NA_real_
  )


  for (
    j in seq_along(
      P12_SENSITIVITY
    )
  ) {

    rr <- coloc.abf(
      dataset1 = D1,
      dataset2 = D2,
      p1 = P1,
      p2 = P2,
      p12 = P12_SENSITIVITY[j]
    )


    prior_sensitivity[
      j,
      PP_H3 :=
        as.numeric(
          rr$summary[
            "PP.H3.abf"
          ]
        )
    ]


    prior_sensitivity[
      j,
      PP_H4 :=
        as.numeric(
          rr$summary[
            "PP.H4.abf"
          ]
        )
    ]
  }


  # ===========================================================================
  # 7O. SNP-LEVEL OUTPUT TABLE
  # ===========================================================================

  overlap_table <- m[
    ,
    .(
      snp = rsid,
      chr = chromosome,
      pos = position,

      gwas_beta = beta_ANM_harmonized,
      gwas_se = standard_error,
      gwas_n = ANM_N,
      gwas_maf = ANM_MAF,

      qtl_beta = beta.x,
      qtl_se = se,
      qtl_n = GTEX_N,
      qtl_maf = GTEx_MAF
    )
  ]


  posterior_table <- snp_results[
    ,
    .(
      snp,
      SNP.PP.H4,
      rank
    )
  ]


  overlap_table <- merge(
    overlap_table,
    posterior_table,
    by = "snp",
    all.x = TRUE,
    sort = FALSE
  )


  setorder(
    overlap_table,
    chr,
    pos
  )


  # ===========================================================================
  # 7P. REGIONAL SUMMARY
  # ===========================================================================

  regional_summary <- data.table(

    Gene =
      gene,

    Ensembl_ID =
      ensg,

    Tissue =
      tissue,

    Genome_build =
      "GRCh38",

    Representative_SNP =
      representative_snp,

    Representative_SNP_position_GRCh38 =
      center_pos,

    Window_kb =
      500,

    Region_start_GRCh38 =
      region_start,

    Region_end_GRCh38 =
      region_end,

    N_ANM_regional_variants =
      N_ANM_REGION,

    N_GTEx_regional_rows =
      N_GTEX_REGION,

    GTEx_N =
      GTEX_N,

    N_position_overlap =
      N_POSITION_OVERLAP,

    N_exact_rsID_match =
      N_RSID_MATCH,

    N_allele_same =
      N_SAME,

    N_allele_reversed =
      N_REVERSED,

    N_allele_mismatch =
      N_MISMATCH,

    N_removed_MAF_QC =
      N_REMOVED_MAF,

    N_coloc_SNPs =
      N_FINAL_SNPS,

    PP_H0 =
      PP_H0,

    PP_H1 =
      PP_H1,

    PP_H2 =
      PP_H2,

    PP_H3 =
      PP_H3,

    PP_H4 =
      PP_H4,

    Top_shared_SNP =
      top_snp$snp,

    Top_shared_SNP_position_GRCh38 =
      top_snp$position,

    Top_SNP_PP_H4 =
      top_snp$SNP.PP.H4,

    Estimated_sdY_ANM =
      sdY_ANM,

    Estimated_sdY_GTEx =
      sdY_GTEx,

    Conservative_p12_1e_6_PP_H4 =
      prior_sensitivity[
        p12 == 1e-6,
        PP_H4
      ],

    MR_Best_shared_tissue =
      mr_row$MR_Best_shared_tissue,

    MR_beta =
      mr_row$MR_beta,

    MR_SE =
      mr_row$MR_SE,

    MR_P =
      mr_row$MR_P,

    MR_FDR =
      mr_row$MR_FDR,

    MR_CI_lower =
      mr_row$MR_CI_lower,

    MR_CI_upper =
      mr_row$MR_CI_upper,

    MR_method =
      mr_row$MR_method,

    MR_N_instruments =
      mr_row$MR_N_instruments,

    MR_instrument_rsIDs =
      mr_row$MR_instrument_rsIDs
  )


  # ===========================================================================
  # 7Q. REGIONAL INTERPRETATION
  # ===========================================================================

  regional_summary[
    ,
    Coloc_interpretation :=
      fifelse(

        PP_H4 >= 0.80,

        "Strong H4 support",

        fifelse(

          PP_H4 >= 0.50,

          "Moderate H4 support",

          fifelse(

            PP_H3 > PP_H4,

            "H3 favored / distinct signals",

            "Limited colocalization evidence"
          )
        )
      )
  ]


  # ===========================================================================
  # 7R. SAVE INDIVIDUAL RESULT WORKBOOK
  # ===========================================================================

  individual_file <- file.path(
    OUT_DIR,
    paste0(
      gene,
      "_",
      tissue,
      "_coloc.xlsx"
    )
  )


  wb <- createWorkbook()


  addWorksheet(
    wb,
    "overlap_snps"
  )


  writeData(
    wb,
    "overlap_snps",
    overlap_table,
    withFilter = TRUE
  )


  freezePane(
    wb,
    "overlap_snps",
    firstRow = TRUE
  )


  addWorksheet(
    wb,
    "regional_summary"
  )


  writeData(
    wb,
    "regional_summary",
    regional_summary,
    withFilter = TRUE
  )


  saveWorkbook(
    wb,
    individual_file,
    overwrite = TRUE
  )


  summary_list[[
    length(summary_list) + 1
  ]] <- regional_summary


  status_list[[
    length(status_list) + 1
  ]] <- data.table(
    Gene = gene,
    Ensembl_ID = ensg,
    Tissue = tissue,
    Status = "OK"
  )


  # Checkpoint.
  fwrite(
    rbindlist(
      summary_list,
      fill = TRUE
    ),
    file.path(
      OUT_DIR,
      "CHECKPOINT_coloc_summary.tsv"
    ),
    sep = "\t"
  )
}


# =============================================================================
# 8. COMBINE RESULTS
# =============================================================================

if (length(summary_list) == 0) {

  stop(
    "No colocalization analyses completed."
  )
}


final_summary <- rbindlist(
  summary_list,
  fill = TRUE
)


setorder(
  final_summary,
  -PP_H4
)


status_table <- rbindlist(
  status_list,
  fill = TRUE
)


# =============================================================================
# 9. SUMMARY COUNTS
# =============================================================================

summary_counts <- data.table(

  Metric = c(
    "Successfully completed",
    "Strong H4 support (PP.H4 >= 0.80)",
    "Moderate H4 support (0.50 <= PP.H4 < 0.80)",
    "H3 favored",
    "Other / limited evidence"
  ),

  N = c(

    nrow(
      final_summary
    ),

    sum(
      final_summary$PP_H4 >= 0.80
    ),

    sum(
      final_summary$PP_H4 >= 0.50 &
      final_summary$PP_H4 < 0.80
    ),

    sum(
      final_summary$PP_H3 >
      final_summary$PP_H4
    ),

    sum(
      final_summary$PP_H4 < 0.50 &
      final_summary$PP_H3 <=
      final_summary$PP_H4
    )
  )
)


# =============================================================================
# 10. SAVE COMBINED RESULTS
# =============================================================================

FINAL_FILE <- file.path(
  OUT_DIR,
  "ANM_MR_prioritized_gene_colocalization_summary.xlsx"
)


wb_final <- createWorkbook()


addWorksheet(
  wb_final,
  "Coloc_results"
)


writeData(
  wb_final,
  "Coloc_results",
  final_summary,
  withFilter = TRUE
)


freezePane(
  wb_final,
  "Coloc_results",
  firstRow = TRUE
)


addWorksheet(
  wb_final,
  "Summary_counts"
)


writeData(
  wb_final,
  "Summary_counts",
  summary_counts
)


addWorksheet(
  wb_final,
  "Batch_status"
)


writeData(
  wb_final,
  "Batch_status",
  status_table,
  withFilter = TRUE
)


saveWorkbook(
  wb_final,
  FINAL_FILE,
  overwrite = TRUE
)


# =============================================================================
# 11. FINAL REPORT
# =============================================================================

message("")
message("============================================================")
message("COLOCALIZATION ANALYSIS COMPLETE")
message("============================================================")

message(
  "Successful analyses: ",
  nrow(final_summary)
)

message(
  "Strong H4 (PP.H4 >= 0.80): ",
  sum(
    final_summary$PP_H4 >= 0.80
  )
)

message(
  "Moderate H4 (0.50 <= PP.H4 < 0.80): ",
  sum(
    final_summary$PP_H4 >= 0.50 &
    final_summary$PP_H4 < 0.80
  )
)

message(
  "H3 favored: ",
  sum(
    final_summary$PP_H3 >
    final_summary$PP_H4
  )
)

message("")
message("Final output:")
message(FINAL_FILE)

message("============================================================")
