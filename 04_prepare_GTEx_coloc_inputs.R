# =============================================================================
# 04_prepare_GTEx_coloc_inputs.R
#
# Purpose
# -------
# Prepare regional GTEx v8 cis-eQTL summary-statistic files for systematic
# colocalization of MR-prioritized genes with age at natural menopause (ANM).
#
# For each gene:
#   1. Read the MR-best GTEx tissue and ANM representative SNP.
#   2. Locate the representative SNP in the ANM GRCh38 GWAS.
#   3. Define a symmetric +/-500 kb genomic window.
#   4. Query only that interval from the tabix-indexed GTEx v8 resources
#      distributed through the eQTL Catalogue.
#   5. Retain cis-eQTL rows corresponding to the target Ensembl gene.
#   6. Save one small gene-by-tissue regional file for downstream coloc.abf.
#
# Data source
# -----------
# eQTL Catalogue imported GTEx v8 summary statistics:
#
# https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/
# blob/master/tabix/tabix_ftp_paths_imported.tsv
#
# Genome build
# ------------
# GRCh38
#
# Author
# ------
# Suzan Farhang-Sardroodi
# =============================================================================


suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(seqminer)
})


# =============================================================================
# 1. PROJECT PATHS
#
# Change PROJECT_DIR once if your repository is located elsewhere.
# =============================================================================

PROJECT_DIR <- "."

MASTER_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "processed",
  "ANM_ROBUST_MR_GENES_FOR_COLOC.xlsx"
)

MASTER_SHEET <- "COLOC_ready"

ANM_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "gwas",
  "ANM_GCST90292553_GRCh38.txt"
)

OUT_DIR <- file.path(
  PROJECT_DIR,
  "data",
  "processed",
  "GTEx_coloc_inputs"
)

LOG_FILE <- file.path(
  OUT_DIR,
  "GTEx_coloc_input_preparation_log.xlsx"
)


dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 2. ANALYSIS SETTINGS
# =============================================================================

WINDOW_BP <- 500000L

MAX_RETRIES <- 5L

RETRY_WAIT_SECONDS <- 5L


# =============================================================================
# 3. GTEx v8 TABIX RESOURCE
#
# GTEx v8 tissue-specific cis-eQTL summary statistics are distributed as
# tabix-indexed files through the eQTL Catalogue.
#
# The repository file tabix_ftp_paths_imported.tsv documents the imported
# dataset locations.
# =============================================================================

GTEX_BASE <- paste0(
  "ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/",
  "imported/GTEx_V8/ge/"
)


GTEX_COLUMNS <- c(
  "variant",
  "r2",
  "pvalue",
  "molecular_trait_object_id",
  "molecular_trait_id",
  "maf",
  "gene_id",
  "median_tpm",
  "beta",
  "se",
  "an",
  "ac",
  "chromosome",
  "position",
  "ref",
  "alt",
  "type",
  "rsid"
)


# =============================================================================
# 4. READ MR-PRIORITIZED GENE MANIFEST
#
# Required fields:
#   Gene
#   Ensembl_ID
#   Representative_SNP
#   MR_Best_shared_tissue
#
# Colocalization is deliberately restricted to the MR-best tissue for each
# gene to provide a systematic validation layer for the reported MR signal.
# =============================================================================

message("Reading MR-prioritized gene manifest...")

master <- as.data.table(
  read.xlsx(
    MASTER_FILE,
    sheet = MASTER_SHEET
  )
)


required_master <- c(
  "Gene",
  "Ensembl_ID",
  "Representative_SNP",
  "MR_Best_shared_tissue"
)


missing_master <- setdiff(
  required_master,
  names(master)
)


if (length(missing_master) > 0) {

  stop(
    "Missing required manifest columns: ",
    paste(
      missing_master,
      collapse = ", "
    )
  )
}


message(
  "Genes available for preparation: ",
  nrow(master)
)


# =============================================================================
# 5. READ ANM GRCh38 GWAS
# =============================================================================

message("Reading ANM GRCh38 GWAS...")

anm <- fread(
  ANM_FILE
)


required_anm <- c(
  "chromosome",
  "base_pair_location",
  "rs_id"
)


missing_anm <- setdiff(
  required_anm,
  names(anm)
)


if (length(missing_anm) > 0) {

  stop(
    "ANM GWAS missing required columns: ",
    paste(
      missing_anm,
      collapse = ", "
    )
  )
}


anm[, chromosome := as.character(chromosome)]

anm[, base_pair_location := as.integer(base_pair_location)]


# =============================================================================
# 6. OUTPUT LOG
# =============================================================================

log_list <- vector(
  "list",
  nrow(master)
)


# =============================================================================
# 7. LOOP OVER GENES
# =============================================================================

for (i in seq_len(nrow(master))) {

  gene <- as.character(
    master$Gene[i]
  )

  ensg <- as.character(
    master$Ensembl_ID[i]
  )

  representative_snp <- as.character(
    master$Representative_SNP[i]
  )

  tissue <- as.character(
    master$MR_Best_shared_tissue[i]
  )


  message("")
  message("============================================================")
  message("Preparing ", i, " of ", nrow(master))
  message("Gene: ", gene)
  message("Tissue: ", tissue)
  message("Representative SNP: ", representative_snp)
  message("============================================================")


  # ---------------------------------------------------------------------------
  # Locate the representative ANM SNP in GRCh38.
  # ---------------------------------------------------------------------------

  anchor <- anm[
    rs_id == representative_snp
  ]


  if (nrow(anchor) == 0) {

    log_list[[i]] <- data.table(
      Gene = gene,
      Ensembl_ID = ensg,
      Tissue = tissue,
      Representative_SNP = representative_snp,
      Chromosome = NA_character_,
      Center_GRCh38 = NA_integer_,
      Region_start = NA_integer_,
      Region_end = NA_integer_,
      GTEx_rows_region = NA_integer_,
      GTEx_rows_gene = NA_integer_,
      Attempts = 0L,
      Status = "REP_SNP_NOT_FOUND",
      Output_file = NA_character_
    )

    next
  }


  anchor <- anchor[1]


  chr <- as.character(
    anchor$chromosome
  )


  center_pos <- as.integer(
    anchor$base_pair_location
  )


  region_start <- max(
    1L,
    center_pos - WINDOW_BP
  )


  region_end <- center_pos + WINDOW_BP


  region_string <- paste0(
    chr,
    ":",
    region_start,
    "-",
    region_end
  )


  message(
    "GRCh38 region: ",
    region_string
  )


  # ---------------------------------------------------------------------------
  # Local output.
  # ---------------------------------------------------------------------------

  output_file <- file.path(
    OUT_DIR,
    paste0(
      gene,
      "_",
      tissue,
      "_GTEx_v8_GRCh38_symmetric.tsv"
    )
  )


  # Restart-safe behavior.
  if (
    file.exists(output_file) &&
    file.info(output_file)$size > 0
  ) {

    message("Existing prepared file found; skipping.")

    existing_n <- try(
      nrow(
        fread(
          output_file,
          select = 1
        )
      ),
      silent = TRUE
    )


    if (inherits(existing_n, "try-error")) {
      existing_n <- NA_integer_
    }


    log_list[[i]] <- data.table(
      Gene = gene,
      Ensembl_ID = ensg,
      Tissue = tissue,
      Representative_SNP = representative_snp,
      Chromosome = chr,
      Center_GRCh38 = center_pos,
      Region_start = region_start,
      Region_end = region_end,
      GTEx_rows_region = NA_integer_,
      GTEx_rows_gene = existing_n,
      Attempts = 0L,
      Status = "EXISTING_OK",
      Output_file = output_file
    )

    next
  }


  # ---------------------------------------------------------------------------
  # GTEx tissue file.
  # ---------------------------------------------------------------------------

  gtex_url <- paste0(
    GTEX_BASE,
    tissue,
    ".tsv.gz"
  )


  # ---------------------------------------------------------------------------
  # Regional tabix query with retries.
  # ---------------------------------------------------------------------------

  eq <- NULL

  query_success <- FALSE

  attempts_used <- 0L


  for (attempt in seq_len(MAX_RETRIES)) {

    attempts_used <- attempt


    message(
      "GTEx query attempt ",
      attempt,
      " of ",
      MAX_RETRIES,
      "..."
    )


    eq_try <- try(
      seqminer::tabix.read.table(
        tabixFile = gtex_url,
        tabixRange = region_string,
        stringsAsFactors = FALSE
      ),
      silent = TRUE
    )


    if (
      !inherits(eq_try, "try-error") &&
      !is.null(eq_try) &&
      nrow(eq_try) > 0 &&
      ncol(eq_try) == 18
    ) {

      eq <- as.data.table(
        eq_try
      )

      query_success <- TRUE

      break
    }


    if (attempt < MAX_RETRIES) {

      message(
        "Query failed; retrying in ",
        RETRY_WAIT_SECONDS,
        " seconds..."
      )

      Sys.sleep(
        RETRY_WAIT_SECONDS
      )
    }
  }


  if (!query_success) {

    log_list[[i]] <- data.table(
      Gene = gene,
      Ensembl_ID = ensg,
      Tissue = tissue,
      Representative_SNP = representative_snp,
      Chromosome = chr,
      Center_GRCh38 = center_pos,
      Region_start = region_start,
      Region_end = region_end,
      GTEx_rows_region = 0L,
      GTEx_rows_gene = 0L,
      Attempts = attempts_used,
      Status = "GTEX_QUERY_FAILED",
      Output_file = NA_character_
    )

    next
  }


  # ---------------------------------------------------------------------------
  # Assign eQTL Catalogue column names.
  # ---------------------------------------------------------------------------

  setnames(
    eq,
    GTEX_COLUMNS
  )


  eq[, chromosome := as.character(chromosome)]

  eq[, position := as.integer(position)]


  n_region <- nrow(eq)


  # ---------------------------------------------------------------------------
  # Retain the target gene.
  #
  # These are regional cis-eQTL summary statistics, not merely the
  # conditionally independent eQTL instruments used for MR.
  # ---------------------------------------------------------------------------

  eq_gene <- eq[
    gene_id == ensg
  ]


  n_gene <- nrow(
    eq_gene
  )


  message(
    "Regional GTEx rows: ",
    n_region
  )

  message(
    "Target-gene rows retained: ",
    n_gene
  )


  if (n_gene == 0) {

    log_list[[i]] <- data.table(
      Gene = gene,
      Ensembl_ID = ensg,
      Tissue = tissue,
      Representative_SNP = representative_snp,
      Chromosome = chr,
      Center_GRCh38 = center_pos,
      Region_start = region_start,
      Region_end = region_end,
      GTEx_rows_region = n_region,
      GTEx_rows_gene = 0L,
      Attempts = attempts_used,
      Status = "NO_GENE_ROWS",
      Output_file = NA_character_
    )

    next
  }


  # ---------------------------------------------------------------------------
  # Save gene-specific regional summary statistics.
  # ---------------------------------------------------------------------------

  fwrite(
    eq_gene,
    output_file,
    sep = "\t",
    quote = FALSE
  )


  log_list[[i]] <- data.table(
    Gene = gene,
    Ensembl_ID = ensg,
    Tissue = tissue,
    Representative_SNP = representative_snp,
    Chromosome = chr,
    Center_GRCh38 = center_pos,
    Region_start = region_start,
    Region_end = region_end,
    GTEx_rows_region = n_region,
    GTEx_rows_gene = n_gene,
    Attempts = attempts_used,
    Status = "OK",
    Output_file = output_file
  )


  # ---------------------------------------------------------------------------
  # Checkpoint log.
  # ---------------------------------------------------------------------------

  current_log <- rbindlist(
    log_list[
      !vapply(
        log_list,
        is.null,
        logical(1)
      )
    ],
    fill = TRUE
  )


  write.xlsx(
    current_log,
    LOG_FILE,
    sheetName = "Preparation_log",
    overwrite = TRUE
  )
}


# =============================================================================
# 8. FINAL LOG
# =============================================================================

final_log <- rbindlist(
  log_list[
    !vapply(
      log_list,
      is.null,
      logical(1)
    )
  ],
  fill = TRUE
)


status_summary <- final_log[
  ,
  .N,
  by = Status
][
  order(Status)
]


wb <- createWorkbook()

addWorksheet(
  wb,
  "Preparation_log"
)

writeData(
  wb,
  "Preparation_log",
  final_log,
  withFilter = TRUE
)

freezePane(
  wb,
  "Preparation_log",
  firstRow = TRUE
)


addWorksheet(
  wb,
  "Status_summary"
)

writeData(
  wb,
  "Status_summary",
  status_summary
)


saveWorkbook(
  wb,
  LOG_FILE,
  overwrite = TRUE
)


message("")
message("============================================================")
message("GTEx COLOCALIZATION INPUT PREPARATION COMPLETE")
message("============================================================")

print(
  status_summary
)

message("")
message("Prepared files:")
message(OUT_DIR)

message("")
message("Preparation log:")
message(LOG_FILE)

message("============================================================")