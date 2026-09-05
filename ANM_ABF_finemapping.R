# ============================================================
# Age at Natural Menopause (ANM) Fine-Mapping Analysis
# ============================================================
#
# Author: Suzan Farhang-Sardroodi
#
# This code was developed by Suzan Farhang-Sardroodi for the
# fine-mapping analyses of age at natural menopause (ANM)
# reported in this study.
#
# Please cite the associated publication when using or adapting
# this code.
#
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})


# =========================================================
# INSTALL AND LOAD OPENXLSX
# =========================================================

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages(
    "openxlsx",
    repos = "https://cloud.r-project.org"
  )
}

library(openxlsx)


# =========================================================
# 1. USER INPUTS
# =========================================================

# Full ANM GWAS summary-statistics file
ANM_FILE <- file.path(
  path.expand("~/Desktop"),
  "ANM_GCST90292553.txt"
)


# FUMA genomic-risk-loci file
RISK_LOCI_FILE <- file.path(
  path.expand("~/Desktop"),
  "GenomicRiskLoci.txt"
)


# New output folder for the FUMA-locus-based analysis
OUT_DIR <- file.path(
  path.expand("~/Desktop"),
  "ANM_ABF_finemapping_FUMA_loci"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ANM GWAS discovery sample size
ANM_SAMPLE_SIZE <- 173424


# Credible-set probability
CRED_P <- 0.90


# Prior probability that any individual SNP is causal
PRIOR_P <- 1e-4


# Minimum number of valid SNPs required for a locus
MIN_SNPS <- 10L


# ---------------------------------------------------------
# Complete output files
# ---------------------------------------------------------

OUT_ALL_CSV <- file.path(
  OUT_DIR,
  "ANM_ABF_FUMA_loci_all_SNPs.csv"
)

OUT_CS_CSV <- file.path(
  OUT_DIR,
  "ANM_ABF_FUMA_loci_90pct_credible_sets.csv"
)

OUT_SUMMARY_CSV <- file.path(
  OUT_DIR,
  "ANM_ABF_FUMA_locus_summary.csv"
)

OUT_SKIPPED_CSV <- file.path(
  OUT_DIR,
  "ANM_ABF_FUMA_skipped_loci.csv"
)


# ---------------------------------------------------------
# Clean Excel workbook
# ---------------------------------------------------------

OUT_XLSX <- file.path(
  OUT_DIR,
  "ANM_ABF_finemapping_FUMA_Final.xlsx"
)


# =========================================================
# 2. READ AND STANDARDIZE THE ANM GWAS
# =========================================================

read_anm <- function(path) {

  if (!file.exists(path)) {
    stop(
      "ANM GWAS file not found: ",
      path
    )
  }


  message("Reading ANM GWAS file...")


  dt <- fread(
    path,
    showProgress = TRUE
  )


  required_columns <- c(
    "chromosome",
    "base_pair_location",
    "effect_allele",
    "other_allele",
    "beta",
    "standard_error",
    "effect_allele_frequency",
    "p_value",
    "rs_id"
  )


  missing_columns <- setdiff(
    required_columns,
    names(dt)
  )


  if (length(missing_columns) > 0) {
    stop(
      "Missing required ANM columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }


  dt <- dt[
    ,
    .(
      snp = as.character(rs_id),

      chr = as.character(
        chromosome
      ),

      position = as.integer(
        base_pair_location
      ),

      a1 = toupper(
        as.character(
          effect_allele
        )
      ),

      a2 = toupper(
        as.character(
          other_allele
        )
      ),

      beta = as.numeric(
        beta
      ),

      se = as.numeric(
        standard_error
      ),

      eaf = as.numeric(
        effect_allele_frequency
      ),

      p = as.numeric(
        p_value
      )
    )
  ]


  # Remove optional chromosome prefix
  dt[
    ,
    chr := sub(
      "^chr",
      "",
      chr,
      ignore.case = TRUE
    )
  ]


  # Represent chromosome X internally as chromosome 23
  dt[
    toupper(chr) == "X",
    chr := "23"
  ]


  dt[
    ,
    chr_num := suppressWarnings(
      as.integer(chr)
    )
  ]


  # Calculate minor-allele frequency
  dt[
    ,
    maf := pmin(
      eaf,
      1 - eaf
    )
  ]


  # Retain valid GWAS variants
  dt <- dt[
    !is.na(snp) &
    snp != "" &

    !is.na(chr_num) &
    chr_num >= 1 &
    chr_num <= 23 &

    !is.na(position) &
    position > 0 &

    is.finite(beta) &

    is.finite(se) &
    se > 0 &

    is.finite(eaf) &
    eaf > 0 &
    eaf < 1 &

    is.finite(maf) &
    maf > 0 &
    maf <= 0.5
  ]


  # If duplicated rsIDs exist, retain the row
  # with the smallest standard error
  setorder(
    dt,
    snp,
    se
  )


  duplicated_count <- sum(
    duplicated(dt$snp)
  )


  if (duplicated_count > 0) {

    warning(
      duplicated_count,
      " duplicated SNP IDs were found. ",
      "For each duplicate, the row with the ",
      "smallest standard error will be retained."
    )


    dt <- dt[
      !duplicated(snp)
    ]
  }


  # Calculate GWAS Z-scores
  dt[
    ,
    z := beta / se
  ]


  setorder(
    dt,
    chr_num,
    position
  )


  message(
    "ANM variants retained: ",
    format(
      nrow(dt),
      big.mark = ","
    )
  )


  dt
}


# =========================================================
# 3. READ AND STANDARDIZE FUMA GENOMIC RISK LOCI
# =========================================================

read_fuma_loci <- function(path) {

  if (!file.exists(path)) {
    stop(
      "FUMA GenomicRiskLoci file not found: ",
      path
    )
  }


  message("Reading FUMA genomic risk loci...")


  loci <- fread(
    path,
    header = TRUE,
    showProgress = TRUE
  )


  required_columns <- c(
    "GenomicLocus",
    "uniqID",
    "rsID",
    "chr",
    "pos",
    "p",
    "start",
    "end",
    "nSNPs",
    "nGWASSNPs",
    "nIndSigSNPs",
    "IndSigSNPs",
    "nLeadSNPs",
    "LeadSNPs"
  )


  missing_columns <- setdiff(
    required_columns,
    names(loci)
  )


  if (length(missing_columns) > 0) {
    stop(
      "Missing required FUMA columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }


  loci <- loci[
    ,
    .(
      genomic_locus =
        as.integer(
          GenomicLocus
        ),

      fuma_unique_id =
        as.character(
          uniqID
        ),

      representative_snp =
        as.character(
          rsID
        ),

      chr_raw =
        as.character(
          chr
        ),

      representative_position =
        as.integer(
          pos
        ),

      representative_p =
        as.numeric(
          p
        ),

      region_start =
        as.integer(
          start
        ),

      region_end =
        as.integer(
          end
        ),

      fuma_n_snps =
        as.integer(
          nSNPs
        ),

      fuma_n_gwas_snps =
        as.integer(
          nGWASSNPs
        ),

      n_independent_significant_snps =
        as.integer(
          nIndSigSNPs
        ),

      independent_significant_snps =
        as.character(
          IndSigSNPs
        ),

      n_lead_snps =
        as.integer(
          nLeadSNPs
        ),

      lead_snps =
        as.character(
          LeadSNPs
        )
    )
  ]


  # Remove an optional chromosome prefix
  loci[
    ,
    chr_raw := sub(
      "^chr",
      "",
      chr_raw,
      ignore.case = TRUE
    )
  ]


  # Standardize chromosome X
  loci[
    toupper(chr_raw) == "X",
    chr_raw := "23"
  ]


  loci[
    ,
    chr_num := suppressWarnings(
      as.integer(chr_raw)
    )
  ]


  # Retain valid genomic risk loci
  loci <- loci[
    !is.na(genomic_locus) &

    !is.na(chr_num) &
    chr_num >= 1 &
    chr_num <= 23 &

    !is.na(region_start) &
    !is.na(region_end) &

    region_start > 0 &
    region_end >= region_start
  ]


  # Create readable chromosome labels
  loci[
    ,
    chr_label := fifelse(
      chr_num == 23,
      "X",
      as.character(
        chr_num
      )
    )
  ]


  # Create one unique analysis label per FUMA locus
  loci[
    ,
    locus := sprintf(
      "FUMA_locus_%03d_chr%s_%d_%d",
      genomic_locus,
      chr_label,
      region_start,
      region_end
    )
  ]


  # Remove duplicated genomic-locus records, if present
  duplicated_loci <- sum(
    duplicated(
      loci$genomic_locus
    )
  )


  if (duplicated_loci > 0) {

    warning(
      duplicated_loci,
      " duplicated FUMA locus numbers were found. ",
      "The first row for each locus will be retained."
    )


    loci <- loci[
      !duplicated(
        genomic_locus
      )
    ]
  }


  setorder(
    loci,
    chr_num,
    region_start
  )


  message(
    "FUMA genomic risk loci retained: ",
    format(
      nrow(loci),
      big.mark = ","
    )
  )


  loci
}


# =========================================================
# 4. BUILD ONE ABF DATASET PER FUMA RISK LOCUS
# =========================================================

build_anm_inputs <- function(
  gwas,
  loci,
  sample_size = 173424,
  min_snps = 10L
) {

  inputs <- list()
  metadata <- list()
  skipped <- list()


  for (i in seq_len(nrow(loci))) {

    current_locus <- loci[i]


    locus_name <-
      current_locus$locus

    locus_number <-
      current_locus$genomic_locus

    locus_chr <-
      current_locus$chr_num

    region_start <-
      current_locus$region_start

    region_end <-
      current_locus$region_end


    message(
      "Preparing ",
      locus_name,
      "..."
    )


    # Extract all GWAS SNPs inside the exact
    # FUMA-defined start and end positions
    locus_snps <- gwas[
      chr_num == locus_chr &
      position >= region_start &
      position <= region_end
    ]


    # Retain variants suitable for ABF fine-mapping
    locus_snps <- locus_snps[
      is.finite(beta) &

      is.finite(se) &
      se > 0 &

      is.finite(maf) &
      maf > 0 &
      maf <= 0.5 &

      !is.na(snp) &
      snp != "" &

      !is.na(position) &
      position > 0
    ]


    # Skip loci with too few valid GWAS variants
    if (nrow(locus_snps) < min_snps) {

      warning(
        "Skipping ",
        locus_name,
        ": only ",
        nrow(locus_snps),
        " valid GWAS SNPs were found."
      )


      skipped[[locus_name]] <- data.table(
        locus = locus_name,

        genomic_locus =
          locus_number,

        chr =
          locus_chr,

        region_start =
          region_start,

        region_end =
          region_end,

        representative_snp =
          current_locus$representative_snp,

        number_of_valid_snps =
          nrow(locus_snps),

        reason = paste0(
          "Fewer than ",
          min_snps,
          " valid SNPs"
        )
      )


      next
    }


    # ABF dataset for a quantitative phenotype
    inp <- list(
      beta =
        locus_snps$beta,

      varbeta =
        locus_snps$se^2,

      MAF =
        locus_snps$maf,

      N =
        sample_size,

      snp =
        locus_snps$snp,

      position =
        locus_snps$position,

      type =
        "quant"
    )


    inputs[[locus_name]] <- inp


    # Retain locus-level metadata
    metadata[[locus_name]] <- data.table(
      locus =
        locus_name,

      genomic_locus =
        locus_number,

      chr =
        locus_chr,

      chr_label = ifelse(
        locus_chr == 23,
        "X",
        as.character(
          locus_chr
        )
      ),

      fuma_unique_id =
        current_locus$fuma_unique_id,

      representative_snp =
        current_locus$representative_snp,

      representative_position =
        current_locus$representative_position,

      representative_p =
        current_locus$representative_p,

      region_start =
        region_start,

      region_end =
        region_end,

      region_width_bp =
        region_end -
        region_start +
        1L,

      fuma_n_snps =
        current_locus$fuma_n_snps,

      fuma_n_gwas_snps =
        current_locus$fuma_n_gwas_snps,

      n_independent_significant_snps =
        current_locus$n_independent_significant_snps,

      independent_significant_snps =
        current_locus$independent_significant_snps,

      n_lead_snps =
        current_locus$n_lead_snps,

      lead_snps =
        current_locus$lead_snps,

      number_of_snps_analyzed =
        nrow(locus_snps)
    )
  }


  metadata_dt <- if (
    length(metadata) > 0
  ) {

    rbindlist(
      metadata,
      use.names = TRUE,
      fill = TRUE
    )

  } else {

    data.table()
  }


  skipped_dt <- if (
    length(skipped) > 0
  ) {

    rbindlist(
      skipped,
      use.names = TRUE,
      fill = TRUE
    )

  } else {

    data.table()
  }


  list(
    inputs = inputs,
    metadata = metadata_dt,
    skipped = skipped_dt
  )
}


# =========================================================
# 5. RUN coloc::finemap.abf()
# =========================================================

run_abf_finemapping <- function(
  input_list,
  prior = 1e-4
) {

  results <- vector(
    mode = "list",
    length = length(
      input_list
    )
  )


  names(results) <- names(
    input_list
  )


  for (
    locus_name in names(
      input_list
    )
  ) {

    message(
      "Fine-mapping ",
      locus_name,
      "..."
    )


    result <- tryCatch(
      {

        finemap.abf(
          dataset =
            input_list[[locus_name]],

          p1 =
            prior
        )

      },
      error = function(e) {

        warning(
          "Fine-mapping failed for ",
          locus_name,
          ": ",
          conditionMessage(e)
        )


        NULL
      }
    )


    results[[locus_name]] <-
      result
  }


  # Remove failed loci
  results[
    !vapply(
      results,
      is.null,
      logical(1)
    )
  ]
}


# =========================================================
# 6. COMBINE ALL LOCUS RESULTS
# =========================================================

finemap_to_dt <- function(
  fm_list
) {

  if (length(fm_list) == 0) {
    return(
      data.table()
    )
  }


  rbindlist(
    lapply(
      names(fm_list),
      function(locus_name) {

        x <- as.data.table(
          fm_list[[locus_name]]
        )


        x[
          ,
          locus := locus_name
        ]


        x
      }
    ),

    use.names = TRUE,
    fill = TRUE
  )
}


# =========================================================
# 7. CREATE LOCUS-SPECIFIC 90% CREDIBLE SETS
# =========================================================

add_credible_set <- function(
  dt,
  credible_probability = 0.90
) {

  setDT(dt)


  posterior_column <- if (
    "SNP.PP" %in% names(dt)
  ) {

    "SNP.PP"

  } else if (
    "SNP.PP.H4" %in% names(dt)
  ) {

    "SNP.PP.H4"

  } else {

    stop(
      "No posterior-probability column was found."
    )
  }


  dt[
    ,
    credible_set := 0L
  ]


  dt[
    ,
    cumulative_PP := NA_real_
  ]


  result <- dt[
    ,
    {

      x <- copy(.SD)


      # Rank variants from highest to lowest PIP
      setorderv(
        x,
        posterior_column,
        order = -1L
      )


      # Calculate cumulative posterior probability
      x[
        ,
        cumulative_PP :=
          cumsum(
            get(
              posterior_column
            )
          )
      ]


      # Include the first variant that causes the
      # cumulative PIP to reach or exceed 90%
      k <- which(
        x$cumulative_PP >=
        credible_probability
      )[1]


      if (is.na(k)) {
        k <- nrow(x)
      }


      x[
        seq_len(k),
        credible_set := 1L
      ]


      x
    },

    by = locus
  ]


  setnames(
    result,
    "credible_set",
    "CS90"
  )


  result
}


# =========================================================
# 8. MERGE ORIGINAL ANM GWAS INFORMATION
# =========================================================

merge_anm_information <- function(
  fm,
  gwas
) {

  annotation <- unique(
    gwas[
      ,
      .(
        snp,

        chr =
          chr_num,

        position_gwas =
          position,

        a1,
        a2,
        eaf,
        maf,

        beta_original =
          beta,

        se_original =
          se,

        p_original =
          p
      )
    ],

    by = "snp"
  )


  merge(
    fm,
    annotation,
    by = "snp",
    all.x = TRUE,
    sort = FALSE
  )
}


# =========================================================
# 9. MERGE FUMA LOCUS INFORMATION
# =========================================================

merge_fuma_information <- function(
  fm,
  metadata
) {

  locus_annotation <- unique(
    metadata[
      ,
      .(
        locus,
        genomic_locus,
        fuma_unique_id,
        representative_snp,
        representative_position,
        representative_p,
        region_start,
        region_end,
        region_width_bp,
        n_independent_significant_snps,
        independent_significant_snps,
        n_lead_snps,
        lead_snps
      )
    ],

    by = "locus"
  )


  merge(
    fm,
    locus_annotation,
    by = "locus",
    all.x = TRUE,
    sort = FALSE
  )
}


# =========================================================
# 10. CREATE LOCUS-LEVEL SUMMARY
# =========================================================

create_locus_summary <- function(
  fm,
  metadata
) {

  summary <- fm[
    snp != "null",
    {

      ordered <- .SD[
        order(
          -SNP.PP
        )
      ]


      cs <- ordered[
        CS90 == 1L
      ]


      .(
        top_posterior_snp =
          ordered$snp[1],

        top_SNP_PP =
          ordered$SNP.PP[1],

        top_snp_position =
          ordered$position[1],

        top_snp_p_original =
          ordered$p_original[1],

        number_of_finemapped_snps =
          .N,

        CS90_size =
          nrow(cs),

        CS90_total_PP =
          sum(
            cs$SNP.PP,
            na.rm = TRUE
          ),

        CS90_min_PP = if (
          nrow(cs) > 0
        ) {

          min(
            cs$SNP.PP,
            na.rm = TRUE
          )

        } else {

          NA_real_
        },

        CS90_max_PP = if (
          nrow(cs) > 0
        ) {

          max(
            cs$SNP.PP,
            na.rm = TRUE
          )

        } else {

          NA_real_
        }
      )
    },

    by = locus
  ]


  merge(
    metadata,
    summary,
    by = "locus",
    all.x = TRUE,
    sort = FALSE
  )
}


# =========================================================
# 11. RUN THE COMPLETE FINE-MAPPING ANALYSIS
# =========================================================

anm <- read_anm(
  ANM_FILE
)


anm_loci <- read_fuma_loci(
  RISK_LOCI_FILE
)


message(
  "Building fine-mapping datasets using ",
  "the exact FUMA start/end boundaries..."
)


anm_pack <- build_anm_inputs(
  gwas = anm,
  loci = anm_loci,
  sample_size = ANM_SAMPLE_SIZE,
  min_snps = MIN_SNPS
)


message(
  "FUMA loci prepared for fine-mapping: ",
  length(
    anm_pack$inputs
  )
)


if (
  length(
    anm_pack$inputs
  ) == 0
) {

  stop(
    "No FUMA loci contained enough valid SNPs ",
    "for fine-mapping."
  )
}


message(
  "Running coloc::finemap.abf()..."
)


anm_fm_list <- run_abf_finemapping(
  input_list =
    anm_pack$inputs,

  prior =
    PRIOR_P
)


anm_raw <- finemap_to_dt(
  anm_fm_list
)


if (nrow(anm_raw) == 0) {
  stop(
    "No fine-mapping results were produced."
  )
}


# Remove the null model before building SNP credible sets
anm_raw <- anm_raw[
  snp != "null"
]


# Create locus-specific 90% credible sets
anm_out <- add_credible_set(
  anm_raw,
  credible_probability = CRED_P
)


# Add original GWAS annotations
anm_out <- merge_anm_information(
  anm_out,
  anm
)


# Add FUMA locus annotations
anm_out <- merge_fuma_information(
  anm_out,
  anm_pack$metadata
)


# =========================================================
# 12. ARRANGE COMPLETE SNP-LEVEL OUTPUT
# =========================================================

preferred_columns <- c(
  "snp",
  "locus",
  "genomic_locus",
  "chr",
  "position",
  "position_gwas",
  "a1",
  "a2",
  "eaf",
  "maf",
  "beta_original",
  "se_original",
  "p_original",
  "V.",
  "z.",
  "lABF.",
  "prior",
  "SNP.PP",
  "cumulative_PP",
  "CS90",
  "fuma_unique_id",
  "representative_snp",
  "representative_position",
  "representative_p",
  "region_start",
  "region_end",
  "region_width_bp",
  "n_independent_significant_snps",
  "independent_significant_snps",
  "n_lead_snps",
  "lead_snps"
)


available_columns <- intersect(
  preferred_columns,
  names(anm_out)
)


other_columns <- setdiff(
  names(anm_out),
  available_columns
)


setcolorder(
  anm_out,
  c(
    available_columns,
    other_columns
  )
)


# =========================================================
# 13. CREATE COMPLETE RESULT TABLES
# =========================================================

anm_cs90 <- anm_out[
  CS90 == 1L
]


anm_summary <- create_locus_summary(
  fm = anm_out,
  metadata = anm_pack$metadata
)


setorder(
  anm_summary,
  chr,
  region_start
)


setorder(
  anm_out,
  chr,
  position,
  -SNP.PP
)


setorder(
  anm_cs90,
  chr,
  position,
  -SNP.PP
)


setorder(
  anm_loci,
  chr_num,
  region_start
)


# =========================================================
# 14. WRITE COMPLETE CSV OUTPUTS
# =========================================================

message(
  "Writing complete SNP-level CSV..."
)


fwrite(
  anm_out,
  OUT_ALL_CSV
)


message(
  "Writing 90% credible-set CSV..."
)


fwrite(
  anm_cs90,
  OUT_CS_CSV
)


message(
  "Writing locus-summary CSV..."
)


fwrite(
  anm_summary,
  OUT_SUMMARY_CSV
)


if (
  nrow(
    anm_pack$skipped
  ) > 0
) {

  fwrite(
    anm_pack$skipped,
    OUT_SKIPPED_CSV
  )
}


# =========================================================
# 15. CREATE CLEAN JENNIE-STYLE EXCEL TABLES
# =========================================================

excel_columns <- c(
  "snp",
  "locus",
  "V.",
  "z.",
  "lABF.",
  "position",
  "prior",
  "SNP.PP",
  "CS90",
  "maf",
  "chr",
  "a1",
  "a2"
)


missing_excel_columns <- setdiff(
  excel_columns,
  names(anm_out)
)


if (
  length(
    missing_excel_columns
  ) > 0
) {

  stop(
    "The following columns required for the Excel ",
    "workbook are missing: ",
    paste(
      missing_excel_columns,
      collapse = ", "
    )
  )
}


# ---------------------------------------------------------
# All fine-mapped SNP-locus results
# ---------------------------------------------------------

excel_all <- copy(
  anm_out[
    ,
    ..excel_columns
  ]
)


setnames(
  excel_all,
  "CS90",
  "90_credset"
)


# Ensure the credible-set column is integer 0/1
excel_all[
  ,
  `90_credset` :=
    as.integer(
      `90_credset`
    )
]


# ---------------------------------------------------------
# 90% credible-set SNPs
# ---------------------------------------------------------

excel_cs90 <- copy(
  excel_all[
    `90_credset` == 1L
  ]
)


# ---------------------------------------------------------
# Highest-PIP SNP from each FUMA genomic risk locus
# ---------------------------------------------------------

excel_top <- excel_all[
  order(
    locus,
    -SNP.PP
  ),
  .SD[1],
  by = locus
]


# ---------------------------------------------------------
# SNPs with PIP > 0.03
# ---------------------------------------------------------

excel_pip <- copy(
  excel_all[
    SNP.PP > 0.03
  ]
)


# ---------------------------------------------------------
# Credible-set SNPs with PIP > 0.03
# ---------------------------------------------------------

excel_cs90_pip <- copy(
  excel_all[
    `90_credset` == 1L &
    SNP.PP > 0.03
  ]
)


# Sort the clean tables
setorder(
  excel_all,
  chr,
  position,
  -SNP.PP
)


setorder(
  excel_cs90,
  chr,
  position,
  -SNP.PP
)


setorder(
  excel_top,
  chr,
  position
)


setorder(
  excel_pip,
  chr,
  position,
  -SNP.PP
)


setorder(
  excel_cs90_pip,
  chr,
  position,
  -SNP.PP
)


# =========================================================
# 16. WRITE THE CLEAN EXCEL WORKBOOK
# =========================================================

message(
  "Writing clean Excel workbook: ",
  OUT_XLSX
)


wb <- createWorkbook()


sheet_names <- c(
  "All_finemapped_SNPs",
  "90pct_credible_sets",
  "Top_SNP_per_locus",
  "PIP_above_0.03",
  "CS90_PIP_above_0.03"
)


for (
  sheet_name in sheet_names
) {

  addWorksheet(
    wb,
    sheet_name,
    gridLines = TRUE
  )
}


# ---------------------------------------------------------
# Write data tables
# ---------------------------------------------------------

writeDataTable(
  wb,
  "All_finemapped_SNPs",
  excel_all,
  tableStyle = "TableStyleMedium2"
)


writeDataTable(
  wb,
  "90pct_credible_sets",
  excel_cs90,
  tableStyle = "TableStyleMedium2"
)


writeDataTable(
  wb,
  "Top_SNP_per_locus",
  excel_top,
  tableStyle = "TableStyleMedium2"
)


writeDataTable(
  wb,
  "PIP_above_0.03",
  excel_pip,
  tableStyle = "TableStyleMedium2"
)


writeDataTable(
  wb,
  "CS90_PIP_above_0.03",
  excel_cs90_pip,
  tableStyle = "TableStyleMedium2"
)


# ---------------------------------------------------------
# Freeze the header row and retain gridlines
# ---------------------------------------------------------

for (
  sheet_name in sheet_names
) {

  freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )


  showGridLines(
    wb,
    sheet = sheet_name,
    showGridLines = TRUE
  )
}


# ---------------------------------------------------------
# Set practical column widths
# ---------------------------------------------------------

column_widths <- c(
  16,  # snp
  36,  # locus
  15,  # V.
  15,  # z.
  16,  # lABF.
  14,  # position
  12,  # prior
  16,  # SNP.PP
  15,  # 90_credset
  12,  # maf
  10,  # chr
  8,   # a1
  8    # a2
)


for (
  sheet_name in sheet_names
) {

  setColWidths(
    wb,
    sheet = sheet_name,
    cols = seq_along(
      column_widths
    ),
    widths = column_widths
  )
}


# ---------------------------------------------------------
# Create numeric styles
# ---------------------------------------------------------

scientific_style <- createStyle(
  numFmt = "0.000000E+00"
)


decimal_style <- createStyle(
  numFmt = "0.000000000"
)


integer_style <- createStyle(
  numFmt = "0"
)


prior_style <- createStyle(
  numFmt = "0.0000"
)


maf_style <- createStyle(
  numFmt = "0.000000"
)


# ---------------------------------------------------------
# Apply numeric formatting to each sheet
# ---------------------------------------------------------

excel_tables <- list(
  All_finemapped_SNPs =
    excel_all,

  `90pct_credible_sets` =
    excel_cs90,

  Top_SNP_per_locus =
    excel_top,

  `PIP_above_0.03` =
    excel_pip,

  `CS90_PIP_above_0.03` =
    excel_cs90_pip
)


for (
  sheet_name in names(
    excel_tables
  )
) {

  number_of_rows <- nrow(
    excel_tables[[sheet_name]]
  )


  if (number_of_rows == 0) {
    next
  }


  data_rows <- 2:(
    number_of_rows + 1
  )


  # V.
  addStyle(
    wb,
    sheet = sheet_name,
    style = scientific_style,
    rows = data_rows,
    cols = 3,
    gridExpand = TRUE,
    stack = TRUE
  )


  # z.
  addStyle(
    wb,
    sheet = sheet_name,
    style = decimal_style,
    rows = data_rows,
    cols = 4,
    gridExpand = TRUE,
    stack = TRUE
  )


  # lABF.
  addStyle(
    wb,
    sheet = sheet_name,
    style = decimal_style,
    rows = data_rows,
    cols = 5,
    gridExpand = TRUE,
    stack = TRUE
  )


  # position
  addStyle(
    wb,
    sheet = sheet_name,
    style = integer_style,
    rows = data_rows,
    cols = 6,
    gridExpand = TRUE,
    stack = TRUE
  )


  # prior
  addStyle(
    wb,
    sheet = sheet_name,
    style = prior_style,
    rows = data_rows,
    cols = 7,
    gridExpand = TRUE,
    stack = TRUE
  )


  # SNP.PP
  addStyle(
    wb,
    sheet = sheet_name,
    style = scientific_style,
    rows = data_rows,
    cols = 8,
    gridExpand = TRUE,
    stack = TRUE
  )


  # 90_credset
  addStyle(
    wb,
    sheet = sheet_name,
    style = integer_style,
    rows = data_rows,
    cols = 9,
    gridExpand = TRUE,
    stack = TRUE
  )


  # maf
  addStyle(
    wb,
    sheet = sheet_name,
    style = maf_style,
    rows = data_rows,
    cols = 10,
    gridExpand = TRUE,
    stack = TRUE
  )
}


# ---------------------------------------------------------
# Save Excel workbook
# ---------------------------------------------------------

saveWorkbook(
  wb,
  OUT_XLSX,
  overwrite = TRUE
)


# =========================================================
# 17. FINAL QUALITY CHECKS
# =========================================================

if (
  !all(
    excel_cs90$`90_credset` == 1L
  )
) {

  stop(
    "Quality-control error: the 90pct_credible_sets ",
    "worksheet contains rows with 90_credset != 1."
  )
}


if (
  !all(
    excel_pip$SNP.PP > 0.03
  )
) {

  stop(
    "Quality-control error: the PIP_above_0.03 ",
    "worksheet contains SNPs with PIP <= 0.03."
  )
}


if (
  !all(
    excel_cs90_pip$`90_credset` == 1L &
    excel_cs90_pip$SNP.PP > 0.03
  )
) {

  stop(
    "Quality-control error: the CS90_PIP_above_0.03 ",
    "worksheet contains rows that do not satisfy ",
    "both required conditions."
  )
}


# =========================================================
# 18. FINAL SUMMARY
# =========================================================

message("")
message("Done.")
message("")


message(
  "FUMA loci supplied: ",
  format(
    nrow(anm_loci),
    big.mark = ","
  )
)


message(
  "FUMA loci successfully fine-mapped: ",
  format(
    uniqueN(
      anm_out$locus
    ),
    big.mark = ","
  )
)


message(
  "Fine-mapped SNP-locus rows: ",
  format(
    nrow(anm_out),
    big.mark = ","
  )
)


message(
  "SNP-locus rows in 90% credible sets: ",
  format(
    nrow(anm_cs90),
    big.mark = ","
  )
)


message(
  "Top SNPs retained, one per fine-mapped locus: ",
  format(
    nrow(excel_top),
    big.mark = ","
  )
)


message(
  "SNP-locus rows with PIP > 0.03: ",
  format(
    nrow(excel_pip),
    big.mark = ","
  )
)


message(
  "CS90 SNP-locus rows with PIP > 0.03: ",
  format(
    nrow(excel_cs90_pip),
    big.mark = ","
  )
)


message(
  "Skipped loci: ",
  format(
    nrow(
      anm_pack$skipped
    ),
    big.mark = ","
  )
)


message(
  "Output directory: ",
  OUT_DIR
)


message(
  "Complete results CSV: ",
  OUT_ALL_CSV
)


message(
  "90% credible-set CSV: ",
  OUT_CS_CSV
)


message(
  "Locus-summary CSV: ",
  OUT_SUMMARY_CSV
)


message(
  "Clean Excel workbook: ",
  OUT_XLSX
)
