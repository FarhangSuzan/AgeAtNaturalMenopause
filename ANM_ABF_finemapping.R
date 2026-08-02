suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

# ---------------------------------------------------------
# Install/load openxlsx for Excel output
# ---------------------------------------------------------

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

ANM_FILE <- file.path(
  path.expand("~/Desktop"),
  "ANM_GCST90292553.txt"
)

LEAD_FILE <- file.path(
  path.expand("~/Desktop"),
  "ANM_chromosome_plots",
  "ANM_lead_SNPs_500kb_windows.csv"
)

OUT_DIR <- file.path(
  path.expand("~/Desktop"),
  "ANM_ABF_finemapping"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Locus boundaries: ±500 kb from each lead SNP
FLANK <- 500000L

# Credible-set threshold
CRED_P <- 0.90

# Prior probability that any one SNP is causal
PRIOR_P <- 1e-4

# Minimum number of valid SNPs needed to analyze a locus
MIN_SNPS <- 10L

OUT_XLSX <- file.path(
  OUT_DIR,
  "ANM_ABF_finemapping_90pct.xlsx"
)

OUT_ALL_CSV <- file.path(
  OUT_DIR,
  "ANM_ABF_all_loci.csv"
)

OUT_CS_CSV <- file.path(
  OUT_DIR,
  "ANM_ABF_90pct_credible_sets.csv"
)

OUT_SUMMARY_CSV <- file.path(
  OUT_DIR,
  "ANM_ABF_locus_summary.csv"
)


# =========================================================
# 2. READ AND STANDARDIZE ANM GWAS
# =========================================================

read_anm <- function(path) {

  if (!file.exists(path)) {
    stop("ANM file not found: ", path)
  }

  message("Reading ANM GWAS file...")
  dt <- fread(path)

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
      paste(missing_columns, collapse = ", ")
    )
  }

  dt <- dt[
    ,
    .(
      snp = as.character(rs_id),
      chr = as.character(chromosome),
      position = as.integer(base_pair_location),
      a1 = toupper(as.character(effect_allele)),
      a2 = toupper(as.character(other_allele)),
      beta = as.numeric(beta),
      se = as.numeric(standard_error),
      eaf = as.numeric(effect_allele_frequency),
      p = as.numeric(p_value)
    )
  ]

  # Remove "chr" prefix if present
  dt[
    ,
    chr := sub(
      "^chr",
      "",
      chr,
      ignore.case = TRUE
    )
  ]

  # Standardize chromosome X to 23
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

  # Convert EAF to MAF
  dt[
    ,
    maf := pmin(
      eaf,
      1 - eaf
    )
  ]

  # Remove invalid rows
  dt <- dt[
    !is.na(snp) &
    snp != "" &
    !is.na(chr_num) &
    chr_num >= 1 &
    chr_num <= 23 &
    !is.na(position) &
    position > 0 &
    !is.na(beta) &
    !is.na(se) &
    se > 0
  ]

  # Remove duplicated SNP IDs.
  # If duplicates exist, keep the row with the smallest SE.
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
      " duplicated SNP IDs found; keeping the row ",
      "with the smallest standard error."
    )

    dt <- dt[
      !duplicated(snp)
    ]
  }

  # Create Z score
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
    format(nrow(dt), big.mark = ",")
  )

  dt
}


# =========================================================
# 3. READ LEAD-SNP TABLE
# =========================================================

read_anm_leads <- function(path) {

  if (!file.exists(path)) {
    stop("Lead-SNP file not found: ", path)
  }

  message("Reading ANM lead-SNP table...")
  leads <- fread(path)

  # Accept the exact names from our Python output
  if (
    all(
      c("CHR", "POS", "SNPID") %in%
      names(leads)
    )
  ) {

    leads <- leads[
      ,
      .(
        chr_num = as.integer(CHR),
        position = as.integer(POS),
        lead_snp = as.character(SNPID)
      )
    ]

  } else if (
    all(
      c(
        "Chromosome_label",
        "POS",
        "SNPID"
      ) %in% names(leads)
    )
  ) {

    leads <- leads[
      ,
      .(
        chr_num = fifelse(
          toupper(
            as.character(
              Chromosome_label
            )
          ) == "X",
          23L,
          as.integer(
            Chromosome_label
          )
        ),
        position = as.integer(POS),
        lead_snp = as.character(SNPID)
      )
    ]

  } else {

    stop(
      "Lead file must contain CHR, POS, SNPID, ",
      "or Chromosome_label, POS, SNPID."
    )
  }

  leads <- unique(
    leads[
      !is.na(chr_num) &
      chr_num >= 1 &
      chr_num <= 23 &
      !is.na(position) &
      position > 0 &
      !is.na(lead_snp) &
      lead_snp != ""
    ]
  )

  setorder(
    leads,
    chr_num,
    position
  )

  # Create readable locus IDs
  leads[
    ,
    chr_label := fifelse(
      chr_num == 23,
      "X",
      as.character(chr_num)
    )
  ]

  leads[
    ,
    locus := sprintf(
      "chr%s_%s_%d",
      chr_label,
      lead_snp,
      position
    )
  ]

  message(
    "Lead SNPs retained: ",
    nrow(leads)
  )

  leads
}


# =========================================================
# 4. BUILD ONE ABF DATASET PER LOCUS
# =========================================================

build_anm_inputs <- function(
  gwas,
  leads,
  flank = 500000L,
  min_snps = 10L
) {

  inputs <- list()
  metadata <- list()
  skipped <- list()

  for (i in seq_len(nrow(leads))) {

    lead <- leads[i]

    locus_name <- lead$locus
    lead_snp <- lead$lead_snp
    lead_chr <- lead$chr_num
    lead_pos <- lead$position

    # Prefer the lead SNP's position in the GWAS file
    lead_match <- gwas[
      snp == lead_snp
    ]

    if (nrow(lead_match) > 0) {

      lead_match <- lead_match[1]

      if (
        lead_match$chr_num != lead_chr ||
        lead_match$position != lead_pos
      ) {
        warning(
          "Lead-file coordinate differs from GWAS for ",
          lead_snp,
          ". Using the GWAS coordinate."
        )
      }

      lead_chr <- lead_match$chr_num
      lead_pos <- lead_match$position

    } else {

      warning(
        "Lead SNP not found by rsID in the ANM GWAS: ",
        lead_snp,
        ". Using coordinates from lead table."
      )
    }

    region_start <- max(
      1L,
      lead_pos - flank
    )

    region_end <- lead_pos + flank

    locus_snps <- gwas[
      chr_num == lead_chr &
      position >= region_start &
      position <= region_end
    ]

    # Ensure required values are finite
    locus_snps <- locus_snps[
      is.finite(beta) &
      is.finite(se) &
      se > 0 &
      !is.na(snp) &
      snp != ""
    ]

    if (nrow(locus_snps) < min_snps) {

      warning(
        "Skipping ",
        locus_name,
        ": only ",
        nrow(locus_snps),
        " valid SNPs."
      )

      skipped[[locus_name]] <- data.table(
        locus = locus_name,
        lead_snp = lead_snp,
        chr = lead_chr,
        lead_position = lead_pos,
        reason = paste0(
          "Fewer than ",
          min_snps,
          " valid SNPs"
        )
      )

      next
    }

    # Beta and varbeta are sufficient for ABF fine-mapping.
    # ANM is a quantitative phenotype.
    inp <- list(
  beta = locus_snps$beta,
  varbeta = locus_snps$se^2,
  MAF = locus_snps$maf,
  N = 243506,
  snp = locus_snps$snp,
  position = locus_snps$position,
  type = "quant"
)

    inputs[[locus_name]] <- inp

    metadata[[locus_name]] <- data.table(
      locus = locus_name,
      lead_snp = lead_snp,
      chr = lead_chr,
      chr_label = ifelse(
        lead_chr == 23,
        "X",
        as.character(lead_chr)
      ),
      lead_position = lead_pos,
      region_start = region_start,
      region_end = region_end,
      number_of_snps = nrow(locus_snps)
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
# 5. RUN finemap.abf()
# =========================================================

run_abf_finemapping <- function(
  input_list,
  prior = 1e-4
) {

  results <- vector(
    "list",
    length(input_list)
  )

  names(results) <- names(input_list)

  for (locus_name in names(input_list)) {

    message(
      "Fine-mapping ",
      locus_name,
      "..."
    )

    result <- tryCatch(
      {
        finemap.abf(
          dataset = input_list[[locus_name]],
          p1 = prior
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

    results[[locus_name]] <- result
  }

  results[
    !vapply(
      results,
      is.null,
      logical(1)
    )
  ]
}


# =========================================================
# 6. CONVERT RESULTS TO ONE TABLE
# =========================================================

finemap_to_dt <- function(fm_list) {

  if (length(fm_list) == 0) {
    return(data.table())
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
# 7. CREATE 90% CREDIBLE SETS
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
      "No posterior-probability column found."
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

      setorderv(
        x,
        posterior_column,
        order = -1L
      )

      x[
        ,
        cumulative_PP :=
          cumsum(
            get(posterior_column)
          )
      ]

      # Include the first SNP that makes cumulative PP
      # reach or exceed the target threshold.
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
# 8. ADD ORIGINAL ANM INFORMATION
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
        chr = chr_num,
        position_gwas = position,
        a1,
        a2,
        eaf,
        maf,
        beta_original = beta,
        se_original = se,
        p_original = p
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
# 9. CREATE LOCUS-LEVEL SUMMARY
# =========================================================

create_locus_summary <- function(
  fm,
  metadata
) {

  summary <- fm[
    snp != "null",
    {
      ordered <- .SD[
        order(-SNP.PP)
      ]

      cs <- ordered[
        CS90 == 1L
      ]

      .(
        top_posterior_snp =
          ordered$snp[1],

        top_SNP_PP =
          ordered$SNP.PP[1],

        number_of_finemapped_snps =
          .N,

        CS90_size =
          nrow(cs),

        CS90_total_PP =
          sum(
            cs$SNP.PP,
            na.rm = TRUE
          ),

        CS90_min_PP =
          if (
            nrow(cs) > 0
          ) {
            min(
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
# 10. RUN THE COMPLETE ANALYSIS
# =========================================================

anm <- read_anm(
  ANM_FILE
)

anm_leads <- read_anm_leads(
  LEAD_FILE
)

message(
  "Building ±",
  FLANK / 1000,
  "-kb loci..."
)

anm_pack <- build_anm_inputs(
  gwas = anm,
  leads = anm_leads,
  flank = FLANK,
  min_snps = MIN_SNPS
)

message(
  "Loci prepared for fine-mapping: ",
  length(anm_pack$inputs)
)

message(
  "Running coloc::finemap.abf()..."
)

anm_fm_list <- run_abf_finemapping(
  input_list = anm_pack$inputs,
  prior = PRIOR_P
)

anm_raw <- finemap_to_dt(
  anm_fm_list
)

if (nrow(anm_raw) == 0) {
  stop(
    "No fine-mapping results were produced."
  )
}

# Remove the null model before creating SNP credible sets
anm_raw <- anm_raw[
  snp != "null"
]

anm_out <- add_credible_set(
  anm_raw,
  credible_probability = CRED_P
)

anm_out <- merge_anm_information(
  anm_out,
  anm
)


# =========================================================
# 11. COLUMN ORDER
# =========================================================

preferred_columns <- c(
  "snp",
  "locus",
  "chr",
  "position",
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
  "CS90"
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
# 12. CREATE CREDIBLE-SET AND SUMMARY TABLES
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
  lead_position
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


# =========================================================
# 13. WRITE CSV OUTPUTS
# =========================================================

fwrite(
  anm_out,
  OUT_ALL_CSV
)

fwrite(
  anm_cs90,
  OUT_CS_CSV
)

fwrite(
  anm_summary,
  OUT_SUMMARY_CSV
)

if (nrow(anm_pack$skipped) > 0) {

  fwrite(
    anm_pack$skipped,
    file.path(
      OUT_DIR,
      "ANM_ABF_skipped_loci.csv"
    )
  )
}


# =========================================================
# 14. WRITE EXCEL WORKBOOK
# =========================================================

message(
  "Writing Excel workbook: ",
  OUT_XLSX
)

wb <- createWorkbook()

addWorksheet(
  wb,
  "Locus_summary"
)

addWorksheet(
  wb,
  "CS90"
)

addWorksheet(
  wb,
  "All_SNPs"
)

addWorksheet(
  wb,
  "Lead_SNPs"
)

writeDataTable(
  wb,
  "Locus_summary",
  anm_summary
)

writeDataTable(
  wb,
  "CS90",
  anm_cs90
)

writeDataTable(
  wb,
  "All_SNPs",
  anm_out
)

writeDataTable(
  wb,
  "Lead_SNPs",
  anm_leads
)

if (nrow(anm_pack$skipped) > 0) {

  addWorksheet(
    wb,
    "Skipped_loci"
  )

  writeDataTable(
    wb,
    "Skipped_loci",
    anm_pack$skipped
  )
}

freezePane(
  wb,
  "Locus_summary",
  firstRow = TRUE
)

freezePane(
  wb,
  "CS90",
  firstRow = TRUE
)

freezePane(
  wb,
  "All_SNPs",
  firstRow = TRUE
)

setColWidths(
  wb,
  sheet = "Locus_summary",
  cols = 1:ncol(anm_summary),
  widths = "auto"
)

setColWidths(
  wb,
  sheet = "CS90",
  cols = 1:ncol(anm_cs90),
  widths = "auto"
)

saveWorkbook(
  wb,
  OUT_XLSX,
  overwrite = TRUE
)


# =========================================================
# 15. FINAL SUMMARY
# =========================================================

message("Done.")

message(
  "Fine-mapped SNP rows: ",
  format(
    nrow(anm_out),
    big.mark = ","
  )
)

message(
  "Fine-mapped loci: ",
  uniqueN(anm_out$locus)
)

message(
  "SNPs included in 90% credible sets: ",
  format(
    nrow(anm_cs90),
    big.mark = ","
  )
)

message(
  "Output directory: ",
  OUT_DIR
)

message(
  "Excel workbook: ",
  OUT_XLSX
)
