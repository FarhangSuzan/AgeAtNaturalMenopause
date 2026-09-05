#!/usr/bin/env python3

# ============================================================
# Age at Natural Menopause (ANM) Mendelian Randomization
# Sensitivity Analysis
# ============================================================
#
# Author: Suzan Farhang-Sardroodi
#
# This code was developed by Suzan Farhang-Sardroodi for the
# sensitivity Mendelian randomization analyses of age at natural
# menopause (ANM) reported in this study.
#
# The analysis applies a more stringent cis-eQTL instrument
# selection criterion while retaining the same MR framework used
# in the primary analysis.
#
# Please cite the associated publication when using or adapting
# this code.
#
# ============================================================
#!/usr/bin/env python3

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm


# ============================================================
# CONFIGURATION
# ============================================================

HOME = Path.home()
WORKDIR = HOME / "Desktop"

HARMONIZED_FILE = (
    WORKDIR /
    "ANM_GTEx_v8_MR_matched_harmonized.txt"
)

GENE_FILE = (
    WORKDIR /
    "ANM_gene_prioritization_FINAL_INTEGRATED_CLEAN.xlsx"
)

LD_DIR = (
    WORKDIR /
    "ANM_1000G_EUR_MR_PLINK"
)

PLINK2 = HOME / "bin" / "plink2"

PRIMARY_RESULTS = (
    WORKDIR /
    "ANM_GTEx_v8_MR_FINAL_LDaware.csv"
)

OUTPUT_INSTRUMENTS = (
    WORKDIR /
    "ANM_GTEx_v8_MR_SENSITIVITY_P5e8_instruments.txt"
)

OUTPUT_RESULTS = (
    WORKDIR /
    "ANM_GTEx_v8_MR_SENSITIVITY_P5e8_results.csv"
)

OUTPUT_SIGNIFICANT = (
    WORKDIR /
    "ANM_GTEx_v8_MR_SENSITIVITY_P5e8_SIGNIFICANT_FDR05.csv"
)

OUTPUT_GENE_SUMMARY = (
    WORKDIR /
    "ANM_GTEx_v8_MR_SENSITIVITY_P5e8_GENE_SUMMARY.csv"
)

OUTPUT_COMPARISON = (
    WORKDIR /
    "ANM_GTEx_v8_MR_SENSITIVITY_vs_PRIMARY.csv"
)


# ============================================================
# READ HARMONIZED GTEx–ANM DATA
# ============================================================

print("Reading harmonized GTEx-ANM data...")

df = pd.read_csv(
    HARMONIZED_FILE,
    sep="\t",
    low_memory=False
)

for col in [
    "pval_nominal",
    "pval_beta",
    "slope",
    "slope_se",
    "maf",
    "beta_outcome",
    "se_outcome"
]:
    df[col] = pd.to_numeric(
        df[col],
        errors="coerce"
    )

df["F_stat"] = (
    df["slope"] /
    df["slope_se"]
) ** 2


# ============================================================
# STRICT SENSITIVITY INSTRUMENT FILTER
# ============================================================

inst = df[
    (df["pval_beta"] < 0.05) &
    (df["pval_nominal"] < 5e-8) &
    (df["F_stat"] > 10) &
    (df["maf"] >= 0.01) &
    df["slope"].notna() &
    df["slope_se"].notna() &
    df["beta_outcome"].notna() &
    df["se_outcome"].notna()
].copy()

inst = inst.drop_duplicates(
    [
        "gene_clean",
        "tissue",
        "variant_id"
    ]
)

inst.to_csv(
    OUTPUT_INSTRUMENTS,
    sep="\t",
    index=False
)

pair_counts = (
    inst.groupby(
        ["gene_clean", "tissue"]
    )
    .size()
)

print("\n" + "=" * 70)
print("STRICT SENSITIVITY INSTRUMENT SET")
print("=" * 70)

print(f"Instrument rows:          {len(inst):,}")
print(f"Unique genes:             {inst['gene_clean'].nunique():,}")
print(
    "Gene-tissue pairs:        "
    f"{inst[['gene_clean','tissue']].drop_duplicates().shape[0]:,}"
)
print(f"Unique variants:          {inst['variant_id'].nunique():,}")
print(f"Pairs with 1 instrument:  {(pair_counts == 1).sum():,}")
print(f"Pairs with >=2:           {(pair_counts >= 2).sum():,}")


# ============================================================
# STANDARDIZE VARIABLES
# ============================================================

inst["chromosome"] = (
    inst["chromosome"]
    .astype(str)
    .str.replace("chr", "", regex=False)
)

inst["base_pair_location"] = pd.to_numeric(
    inst["base_pair_location"],
    errors="coerce"
).astype("Int64")

# GTEx slope is relative to ALT allele.
inst["exposure_allele"] = (
    inst["alt"]
    .astype(str)
    .str.upper()
)


# ============================================================
# ADD GENE SYMBOLS
# ============================================================

genes = pd.read_excel(
    GENE_FILE,
    sheet_name="Gene_Prioritization"
)

genes["gene_clean"] = (
    genes["Ensembl_ID"]
    .astype(str)
    .str.split(".").str[0]
)

gene_map = (
    genes[
        ["gene_clean", "Gene"]
    ]
    .drop_duplicates("gene_clean")
)

inst = inst.merge(
    gene_map,
    on="gene_clean",
    how="left"
)


# ============================================================
# MAP TO DENSE 1000G EUR LD REFERENCE
# ============================================================

mapped = []

for chrom in range(1, 23):

    g = inst[
        inst["chromosome"] == str(chrom)
    ].copy()

    if g.empty:
        continue

    pvar = LD_DIR / f"1000G_EUR_chr{chrom}.pvar"

    pv = pd.read_csv(
        pvar,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "ld_chr",
            "ld_pos",
            "ld_id",
            "ld_ref",
            "ld_alt"
        ],
        usecols=[0, 1, 2, 3, 4]
    )

    pv["ld_pos"] = pd.to_numeric(
        pv["ld_pos"],
        errors="coerce"
    ).astype("Int64")

    pv["ld_ref"] = (
        pv["ld_ref"]
        .astype(str)
        .str.upper()
    )

    pv["ld_alt"] = (
        pv["ld_alt"]
        .astype(str)
        .str.upper()
    )

    pv = pv[
        ~pv["ld_alt"].str.contains(
            ",",
            regex=False
        )
    ].copy()

    m = g.merge(
        pv,
        left_on="base_pair_location",
        right_on="ld_pos",
        how="left"
    )

    ea = (
        m["effect_allele"]
        .astype(str)
        .str.upper()
    )

    oa = (
        m["other_allele"]
        .astype(str)
        .str.upper()
    )

    allele_match = (
        (
            (ea == m["ld_ref"]) &
            (oa == m["ld_alt"])
        )
        |
        (
            (ea == m["ld_alt"]) &
            (oa == m["ld_ref"])
        )
    )

    m.loc[
        ~allele_match,
        ["ld_id", "ld_ref", "ld_alt"]
    ] = pd.NA

    m["has_ld"] = m["ld_id"].notna()

    m = (
        m.sort_values(
            "has_ld",
            ascending=False
        )
        .drop_duplicates(
            [
                "gene_clean",
                "tissue",
                "variant_id"
            ],
            keep="first"
        )
    )

    mapped.append(m)

inst = pd.concat(
    mapped,
    ignore_index=True
)

print(
    "\nInstrument rows mapped to dense 1000G EUR: "
    f"{inst['ld_id'].notna().sum():,} / {len(inst):,}"
)


# ============================================================
# RUN MR
# ============================================================

groups = list(
    inst.groupby(
        ["gene_clean", "tissue"],
        sort=False
    )
)

temp_dir = Path(
    tempfile.mkdtemp(
        prefix="ANM_SENS_LD_",
        dir=str(WORKDIR)
    )
)

results = []

print(f"\nRunning {len(groups):,} sensitivity MR tests...\n")


for i, ((gene, tissue), g) in enumerate(
    groups,
    start=1
):

    g = (
        g.sort_values("pval_nominal")
        .drop_duplicates("variant_id")
        .copy()
    )

    symbol = (
        g["Gene"].dropna().iloc[0]
        if g["Gene"].notna().any()
        else gene
    )

    n_original = len(g)

    # --------------------------------------------------------
    # ONE SNP → WALD RATIO
    # --------------------------------------------------------

    if n_original == 1:

        row = g.iloc[0]

        bx = float(row["slope"])
        by = float(row["beta_outcome"])
        sey = float(row["se_outcome"])

        beta_mr = by / bx
        se_mr = abs(sey / bx)

        method = "Wald ratio"
        n_used = 1

        max_abs_r = np.nan
        condition_number = np.nan

    else:

        g_ld = g[
            g["ld_id"].notna()
        ].copy()

        # ----------------------------------------------------
        # FALLBACK
        # ----------------------------------------------------

        if len(g_ld) < 2:

            row = (
                g_ld.iloc[0]
                if len(g_ld) == 1
                else g.iloc[0]
            )

            bx = float(row["slope"])
            by = float(row["beta_outcome"])
            sey = float(row["se_outcome"])

            beta_mr = by / bx
            se_mr = abs(sey / bx)

            method = "Wald fallback: 1 LD-mapped SNP"
            n_used = 1

            max_abs_r = np.nan
            condition_number = np.nan

        else:

            chrom = str(
                g_ld["chromosome"].iloc[0]
            )

            pfile = (
                LD_DIR /
                f"1000G_EUR_chr{chrom}"
            )

            ids_file = (
                temp_dir /
                f"ids_{i}.txt"
            )

            out_prefix = (
                temp_dir /
                f"ld_{i}"
            )

            g_ld["ld_id"].to_csv(
                ids_file,
                index=False,
                header=False
            )

            cmd = [
                str(PLINK2),
                "--pfile",
                str(pfile),
                "--extract",
                str(ids_file),
                "--r-unphased",
                "square",
                "ref-based",
                "--out",
                str(out_prefix)
            ]

            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True
            )

            matrix_file = Path(
                str(out_prefix) +
                ".unphased.vcor1"
            )

            vars_file = Path(
                str(out_prefix) +
                ".unphased.vcor1.vars"
            )

            if (
                proc.returncode != 0
                or not matrix_file.exists()
                or not vars_file.exists()
            ):

                row = g.iloc[0]

                bx = float(row["slope"])
                by = float(row["beta_outcome"])
                sey = float(row["se_outcome"])

                beta_mr = by / bx
                se_mr = abs(sey / bx)

                method = "Wald fallback: LD calculation failed"
                n_used = 1

                max_abs_r = np.nan
                condition_number = np.nan

            else:

                ld_ids = pd.read_csv(
                    vars_file,
                    header=None
                )[0].astype(str).tolist()

                rho_ref = np.atleast_2d(
                    np.loadtxt(
                        matrix_file
                    )
                )

                g_indexed = (
                    g_ld
                    .drop_duplicates("ld_id")
                    .set_index("ld_id")
                )

                valid_ids = [
                    snp for snp in ld_ids
                    if snp in g_indexed.index
                ]

                g2 = (
                    g_indexed
                    .loc[valid_ids]
                    .reset_index()
                )

                positions = {
                    snp: j
                    for j, snp
                    in enumerate(ld_ids)
                }

                order = [
                    positions[snp]
                    for snp in valid_ids
                ]

                rho_ref = rho_ref[
                    np.ix_(
                        order,
                        order
                    )
                ]

                # --------------------------------------------
                # Orient LD to GTEx exposure allele
                # --------------------------------------------

                signs = []

                for _, row in g2.iterrows():

                    exposure = str(
                        row["exposure_allele"]
                    ).upper()

                    ref = str(
                        row["ld_ref"]
                    ).upper()

                    alt = str(
                        row["ld_alt"]
                    ).upper()

                    if exposure == ref:
                        signs.append(1.0)

                    elif exposure == alt:
                        signs.append(-1.0)

                    else:
                        signs.append(np.nan)

                signs = np.asarray(signs)

                valid = np.isfinite(signs)

                g2 = (
                    g2.loc[valid]
                    .reset_index(drop=True)
                )

                rho_ref = rho_ref[
                    np.ix_(valid, valid)
                ]

                signs = signs[valid]

                if len(g2) < 2:

                    row = g2.iloc[0]

                    bx = float(row["slope"])
                    by = float(row["beta_outcome"])
                    sey = float(row["se_outcome"])

                    beta_mr = by / bx
                    se_mr = abs(sey / bx)

                    method = "Wald fallback: 1 oriented SNP"
                    n_used = 1

                    max_abs_r = np.nan
                    condition_number = np.nan

                else:

                    rho = (
                        np.outer(
                            signs,
                            signs
                        )
                        * rho_ref
                    )

                    bx = (
                        g2["slope"]
                        .to_numpy(float)
                    )

                    by = (
                        g2["beta_outcome"]
                        .to_numpy(float)
                    )

                    sey = (
                        g2["se_outcome"]
                        .to_numpy(float)
                    )

                    omega = (
                        np.outer(
                            sey,
                            sey
                        )
                        * rho
                    )

                    condition_number = (
                        np.linalg.cond(
                            omega
                        )
                    )

                    omega_inv = (
                        np.linalg.pinv(
                            omega,
                            rcond=1e-10
                        )
                    )

                    denominator = (
                        bx.T
                        @ omega_inv
                        @ bx
                    )

                    if (
                        not np.isfinite(denominator)
                        or denominator <= 0
                    ):

                        row = g.iloc[0]

                        bx0 = float(row["slope"])
                        by0 = float(row["beta_outcome"])
                        sy0 = float(row["se_outcome"])

                        beta_mr = by0 / bx0
                        se_mr = abs(sy0 / bx0)

                        method = "Wald fallback: unstable LD"
                        n_used = 1

                        max_abs_r = np.nan

                    else:

                        beta_mr = (
                            bx.T
                            @ omega_inv
                            @ by
                        ) / denominator

                        se_mr = np.sqrt(
                            1.0 /
                            denominator
                        )

                        method = (
                            "LD-aware generalized IVW"
                        )

                        n_used = len(g2)

                        offdiag = rho[
                            ~np.eye(
                                len(rho),
                                dtype=bool
                            )
                        ]

                        max_abs_r = (
                            np.max(
                                np.abs(offdiag)
                            )
                            if len(offdiag)
                            else 0.0
                        )

    # --------------------------------------------------------
    # MR statistics
    # --------------------------------------------------------

    z = beta_mr / se_mr
    p = 2 * norm.sf(abs(z))

    results.append({
        "gene": gene,
        "symbol": symbol,
        "tissue": tissue,
        "N_instruments_original": n_original,
        "N_instruments_used": n_used,
        "method": method,
        "MR_beta": beta_mr,
        "MR_SE": se_mr,
        "MR_Z": z,
        "MR_P": p,
        "CI_lower": beta_mr - 1.96 * se_mr,
        "CI_upper": beta_mr + 1.96 * se_mr,
        "max_abs_LD_r": max_abs_r,
        "LD_cov_condition_number": condition_number
    })

    print(
        f"\r{i:,}/{len(groups):,}",
        end=""
    )


print("\n")


# ============================================================
# FDR CORRECTION
# ============================================================

res = pd.DataFrame(results)

pvals = res["MR_P"].to_numpy()

order = np.argsort(pvals)
ranked = pvals[order]

m = len(pvals)

adjusted = (
    ranked
    * m
    / np.arange(1, m + 1)
)

adjusted = np.minimum.accumulate(
    adjusted[::-1]
)[::-1]

adjusted = np.minimum(
    adjusted,
    1
)

fdr = np.empty(m)
fdr[order] = adjusted

res["MR_FDR"] = fdr

res = res.sort_values(
    ["MR_FDR", "MR_P"]
)

res.to_csv(
    OUTPUT_RESULTS,
    index=False
)


# ============================================================
# SIGNIFICANT SENSITIVITY RESULTS
# ============================================================

sig = res[
    res["MR_FDR"] < 0.05
].copy()

sig.to_csv(
    OUTPUT_SIGNIFICANT,
    index=False
)


# ============================================================
# GENE-LEVEL SENSITIVITY SUMMARY
# ============================================================

gene_summary = (
    sig.groupby(
        ["gene", "symbol"]
    )
    .agg(
        N_significant_tissues=(
            "tissue",
            "nunique"
        ),
        significant_tissues=(
            "tissue",
            lambda x:
            ";".join(
                sorted(
                    set(
                        x.astype(str)
                    )
                )
            )
        ),
        Best_MR_P=(
            "MR_P",
            "min"
        ),
        Best_MR_FDR=(
            "MR_FDR",
            "min"
        )
    )
    .reset_index()
)

gene_summary.to_csv(
    OUTPUT_GENE_SUMMARY,
    index=False
)


# ============================================================
# COMPARE WITH PRIMARY MR
# ============================================================

primary = pd.read_csv(
    PRIMARY_RESULTS
)

primary_sig = primary[
    primary["MR_FDR"] < 0.05
].copy()

primary_sig["primary_significant"] = True
sig["sensitivity_significant"] = True

comparison = primary_sig[
    [
        "gene",
        "symbol",
        "tissue",
        "MR_beta",
        "MR_P",
        "MR_FDR",
        "primary_significant"
    ]
].rename(
    columns={
        "MR_beta": "Primary_MR_beta",
        "MR_P": "Primary_MR_P",
        "MR_FDR": "Primary_MR_FDR"
    }
)

sensitivity_compare = sig[
    [
        "gene",
        "symbol",
        "tissue",
        "MR_beta",
        "MR_P",
        "MR_FDR",
        "sensitivity_significant"
    ]
].rename(
    columns={
        "MR_beta":
            "Sensitivity_MR_beta",
        "MR_P":
            "Sensitivity_MR_P",
        "MR_FDR":
            "Sensitivity_MR_FDR"
    }
)

comparison = comparison.merge(
    sensitivity_compare,
    on=[
        "gene",
        "symbol",
        "tissue"
    ],
    how="outer"
)

comparison[
    "primary_significant"
] = comparison[
    "primary_significant"
].fillna(False)

comparison[
    "sensitivity_significant"
] = comparison[
    "sensitivity_significant"
].fillna(False)

comparison.to_csv(
    OUTPUT_COMPARISON,
    index=False
)


# ============================================================
# FINAL SUMMARY
# ============================================================

primary_genes = set(
    primary_sig["gene"]
)

sens_genes = set(
    sig["gene"]
)

overlap_genes = (
    primary_genes &
    sens_genes
)

primary_pairs = set(
    zip(
        primary_sig["gene"],
        primary_sig["tissue"]
    )
)

sens_pairs = set(
    zip(
        sig["gene"],
        sig["tissue"]
    )
)

overlap_pairs = (
    primary_pairs &
    sens_pairs
)


print("=" * 72)
print("STRICT MR SENSITIVITY ANALYSIS COMPLETE")
print("=" * 72)

print(f"Sensitivity MR tests:             {len(res):,}")
print(f"Genes tested:                     {res['gene'].nunique():,}")

print(f"\nSensitivity FDR < 0.05 rows:      {len(sig):,}")
print(f"Sensitivity significant genes:    {sig['gene'].nunique():,}")

print(f"\nPrimary significant rows:         {len(primary_sig):,}")
print(f"Primary significant genes:        {primary_sig['gene'].nunique():,}")

print(f"\nSignificant gene-tissue overlap:   {len(overlap_pairs):,}")
print(f"Significant gene overlap:          {len(overlap_genes):,}")

if len(primary_genes):
    print(
        "Primary genes retained in sensitivity: "
        f"{100 * len(overlap_genes) / len(primary_genes):.2f}%"
    )

print("\nSaved:")
print(OUTPUT_INSTRUMENTS)
print(OUTPUT_RESULTS)
print(OUTPUT_SIGNIFICANT)
print(OUTPUT_GENE_SUMMARY)
print(OUTPUT_COMPARISON)

shutil.rmtree(
    temp_dir,
    ignore_errors=True
)

