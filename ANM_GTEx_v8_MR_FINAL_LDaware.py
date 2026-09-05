#!/usr/bin/env python3

# ============================================================
# Age at Natural Menopause (ANM) Mendelian Randomization
# ============================================================
#
# Author: Suzan Farhang-Sardroodi
#
# This code was developed by Suzan Farhang-Sardroodi for the
# Mendelian randomization analyses of age at natural menopause
# (ANM) reported in this study.
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

INSTRUMENT_FILE = (
    WORKDIR /
    "ANM_GTEx_v8_MR_final_instruments.txt"
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

OUTPUT_RESULTS = (
    WORKDIR /
    "ANM_GTEx_v8_MR_FINAL_LDaware.csv"
)

OUTPUT_QC = (
    WORKDIR /
    "ANM_GTEx_v8_MR_FINAL_LDaware_QC.csv"
)


# ============================================================
# READ INSTRUMENT DATA
# ============================================================

df = pd.read_csv(
    INSTRUMENT_FILE,
    sep="\t",
    low_memory=False
)

numeric_columns = [
    "slope",
    "slope_se",
    "beta_outcome",
    "se_outcome",
    "pval_nominal"
]

for col in numeric_columns:
    df[col] = pd.to_numeric(
        df[col],
        errors="coerce"
    )

df["chromosome"] = (
    df["chromosome"]
    .astype(str)
    .str.replace("chr", "", regex=False)
)

df["base_pair_location"] = pd.to_numeric(
    df["base_pair_location"],
    errors="coerce"
).astype("Int64")

# GTEx v8 slope is aligned to the GTEx ALT allele.
df["exposure_allele"] = (
    df["alt"]
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

df = df.merge(
    gene_map,
    on="gene_clean",
    how="left"
)


# ============================================================
# MAP INSTRUMENTS TO DENSE 1000G EUR REFERENCE
# ============================================================

mapped = []

for chrom in range(1, 23):

    g = df[
        df["chromosome"] == str(chrom)
    ].copy()

    if g.empty:
        continue

    pvar = (
        LD_DIR /
        f"1000G_EUR_chr{chrom}.pvar"
    )

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

    # Restrict to biallelic variants.
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

df = pd.concat(
    mapped,
    ignore_index=True
)


# ============================================================
# MR ANALYSIS
# ============================================================

groups = list(
    df.groupby(
        ["gene_clean", "tissue"],
        sort=False
    )
)

temp_dir = Path(
    tempfile.mkdtemp(
        prefix="ANM_MR_LD_",
        dir=str(WORKDIR)
    )
)

results = []
qc = []


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

    # ========================================================
    # SINGLE-INSTRUMENT WALD RATIO
    # ========================================================

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


    # ========================================================
    # MULTI-INSTRUMENT ANALYSIS
    # ========================================================

    else:

        g_ld = (
            g[g["ld_id"].notna()]
            .copy()
        )

        # ----------------------------------------------------
        # FALLBACK IF <2 SNPs ARE PRESENT IN LD REFERENCE
        # ----------------------------------------------------

        if len(g_ld) < 2:

            if len(g_ld) == 1:
                row = g_ld.iloc[0]
            else:
                row = g.iloc[0]

            bx = float(row["slope"])
            by = float(row["beta_outcome"])
            sey = float(row["se_outcome"])

            beta_mr = by / bx
            se_mr = abs(sey / bx)

            method = "Wald fallback: 1 LD-mapped SNP"
            n_used = 1
            max_abs_r = np.nan
            condition_number = np.nan


        # ----------------------------------------------------
        # LD-AWARE GENERALIZED IVW
        # ----------------------------------------------------

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

                method = (
                    "Wald fallback: "
                    "LD calculation failed"
                )

                n_used = 1
                max_abs_r = np.nan
                condition_number = np.nan

            else:

                ld_ids = pd.read_csv(
                    vars_file,
                    header=None
                )[0].astype(str).tolist()

                rho_ref = np.loadtxt(
                    matrix_file
                )

                rho_ref = np.atleast_2d(
                    rho_ref
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

                # ============================================
                # ORIENT LD TO GTEx EXPOSURE ALLELE
                # ============================================

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
                    sey = float(
                        row["se_outcome"]
                    )

                    beta_mr = by / bx
                    se_mr = abs(sey / bx)

                    method = (
                        "Wald fallback: "
                        "1 oriented SNP"
                    )

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
                        sy0 = float(
                            row["se_outcome"]
                        )

                        beta_mr = by0 / bx0
                        se_mr = abs(sy0 / bx0)

                        method = (
                            "Wald fallback: "
                            "unstable LD"
                        )

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


    # ========================================================
    # MR TEST STATISTICS
    # ========================================================

    z = beta_mr / se_mr

    p = 2 * norm.sf(
        abs(z)
    )

    results.append({
        "gene": gene,
        "symbol": symbol,
        "tissue": tissue,
        "N_instruments_original":
            n_original,
        "N_instruments_used":
            n_used,
        "method":
            method,
        "MR_beta":
            beta_mr,
        "MR_SE":
            se_mr,
        "MR_Z":
            z,
        "MR_P":
            p,
        "CI_lower":
            beta_mr - 1.96 * se_mr,
        "CI_upper":
            beta_mr + 1.96 * se_mr,
        "max_abs_LD_r":
            max_abs_r,
        "LD_cov_condition_number":
            condition_number
    })

    qc.append({
        "gene": gene,
        "symbol": symbol,
        "tissue": tissue,
        "N_original": n_original,
        "N_used": n_used,
        "method": method
    })


# ============================================================
# BENJAMINI-HOCHBERG FDR
# ============================================================

res = pd.DataFrame(results)

pvals = res["MR_P"].to_numpy()

order = np.argsort(pvals)
ranked = pvals[order]

m = len(pvals)

adjusted = (
    ranked
    * m
    / np.arange(
        1,
        m + 1
    )
)

adjusted = (
    np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]
)

adjusted = np.minimum(
    adjusted,
    1
)

fdr = np.empty(m)
fdr[order] = adjusted

res["MR_FDR"] = fdr

res = res.sort_values(
    "MR_P"
)

res.to_csv(
    OUTPUT_RESULTS,
    index=False
)

pd.DataFrame(qc).to_csv(
    OUTPUT_QC,
    index=False
)

shutil.rmtree(
    temp_dir,
    ignore_errors=True
)

print("=" * 70)
print("FINAL LD-AWARE MR COMPLETE")
print("=" * 70)

print(
    f"MR tests: {len(res):,}"
)

print(
    f"Unique genes: "
    f"{res['gene'].nunique():,}"
)

print(
    f"FDR < 0.05: "
    f"{(res['MR_FDR'] < 0.05).sum():,}"
)

print("\nSaved:")
print(OUTPUT_RESULTS)
print(OUTPUT_QC)

