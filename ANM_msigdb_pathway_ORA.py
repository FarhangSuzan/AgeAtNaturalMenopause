# ============================================================
# FINAL ANM MSigDB OVER-REPRESENTATION ANALYSIS
#
# Collections:
#   Reactome
#   GO Biological Process
#   GO Cellular Component
#   Hallmark
#
# Statistics:
#   - Collection-specific MSigDB gene universe
#   - One-sided Fisher's exact test
#   - Benjamini-Hochberg FDR correction
#
# Robust feature:
#   Automatically detects pathway / GO-term column names.
# ============================================================

import pandas as pd
import numpy as np
import re
from pathlib import Path

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


# ============================================================
# 1. FILE PATHS
# ============================================================

desktop = Path.home() / "Desktop"
msigdb_folder = desktop / "msigdb,Useful"

anm_file = (
    desktop
    / "ANM_1289_UNIQUE_GENES_FOR_PATHWAY_ANALYSIS.xlsx"
)

output_file = (
    desktop
    / "ANM_MSIGDB_ENRICHMENT_ALL_COLLECTIONS_FINAL.xlsx"
)

collection_files = {

    "Reactome":
        msigdb_folder
        / "c2.cp.reactome.v2026.1.Hs.symbols.xlsx",

    "GO_BP":
        msigdb_folder
        / "c5.go.bp.v2026.1.Hs.symbols.xlsx",

    "GO_CC":
        msigdb_folder
        / "c5.go.cc.v2026.1.Hs.symbols.xlsx",

    "Hallmark":
        msigdb_folder
        / "h.all.v2026.1.Hs.symbols.xlsx"
}


# ============================================================
# 2. CHECK INPUT FILES
# ============================================================

print("=" * 80)
print("INPUT FILE CHECK")
print("=" * 80)

print("ANM file exists:", anm_file.exists())

if not anm_file.exists():
    raise FileNotFoundError(anm_file)

for collection, file_path in collection_files.items():

    print(
        f"{collection}:",
        file_path.exists(),
        "->",
        file_path
    )

    if not file_path.exists():
        raise FileNotFoundError(file_path)


# ============================================================
# 3. LOAD ANM GENES
# ============================================================

anm_df = pd.read_excel(
    anm_file,
    sheet_name="ANM_Unique_Genes"
)

if "Gene" not in anm_df.columns:
    raise ValueError(
        "Could not find 'Gene' column in ANM input."
    )

anm_genes = set(
    anm_df["Gene"]
    .dropna()
    .astype(str)
    .str.strip()
)

anm_genes = {
    gene
    for gene in anm_genes
    if gene not in {
        "",
        "nan",
        "None",
        "NA"
    }
}

print("\nUnique ANM genes:", len(anm_genes))


# ============================================================
# 4. GENE-STRING PARSER
#
# Handles comma, semicolon, or tab-separated gene lists.
# ============================================================

def parse_gene_string(value):

    if pd.isna(value):
        return set()

    values = re.split(
        r"[,;\t]+",
        str(value)
    )

    return {
        gene.strip()
        for gene in values
        if gene.strip()
    }


# ============================================================
# 5. AUTOMATICALLY FIND IMPORTANT COLUMNS
# ============================================================

def detect_columns(df, collection_name):

    # Clean headers
    df.columns = [
        str(c).strip()
        for c in df.columns
    ]

    print(
        f"\n{collection_name} columns:"
    )

    print(df.columns.tolist())


    # --------------------------------------------------------
    # Find Genes column
    # --------------------------------------------------------

    genes_candidates = [
        c for c in df.columns
        if c.lower() in {
            "genes",
            "gene",
            "gene_symbols",
            "gene symbols",
            "members"
        }
    ]

    if not genes_candidates:

        genes_candidates = [
            c for c in df.columns
            if "gene" in c.lower()
            and "count" not in c.lower()
        ]

    if not genes_candidates:

        raise ValueError(
            f"{collection_name}: "
            "Could not identify gene-list column."
        )

    genes_col = genes_candidates[0]


    # --------------------------------------------------------
    # Find pathway / ontology-term column
    # --------------------------------------------------------

    preferred_names = [
        "Pathway",
        "pathway",
        "Gene_Set",
        "Gene Set",
        "Geneset",
        "Name",
        "Term",
        "GO_term",
        "GO Term",
        "GO_term_name",
        "GO Term Name",
        "Ontology_term",
        "Ontology Term"
    ]

    pathway_col = None

    for name in preferred_names:

        if name in df.columns:
            pathway_col = name
            break


    # Search based on text in column name
    if pathway_col is None:

        for c in df.columns:

            low = c.lower()

            if (
                "pathway" in low
                or "term" in low
                or "geneset" in low
                or "gene_set" in low
            ):

                pathway_col = c
                break


    # Final fallback:
    # first column that is not Genes/Description/Gene_Count
    if pathway_col is None:

        excluded = {
            genes_col.lower(),
            "description",
            "gene_count",
            "gene count"
        }

        candidates = [
            c for c in df.columns
            if c.lower() not in excluded
        ]

        if not candidates:

            raise ValueError(
                f"{collection_name}: "
                "Could not identify pathway/term column."
            )

        pathway_col = candidates[0]


    # --------------------------------------------------------
    # Description column (optional)
    # --------------------------------------------------------

    description_col = None

    for c in df.columns:

        if c.lower() == "description":
            description_col = c
            break


    print(
        f"Detected pathway/term column: '{pathway_col}'"
    )

    print(
        f"Detected genes column: '{genes_col}'"
    )

    if description_col:
        print(
            f"Detected description column: '{description_col}'"
        )


    return (
        pathway_col,
        genes_col,
        description_col
    )


# ============================================================
# 6. RUN ORA FOR ONE COLLECTION
# ============================================================

def run_enrichment(
    collection_name,
    file_path,
    anm_genes
):

    print("\n")
    print("=" * 80)
    print("ANALYZING:", collection_name)
    print("=" * 80)


    # --------------------------------------------------------
    # Read MSigDB collection
    # --------------------------------------------------------

    pathways = pd.read_excel(
        file_path
    )


    (
        pathway_col,
        genes_col,
        description_col
    ) = detect_columns(
        pathways,
        collection_name
    )


    print(
        "Source pathways/terms:",
        len(pathways)
    )


    # --------------------------------------------------------
    # Parse gene sets
    # --------------------------------------------------------

    pathways["_Gene_Set"] = (
        pathways[genes_col]
        .apply(parse_gene_string)
    )


    # --------------------------------------------------------
    # Collection-specific background universe
    # --------------------------------------------------------

    background = set()

    for gene_set in pathways["_Gene_Set"]:

        background.update(
            gene_set
        )


    foreground = (
        anm_genes
        & background
    )

    missing_anm = (
        anm_genes
        - background
    )


    N = len(background)
    K = len(foreground)


    print(
        "Background genes:",
        N
    )

    print(
        "ANM genes represented:",
        K
    )

    print(
        "ANM genes not represented:",
        len(missing_anm)
    )


    # --------------------------------------------------------
    # Fisher ORA
    # --------------------------------------------------------

    results = []


    for _, row in pathways.iterrows():


        pathway_name = str(
            row[pathway_col]
        ).strip()


        description = ""

        if description_col is not None:

            value = row[
                description_col
            ]

            if not pd.isna(value):
                description = str(value)


        pathway_genes = (
            row["_Gene_Set"]
            & background
        )


        M = len(pathway_genes)

        if M == 0:
            continue


        overlap = (
            foreground
            & pathway_genes
        )


        a = len(overlap)

        b = K - a

        c = M - a

        d = N - a - b - c


        if min(
            a,
            b,
            c,
            d
        ) < 0:

            raise ValueError(
                f"Invalid contingency table: "
                f"{collection_name} | {pathway_name}"
            )


        # ----------------------------------------------------
        # One-sided Fisher test:
        # Is pathway overrepresented among ANM genes?
        # ----------------------------------------------------

        odds_ratio, p_value = fisher_exact(

            [
                [a, b],
                [c, d]
            ],

            alternative="greater"
        )


        # Expected overlap under null
        expected_overlap = (
            K * M / N
        )


        # Observed / expected
        fold_enrichment = (

            a / expected_overlap

            if expected_overlap > 0

            else np.nan
        )


        results.append({

            "Collection":
                collection_name,

            "Pathway_or_Term":
                pathway_name,

            "Description":
                description,

            "Pathway_size":
                M,

            "ANM_overlap_count":
                a,

            "Expected_overlap":
                expected_overlap,

            "Fold_enrichment":
                fold_enrichment,

            "Odds_ratio":
                odds_ratio,

            "P_value":
                p_value,

            "Overlap_genes":
                "; ".join(
                    sorted(overlap)
                )
        })


    results_df = pd.DataFrame(
        results
    )


    if results_df.empty:

        raise ValueError(
            f"{collection_name}: "
            "No pathways could be tested."
        )


    # ========================================================
    # 7. BH-FDR ACROSS ALL ELIGIBLE PATHWAYS
    #
    # Pathways with 0 overlap remain included.
    # ========================================================

    _, fdr, _, _ = multipletests(

        results_df["P_value"],

        alpha=0.05,

        method="fdr_bh"
    )


    results_df[
        "FDR"
    ] = fdr


    results_df[
        "FDR_LT_0.05"
    ] = (
        results_df["FDR"]
        < 0.05
    )


    results_df[
        "FDR_LT_0.01"
    ] = (
        results_df["FDR"]
        < 0.01
    )


    # --------------------------------------------------------
    # Sort
    # --------------------------------------------------------

    results_df = (

        results_df

        .sort_values(

            by=[
                "FDR",
                "P_value",
                "Fold_enrichment"
            ],

            ascending=[
                True,
                True,
                False
            ]
        )

        .reset_index(
            drop=True
        )
    )


    # --------------------------------------------------------
    # Significant subsets
    # --------------------------------------------------------

    sig_005 = results_df[
        results_df["FDR"] < 0.05
    ]

    sig_001 = results_df[
        results_df["FDR"] < 0.01
    ]


    best = results_df.iloc[0]


    # --------------------------------------------------------
    # Summary row
    # --------------------------------------------------------

    summary = {

        "Collection":
            collection_name,

        "Pathways_tested":
            len(results_df),

        "Background_genes":
            N,

        "Total_ANM_genes":
            len(anm_genes),

        "ANM_genes_in_background":
            K,

        "ANM_genes_not_in_background":
            len(missing_anm),

        "Significant_FDR_LT_0.05":
            len(sig_005),

        "Significant_FDR_LT_0.01":
            len(sig_001),

        "Top_pathway":
            best["Pathway_or_Term"],

        "Top_overlap":
            best["ANM_overlap_count"],

        "Top_fold_enrichment":
            best["Fold_enrichment"],

        "Top_P_value":
            best["P_value"],

        "Top_FDR":
            best["FDR"]
    }


    print("\n----------------------------------------")
    print(collection_name, "SUMMARY")
    print("----------------------------------------")

    print(
        "Pathways tested:",
        len(results_df)
    )

    print(
        "Significant FDR < 0.05:",
        len(sig_005)
    )

    print(
        "Significant FDR < 0.01:",
        len(sig_001)
    )

    print(
        "Top pathway:",
        best["Pathway_or_Term"]
    )

    print(
        "Top P:",
        best["P_value"]
    )

    print(
        "Top FDR:",
        best["FDR"]
    )


    return (
        results_df,
        summary
    )


# ============================================================
# 8. RUN ALL COLLECTIONS
# ============================================================

all_results = {}

summary_rows = []


for (
    collection_name,
    file_path
) in collection_files.items():


    results_df, summary = (
        run_enrichment(

            collection_name,
            file_path,
            anm_genes
        )
    )


    all_results[
        collection_name
    ] = results_df


    summary_rows.append(
        summary
    )


# ============================================================
# 9. CREATE SUMMARY TABLE
# ============================================================

summary_df = pd.DataFrame(
    summary_rows
)


summary_df = (

    summary_df

    .sort_values(

        by=[
            "Significant_FDR_LT_0.05",
            "Top_FDR"
        ],

        ascending=[
            False,
            True
        ]
    )

    .reset_index(
        drop=True
    )
)


# ============================================================
# 10. CREATE COMBINED SIGNIFICANT RESULTS
# ============================================================

significant_tables = []


for (
    collection,
    result_df
) in all_results.items():


    temp = result_df[
        result_df["FDR"] < 0.05
    ].copy()


    if not temp.empty:

        significant_tables.append(
            temp
        )


if significant_tables:

    all_significant = pd.concat(
        significant_tables,
        ignore_index=True
    )


    all_significant = (
        all_significant
        .sort_values(
            [
                "FDR",
                "P_value"
            ]
        )
        .reset_index(
            drop=True
        )
    )

else:

    all_significant = pd.DataFrame(
        columns=[
            "Collection",
            "Pathway_or_Term",
            "Description",
            "Pathway_size",
            "ANM_overlap_count",
            "Expected_overlap",
            "Fold_enrichment",
            "Odds_ratio",
            "P_value",
            "FDR",
            "Overlap_genes"
        ]
    )


# ============================================================
# 11. SAVE ONE EXCEL WORKBOOK
# ============================================================

with pd.ExcelWriter(
    output_file,
    engine="openpyxl"
) as writer:


    summary_df.to_excel(

        writer,

        sheet_name="Summary",

        index=False
    )


    all_significant.to_excel(

        writer,

        sheet_name="All_FDR_Significant",

        index=False
    )


    all_results[
        "Reactome"
    ].to_excel(

        writer,

        sheet_name="Reactome_All",

        index=False
    )


    all_results[
        "GO_BP"
    ].to_excel(

        writer,

        sheet_name="GO_BP_All",

        index=False
    )


    all_results[
        "GO_CC"
    ].to_excel(

        writer,

        sheet_name="GO_CC_All",

        index=False
    )


    all_results[
        "Hallmark"
    ].to_excel(

        writer,

        sheet_name="Hallmark_All",

        index=False
    )


# ============================================================
# 12. FINAL SUMMARY ON SCREEN
# ============================================================

print("\n")
print("=" * 100)
print("FINAL MSigDB ENRICHMENT SUMMARY")
print("=" * 100)


print(

    summary_df[
        [
            "Collection",
            "Pathways_tested",
            "Background_genes",
            "ANM_genes_in_background",
            "Significant_FDR_LT_0.05",
            "Significant_FDR_LT_0.01",
            "Top_pathway",
            "Top_fold_enrichment",
            "Top_P_value",
            "Top_FDR"
        ]
    ]

    .to_string(
        index=False
    )
)


# ============================================================
# 13. PRINT SIGNIFICANT RESULTS
# ============================================================

print("\n")
print("=" * 100)
print("ALL FDR < 0.05 RESULTS")
print("=" * 100)


if all_significant.empty:

    print(
        "No MSigDB pathways reached FDR < 0.05."
    )

else:

    print(

        all_significant[
            [
                "Collection",
                "Pathway_or_Term",
                "ANM_overlap_count",
                "Pathway_size",
                "Expected_overlap",
                "Fold_enrichment",
                "Odds_ratio",
                "P_value",
                "FDR"
            ]
        ]

        .to_string(
            index=False
        )
    )


# ============================================================
# 14. OUTPUT
# ============================================================

print("\n")
print("=" * 100)
print("OUTPUT")
print("=" * 100)

print(
    output_file
)

print(
    "File exists:",
    output_file.exists()
)

if output_file.exists():

    print(
        "File size:",
        round(
            output_file.stat().st_size
            / (1024 ** 2),
            2
        ),
        "MB"
    )


print("\nDONE.")
