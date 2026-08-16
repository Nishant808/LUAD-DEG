# Corrected re-run of scripts/02_data_process.ipynb
# CHANGE vs original: the 90%-of-comparisons significance filter is applied
# to DESeq2's BH-adjusted p-value (padj), not the raw p-value. Everything
# else (90% threshold logic, feature matrix construction) is unchanged.
#
# Re-implemented here as a plain script (rather than a notebook) so the
# corrected run can be executed non-interactively and logged.

import pandas as pd
import numpy as np
import os

PROJECT_DIR = "/Users/nishantthalwal/Documents/Data ALL/TCGA_LUAD"
DGE_PATH = os.path.join(PROJECT_DIR, "publication_v1", "data", "LUAD_MERGED_ALL_DEG_RESULT_padj.csv")
NORMAL_DIR = os.path.join(PROJECT_DIR, "data", "counts", "NORMAL")
TUMOR_DIR = os.path.join(PROJECT_DIR, "data", "counts", "TUMOR")
OUT_PATH = os.path.join(PROJECT_DIR, "publication_v1", "data", "FEATURE_DATA_MATRIX.csv")

df = pd.read_csv(DGE_PATH)

# CHANGE: PADJ instead of P_VALUE
PADJ_COL = [col for col in df.columns if col.endswith("_PADJ")]
select_col = ['ID'] + PADJ_COL
df = df[select_col]

filtered_df = df[df[PADJ_COL].le(0.05).sum(axis=1) >= len(PADJ_COL) * 0.9]
print("Total comparisons (tumor samples):", len(PADJ_COL))
print("Total Genes Selected (padj <= 0.05 in >= 90% of comparisons):", filtered_df.shape[0])


def MERGE_DFS(DIR_PATH):
    df_list = []
    for file in os.listdir(DIR_PATH):
        if file.endswith('.tsv'):
            df = pd.read_csv(os.path.join(DIR_PATH, file), sep="\t")
            df_list.append(df)
    merged_df = pd.concat(df_list, axis=1)
    return merged_df


NORMAL_COUNT = MERGE_DFS(NORMAL_DIR)
TUMOR_COUNT = pd.read_csv(os.path.join(TUMOR_DIR, "Tumor_LUAD_Counts.tsv"), sep="\t")

NORMAL_FEA = NORMAL_COUNT.T
TUMOR_FEA = TUMOR_COUNT.T
NORMAL_FEA["Label"] = 0
TUMOR_FEA["Label"] = 1

DATA_MAT = pd.concat([TUMOR_FEA, NORMAL_FEA])

SIG_GENES = list(filtered_df['ID'].values)
DATA_MAT_1 = DATA_MAT[SIG_GENES]
DATA_MAT_1["Label"] = DATA_MAT["Label"]

DATA_MAT_1.to_csv(OUT_PATH)
print("Feature matrix shape:", DATA_MAT_1.shape)
print("Saved to:", OUT_PATH)
