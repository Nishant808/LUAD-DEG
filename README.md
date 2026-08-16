# TCGA LUAD — Differential Expression & ML Biomarker Discovery

Pipeline for identifying LUAD (lung adenocarcinoma) biomarkers from TCGA RNA-seq
counts: DESeq2 differential expression → feature matrix construction →
RandomForest/deep-learning classification with SHAP-based biomarker selection.

## Repository layout

```
scripts/     Pipeline stages, run in numeric order (see below)
data/        Input counts and intermediate DEG/feature-matrix outputs
results/     Final biomarker plots and tables
legacy/      Archived earlier iteration of the pipeline (Feb 2025) — not part
             of the current pipeline; kept for reference only
requirements.txt   Python dependencies for the active pipeline
```

## Execution order

1. **`scripts/01_diff_exp_analysis.R`** — DESeq2 differential expression.
   Reads `data/counts/NORMAL/Normal_LUAD_Counts.tsv` and
   `data/counts/TUMOR/Tumor_LUAD_Counts.tsv`, runs each tumor sample against
   the pooled control (normal) samples, writes per-sample DEG tables to
   `data/DEG_RES_ALL/` and the merged `data/LUAD_MERGED_ALL_DEG_RESULT.csv`.
   Dependencies: R packages `DESeq2`, `tibble`, `dplyr`, `purrr`, `writexl`.

2. **`scripts/02_data_process.ipynb`** — Feature matrix construction.
   Filters `data/LUAD_MERGED_ALL_DEG_RESULT.csv` for genes significant
   (p ≤ 0.05) in ≥ 90% of sample comparisons, merges with the raw count
   matrices, and writes `data/FEATURE_DATA_MATRIX.csv` (samples × significant
   genes, labeled tumor/normal).
   Dependencies: `pandas`, `numpy`.

3. **`scripts/03_biomarker_shap.py`** — ML biomarker discovery.
   Trains a RandomForest classifier and a Keras deep neural network on
   `data/FEATURE_DATA_MATRIX.csv`, computes SHAP values (`shap.Explainer`),
   and selects top tumor/normal biomarkers. Writes plots and
   `ML_BIOMARKERS_RES.xlsx` to `results/biomarkers/`.
   Dependencies: see `requirements.txt` (`pandas`, `scikit-learn`, `shap`,
   `numpy`, `tensorflow`, `matplotlib`, `openpyxl`).

## Environment setup

```bash
python3 -m venv BIOENV
source BIOENV/bin/activate   # Windows: BIOENV\Scripts\activate
pip install -r requirements.txt
```

R dependencies (for stage 1) must be installed separately via
`install.packages()` / Bioconductor (`DESeq2`).

## `legacy/step2_2025-02/`

A complete earlier snapshot of this project (Feb 2025), archived as-is.
Its `Data_Process.ipynb` (0.7 significance threshold, vs. 0.9 in the current
pipeline) and `BIOMARKER_SCRIPT_UPDATED.py` (legacy `shap.DeepExplainer`, no
RandomForest stage) are superseded by the scripts in `scripts/`.

It also contains `SCRIPTS/BIOMARKER_PATH_GO_SCRIPT.R`, a GO/KEGG pathway
enrichment stage (biomaRt + enrichR) that has **no counterpart in the current
pipeline** — status undecided. If this analysis is needed for the preprint,
it should be reviewed, have its paths updated to point at
`results/biomarkers/ML_BIOMARKERS_RES.xlsx`, and be promoted to
`scripts/04_pathway_enrichment.R`.
