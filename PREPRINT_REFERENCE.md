# TCGA-LUAD DEG/ML Project — Full Reference for Preprint Drafting

**Purpose of this document:** this is a complete, fact-checked reference to the
entire project — every file, every parameter, every number, every known issue —
written so it can be handed to another LLM (alongside the actual result/data
files, except the raw per-sample count matrices) to draft a bioRxiv preprint.
Every number in this document was directly read from a script, a log, or a
result file during the session that produced it; nothing here is inferred,
estimated, or reconstructed from memory. Where a value could not be directly
verified, that is stated explicitly rather than guessed.

**If you are the LLM drafting the preprint from this document: read the
"Which pipeline run is authoritative" section first.** This project contains
three separate pipeline executions (one current/correct, one deprecated, one
archived/legacy) with different gene counts and different files. Using the
wrong one will produce a preprint with incorrect statistics.

---

## 1. Which pipeline run is authoritative

| Run | Location | Status | Gene-selection basis |
|---|---|---|---|
| **`publication_v1/`** | project root | **✅ USE THIS ONE** | BH-adjusted p-value (padj) — statistically correct |
| Original run | `scripts/`, `data/`, `results/` (project root) | ❌ Deprecated — do not use for the preprint | Raw p-value — statistically flawed (see §6.1) |
| Legacy run | `legacy/step2_2025-02/` | ❌ Archived, superseded — do not use | Raw p-value, different threshold (0.7, not 0.9), different SHAP method |

**All sample sizes, gene counts, model metrics, marker lists, and pathway
enrichment results that should go into the preprint come from
`publication_v1/`.** The original and legacy runs are documented below only
for provenance/completeness and to explain why they must not be cited.

---

## 2. One-paragraph project summary

This project identifies gene-expression differences between lung
adenocarcinoma (LUAD) tumor tissue and normal lung tissue using TCGA RNA-seq
count data (59 normal samples, 539 tumor samples), and asks whether a small
set of genes can both (a) statistically distinguish the two groups (DESeq2
differential expression) and (b) be used to build machine-learning classifiers
(Random Forest, a small feed-forward neural network) that predict tumor vs.
normal status from expression alone, with SHAP used to identify which genes
drive those classifiers' predictions in each direction. A GO/KEGG pathway
enrichment analysis was then run on the SHAP-identified gene markers.

---

## 3. Complete file inventory

### 3.1 `publication_v1/` — the authoritative run

```
publication_v1/
├── README.md                                    (a project-internal README for this run — this document supersedes it for preprint purposes but they should not contradict each other)
├── scripts/
│   ├── 01_diff_exp_analysis_padj.R               DESeq2, padj-based, parallelized (§5.1)
│   ├── 02_data_process_padj.py                   Feature matrix construction, padj-based (§5.2)
│   ├── 03_biomarker_shap_padj.py                 RF + NN + SHAP, with fixed marker selection (§5.3)
│   ├── 04_pathway_go_enrichment.R                GO/KEGG enrichment on SHAP markers (§5.4)
│   └── 05_qc_diagnostics.R                       PCA/MA/volcano/dispersion plots (§5.5)
├── data/
│   ├── LUAD_MERGED_ALL_DEG_RESULT_padj.csv       812 MB. 60,660 gene rows × 1,079 columns (ID + 539 tumor samples × [LOG2FC, PADJ]). Output of script 01.
│   ├── DEG_RES_ALL/                              539 files, one per tumor sample: `<TCGA-barcode>_vs_Control_DEG_RESULT.xlsx`, each with sheets ALL_DEG / SIGNIFICANT / UPREGULATED / DOWNREGULATED for that one tumor-vs-pooled-normal comparison. Output of script 01.
│   └── FEATURE_DATA_MATRIX.csv                   263 KB. 598 rows (samples) × 112 columns (unnamed index + 110 genes + Label). Output of script 02.
├── results/
│   ├── ML_BIOMARKERS_RES.xlsx                    Sheets: RANDOM_FOREST_BIOMARKERS (110 rows, all genes ranked by RF importance), SHAP_BIOMARKERS (20 rows × 2 cols: TUMOR_MARKER, NORMAL_MARKER). Output of script 03.
│   ├── rf_metrics.txt                            Random Forest test-set accuracy + classification_report. Output of script 03.
│   ├── nn_metrics.txt                             Keras NN test-set loss + accuracy. Output of script 03.
│   ├── nn_training_history.csv                   51 rows (header + 50 epochs) × 5 cols (epoch, loss, accuracy, val_loss, val_accuracy). Output of script 03.
│   └── pathway_go/
│       ├── SHAP_TUMOR_PATHWAY_GO_RES.xlsx        4 sheets: KEGG (32 terms), BIOLOGICAL_PRO (320 terms), MOLECULAR_FUN (33 terms), CELLULAR_COMP (49 terms). Output of script 04.
│       ├── SHAP_NORMAL_PATHWAY_GO_RES.xlsx       4 sheets: KEGG (7 terms), BIOLOGICAL_PRO (144 terms), MOLECULAR_FUN (15 terms), CELLULAR_COMP (43 terms). Output of script 04.
│       └── SHAP_Biomarkers_Gene_Symbols.xlsx     Ensembl ID → HGNC symbol mapping for both marker sets (20/20 mapped, 0 dropped). Output of script 04.
├── figures/
│   ├── Top_MARKERS_RandomForest.png              Bar plot, top 20 of 110 genes by RF importance
│   ├── TUMOR_FEATURE_IMP_NN_SELECTED.png         SHAP beeswarm plot, 20 tumor-push markers
│   ├── NORMAL_FEATURE_IMP_NN_SELECTED.png        SHAP beeswarm plot, 20 normal-push markers
│   ├── PCA_plot.png                              PCA of all 598 samples, colored by Group (Control/Treated)
│   ├── MA_plot.png                               DESeq2 MA plot, combined single model across all 598 samples
│   ├── volcano_plot.png                          Volcano plot (log2FC vs. -log10(padj)), combined single model
│   ├── dispersion_plot.png                       DESeq2 dispersion estimates plot, combined single model
│   └── pathway_go/                               16 PNGs: {SHAP_TUMOR, SHAP_NORMAL} × {KEGG, BP, MF, CC} × {DOT, BAR}
└── environment/
    ├── r_sessionInfo.txt                         R sessionInfo() at time of this run (§7)
    ├── python_pip_freeze.txt                     pip freeze from BIOENV at time of this run (§7)
    └── stage1_run.log … stage5_qc_run.log         Raw console output of each script, for provenance
```

### 3.2 Original run (deprecated — project root `scripts/`, `data/`, `results/`)

```
scripts/01_diff_exp_analysis.R          DESeq2, raw-pvalue-based (statistically flawed, see §6.1)
scripts/02_data_process.ipynb           Feature matrix construction, raw-pvalue-based, 90% threshold
scripts/03_biomarker_shap.py            RF + NN + SHAP, original (unfixed) marker-selection logic
data/LUAD_MERGED_ALL_DEG_RESULT.csv     987 MB. Same shape/structure as the padj version, but PVALUE not PADJ columns. Dated 2025-01-23 (pre-dates this project's active session by ~1 year; not regenerated during this session).
data/FEATURE_DATA_MATRIX.csv            700 KB. 598 rows × 311 columns (index + 309 genes + Label). Dated 2025-02-04.
data/counts/NORMAL/Normal_LUAD_Counts.tsv   10 MB, 59 samples × 60,660 genes (raw counts, shared source for every run in this project)
data/counts/TUMOR/Tumor_LUAD_Counts.tsv     86 MB, 539 samples × 60,660 genes (raw counts, shared source for every run in this project)
results/biomarkers/ML_BIOMARKERS_RES.xlsx   Sheets: RANDOM_FOREST_BIOMARKERS (309 rows), SHAP_BIOMARKERS (20 rows — see note below). Dated 2026-06-24.
results/biomarkers/Top_MARKERS_RandomForest.png, TUMOR_FEATURE_IMP_NN_SELECTED.png, NORMAL_FEATURE_IMP_NN_SELECTED.png
```

**Note on the original run's SHAP list:** unlike the first attempt at the
padj-based re-run (which produced 0 markers, see §6.3), the original
309-gene, raw-p-value run's `SHAP_BIOMARKERS` sheet has 20 data rows (i.e.
it is populated) — because the original abs-magnitude + cross-class-eviction
algorithm (see §6.3) happens to work when there are enough candidate genes.
This is *not* a reason to use the original run: its underlying 309-gene list
is itself downstream of the flawed raw-p-value filtering (§6.1), so its SHAP
markers reflect that same statistical error and must not be cited either.
The **only** valid SHAP marker list in this project is in `publication_v1/`.

No per-sample `DEG_RES_ALL`-equivalent directory exists for the original run
— only the final merged CSV and feature matrix survive; no per-tumor-sample
xlsx breakdown was found anywhere under `data/` or `results/` at the project
root.

### 3.3 Legacy run (archived, superseded — `legacy/step2_2025-02/`)

An entire earlier, self-contained iteration of this project, dated
2025-02-11 to 2025-02-15 (~16 months before `publication_v1/`). Uses:
- `Data_Process.ipynb`: **0.7** comparison threshold (not 0.9), raw p-value — produced 2,046 genes (this exact count is printed in that notebook's own saved cell output).
- `BIOMARKER_SCRIPT_UPDATED.py`: `shap.DeepExplainer` (not `shap.Explainer`), **no Random Forest stage**, `epochs=100` (not 50).
- `BIOMARKER_PATH_GO_SCRIPT.R`: the direct ancestor of `publication_v1/scripts/04_pathway_go_enrichment.R` — same analysis logic, adapted in this session (§5.4) to point at `publication_v1`'s marker list instead.
- `DATA/COUNTS/NORMAL/` contains **an extra file, `Normal_LUSC_Counts.txt`** — LUSC (lung *squamous cell* carcinoma) counts, a different cancer type, mixed into this legacy folder's normal-counts directory. This file is not used by any script (only `Normal_LUAD_Counts.txt`/`.tsv` is read), but its presence is worth knowing about if anyone browses this folder — it is not part of any actual analysis in this project.
- GO enrichment output row counts (legacy, for reference only, **not from `publication_v1`**): `GO_BP_Normal.csv` 14 terms, `GO_Tumor_BP.csv` 494 terms, `Normal_KEGG.csv` 14 terms, `Tumor_KEGG.csv` 14 terms. These were computed from the legacy `shap.DeepExplainer`-derived marker list, not from `publication_v1`'s markers.
- `Neural Network.drawio` — a manually-authored architecture diagram source file (not a script output); could be exported as a supplementary figure if wanted, but was not used in this project.
- `README.txt`, `requirements.txt` — original venv-setup notes for this legacy iteration.

### 3.4 Other root-level files

- `README.md` — human-facing project README (execution order, folder layout); written before the padj correction, so does not describe `publication_v1/`.
- `requirements.txt` — Python dependency pins, identical across the whole project (pandas 2.0.3, scikit-learn 1.3.2, shap 0.44.1, numpy 1.24.3, tensorflow 2.13.0, matplotlib 3.7.5, openpyxl 3.1.5).
- `environment/` (root level) — `python_freeze.txt`, `python_interpreter.txt`, `r_packages.csv`: environment snapshot taken when this project's Python/R environment was first set up in this session, before `publication_v1` existed. `publication_v1/environment/` (§3.1) is the environment snapshot for the actual authoritative run and should be preferred if the two ever disagree.
- `BIOENV/` — the Python 3.10.15 virtual environment used to run every Python script in this project (details in §7). Not analysis content.
- `.Rhistory` (root and `publication_v1/`) — R console session history files, not analysis content.
- `.idea/` — JetBrains IDE project config, not analysis content.

---

## 4. Sample sizes and data provenance (verified, shared across every run)

- **N normal (Control) = 59** — `data/counts/NORMAL/Normal_LUAD_Counts.tsv`, 59 sample columns, each a TCGA barcode (e.g. `TCGA-44-2661-11A-01R-1758-07`).
- **N tumor (Treated) = 539** — `data/counts/TUMOR/Tumor_LUAD_Counts.tsv`, 539 sample columns.
- **N total = 598.**
- **60,660 genes** measured per sample (Ensembl gene IDs, versioned, e.g. `ENSG00000000460.17`) — row count of both raw count files, and the gene count DESeq2 itself reports (`N genes = 60660`, printed directly by the combined-model DESeq2 run, script 05 log).
- Both count files are shared, unmodified, raw input across the original run, the legacy run, and `publication_v1` — only the downstream statistical processing differs between runs. **No sample was added, removed, or resampled between runs.**
- **No clinical or sample metadata (age, sex, tumor stage, batch, etc.) exists anywhere in this project as a saved table.** The only sample-level annotation used anywhere is a binary `Group` label (`Control` vs. `Treated`), constructed purely from which count file (`NORMAL` vs. `TUMOR`) a sample came from. If the preprint needs cohort characteristics beyond N and tissue type, that data does not currently exist in this project and would need to be sourced separately from TCGA clinical data.
- Sample identifiers are standard TCGA barcodes (format `TCGA-XX-XXXX-XXX-XXR-XXXX-XX`), present in the raw count file headers and preserved through to the per-sample DEG result filenames in `publication_v1/data/DEG_RES_ALL/`.

---

## 5. Methods — exact parameters, `publication_v1/` (authoritative run)

### 5.1 Stage 1 — Differential expression (`scripts/01_diff_exp_analysis_padj.R`)

- Tool: **DESeq2** (Bioconductor, v1.52.0 in this run's environment).
- Design: `DESeqDataSetFromMatrix(design = ~Group)`, with `dds$Group <- relevel(dds$Group, ref = "Control")`.
- **Structure: one DESeq2 run per tumor sample**, each comparing that single tumor sample against the pooled 59 normal samples (i.e. **539 independent DESeq2 fits**, not one combined model). This is a deliberate design choice inherited unchanged from the original script.
- Significance filter: **`padj <= 0.05`** (DESeq2's Benjamini-Hochberg FDR-adjusted p-value) — this is the corrected threshold; see §6.1 for why the original script's `pvalue <= 0.05` was wrong.
- Direction thresholds: up-regulated = `log2FoldChange >= 1`; down-regulated = `log2FoldChange <= -1`. Unchanged from the original.
- Per-sample output: an xlsx with sheets `ALL_DEG` / `SIGNIFICANT` / `UPREGULATED` / `DOWNREGULATED`, written to `publication_v1/data/DEG_RES_ALL/<sample>_vs_Control_DEG_RESULT.xlsx` (539 files total).
- Merged output: for every sample, `ID`, `log2FoldChange`, and `padj` are kept and column-renamed to `Control_vs_<sample>_LOG2FC` / `Control_vs_<sample>_PADJ`, then inner-joined across all 539 samples on gene `ID` → `publication_v1/data/LUAD_MERGED_ALL_DEG_RESULT_padj.csv` (60,660 rows × 1,079 columns).
- **Execution note (not an analysis change):** this stage was parallelized across the 539 independent per-sample DESeq2 calls using R's `parallel::mclapply(mc.cores = 4)`. The core count was chosen empirically on this machine (Apple M4: 4 Performance cores + 6 slower Efficiency cores) — timed tests on a 20-sample batch gave 154.6s at `mc.cores=4`, 202.6s at `mc.cores=6`, and 174.7s at `mc.cores=10` (P/E-core scheduling contention made the naive "use all cores" setting slower, not faster). This is purely a wall-clock optimization; it does not change DESeq2's inputs, per-sample fitting, or outputs — results are identical to what a serial loop would produce. Total stage 1 runtime for all 539 samples: **102.96 minutes**.
- A technical detail worth knowing for the Methods section: the R script reads the count files with `read.table(file, sep="\t")` without an explicit `header` argument. R's default behavior in this situation (not simply `header=FALSE`) is to auto-detect a header row when the first line has exactly one fewer field than subsequent lines — which is the case here (TCGA count files have no explicit gene-ID column header, so the header row is one field shorter than the data rows). This auto-detection correctly recovers the sample barcode column names; it is standard R behavior, not a script bug, and was left unchanged.

### 5.2 Stage 2 — Feature matrix construction (`scripts/02_data_process_padj.py`)

- Loads the stage-1 merged output, extracts the 539 `*_PADJ` columns.
- **Filter: a gene is kept if `padj <= 0.05` in at least 90% of the 539 tumor-vs-control comparisons** (`.le(0.05).sum(axis=1) >= len(PADJ_COL) * 0.9`). Unchanged threshold logic from the original script; only the column being thresholded changed (padj vs. raw pvalue).
- Result: **110 genes** pass this filter (printed directly by the script: `Total Genes Selected (padj <= 0.05 in >= 90% of comparisons): 110`), down from 309 in the original raw-p-value run (a ~64% reduction — see §6.1 for why).
- These 110 genes' raw counts (from the same two count files as every other run) are transposed, labeled (`Label=0` for normal, `Label=1` for tumor), concatenated, and written to `publication_v1/data/FEATURE_DATA_MATRIX.csv`: **598 rows × 110 gene columns + Label** (112 raw CSV columns including the unnamed sample-ID index column).
- Re-implemented as a plain `.py` script (the original was a Jupyter notebook, `scripts/02_data_process.ipynb`) so this corrected run could be executed non-interactively and logged; the filtering logic itself is unchanged apart from padj vs. raw p-value.

### 5.3 Stage 3 — ML classification & SHAP (`scripts/03_biomarker_shap_padj.py`)

- Input: the 110-gene feature matrix above. Features standardized with `sklearn.preprocessing.StandardScaler` (zero mean, unit variance).
- Train/test split: `train_test_split(test_size=0.2, random_state=42)` → **120 test samples (16 normal, 104 tumor)**, 478 training samples.

**Random Forest:**
- `RandomForestClassifier(n_estimators=100, random_state=52)`.
- Test-set performance (`publication_v1/results/rf_metrics.txt`):
  - **Accuracy: 0.975**
  - Class 0 (normal, n=16): precision 0.88, recall 0.94, F1 0.91
  - Class 1 (tumor, n=104): precision 0.99, recall 0.98, F1 0.99
  - Macro avg F1 0.95, weighted avg F1 0.98
  - `RandomForestClassifier` is deterministic on CPU with a fixed seed — this result is exactly reproducible across runs.
- Full ranked importance list for all 110 genes: `ML_BIOMARKERS_RES.xlsx`, sheet `RANDOM_FOREST_BIOMARKERS`. Top 5 (for illustration; full list is in the file): `ENSG00000150625.16` (0.0948), `ENSG00000182010.11` (0.0547), `ENSG00000224215.1` (0.0503), `ENSG00000169252.6` (0.0433), `ENSG00000101280.8` (0.0426).
- Figure: `figures/Top_MARKERS_RandomForest.png` (top 20 of 110 by importance, horizontal bar plot).

**Neural network:**
- Architecture: `Sequential([Dense(128, relu), Dropout(0.6), Dense(64, relu), Dropout(0.6), Dense(32, relu), Dropout(0.6), Dense(16, relu), Dropout(0.6), Dense(1, sigmoid)])`.
- Compiled with `optimizer='adam'`, `loss='binary_crossentropy'`, `metrics=['accuracy']`.
- Trained with `epochs=50, batch_size=32, validation_split=0.2` (i.e. 20% of the 478 training samples held out for validation during training, separate from the 120-sample test set).
- `tf.random.set_seed(42)` was added to this run (the original script had no TensorFlow seed at all), in addition to the sklearn seeds above.
- Full 50-epoch training/validation loss and accuracy history: `publication_v1/results/nn_training_history.csv`.
- Test-set performance (`publication_v1/results/nn_metrics.txt`, via `model.evaluate()`, never called in the original script): **test loss 0.0957, test accuracy 0.975**.
- **Reproducibility caveat, directly observed in this project:** an earlier run of this script produced **test loss 0.0743, test accuracy 0.9667** — different from the numbers above. The model-definition, `tf.random.set_seed(42)`, training (`model.fit`), and `model.evaluate()` code was byte-identical between that run and the final one (both execute before the code that was later edited); only the downstream SHAP marker-selection logic, which runs *after* `model.evaluate()`, was changed in between (§ SHAP fix below) — so the code that actually produced these two different NN numbers did not change. Same data both times. `RandomForestClassifier` reproduced identically in both runs; the Keras model did not. `tf.random.set_seed()` seeds TensorFlow's global RNG but does not by itself guarantee bit-exact determinism for every op on CPU (thread-level floating-point reduction order can still vary between runs). If exact reproducibility is required for the preprint, further determinism settings (e.g. `tf.config.experimental.enable_op_determinism()`, single-threaded execution) would be needed — these were not added in this project. **The values reported in this document and in the current `publication_v1/results/` files (0.0957 / 0.975) are the final, current ones; treat 96.7–97.5% accuracy as the reproducible range, not a single bit-exact number.**

**SHAP explanation and marker selection — includes a bug found and fixed during this project:**
- SHAP values computed via `shap.Explainer(model, X_train_np)`, applied to the test set via `explainer(X_test_np, max_evals=1000)` (Permutation explainer, generic API, chosen in an earlier revision of this script to work with modern Keras/TF execution styles).
- **Original selection algorithm (used in the very first version of this padj-based run, and unchanged from the original raw-p-value script):** builds two marker lists by walking two SHAP-score-sorted gene rankings — one ranked by mean **absolute** SHAP value on tumor-labeled test samples, one on normal-labeled test samples — in lockstep. Whenever a gene appears near the top of both rankings, it is added to one list then evicted if it shows up in the other (mutual exclusivity). **On this run's 110-gene pool, this produced zero markers in both lists** (`SHAP_BIOMARKERS` sheet had 0 data rows; both SHAP figures were blank).
- **Root cause, established by direct investigation** (re-running the model deterministically and tracing the selection loop iteration-by-iteration): the tumor-subset and normal-subset SHAP magnitude rankings are highly correlated — **Pearson r = 0.87** across all 110 genes, with 70–88% top-K overlap at every K from 10 to 110 (both rankings come from explaining the same single binary classifier, just evaluated on differently-labeled subsets of test points). The traced loop shows both marker lists actually **peaked at 12/12 genes around iteration 63–90** of 110, before the tail-end add/evict churn on lower-ranked, noisier genes erased both lists back down to 0 by the final iteration (110 gene-pairs exhausted without ever reaching the 20-per-class break condition).
- **Fix applied (with explicit sign-off):** rank by **signed** mean SHAP contribution instead of absolute magnitude, and drop the forced cross-class exclusivity. Specifically: `TUMOR_MARKER` = genes with **positive** mean SHAP value among tumor-labeled test samples (i.e. the gene pushes the model's prediction toward "tumor"), ranked by strongest positive push, top 20. `NORMAL_MARKER` = genes with **negative** mean SHAP value among normal-labeled test samples (pushes toward "normal"), ranked by strongest negative push, top 20. A gene can legitimately appear in **both** lists if it's a strong discriminator in both directions — no exclusivity is enforced.
- This produced **20/20 markers**, as intended. Final marker lists (`publication_v1/results/ML_BIOMARKERS_RES.xlsx`, sheet `SHAP_BIOMARKERS`):

  **TUMOR_MARKER (20 genes, ranked by strongest positive SHAP push toward tumor):**
  `ENSG00000234281.5, ENSG00000168484.12, ENSG00000183019.7, ENSG00000078081.8, ENSG00000224215.1, ENSG00000108576.10, ENSG00000170323.9, ENSG00000163815.6, ENSG00000167434.10, ENSG00000150625.16, ENSG00000022267.19, ENSG00000157654.19, ENSG00000135604.10, ENSG00000073754.6, ENSG00000182010.11, ENSG00000134917.10, ENSG00000204305.14, ENSG00000169252.6, ENSG00000105974.13, ENSG00000170989.10`

  **NORMAL_MARKER (20 genes, ranked by strongest negative SHAP push toward normal):**
  `ENSG00000234281.5, ENSG00000168484.12, ENSG00000108576.10, ENSG00000150625.16, ENSG00000224215.1, ENSG00000182010.11, ENSG00000073754.6, ENSG00000022267.19, ENSG00000163815.6, ENSG00000078081.8, ENSG00000135063.19, ENSG00000170323.9, ENSG00000183019.7, ENSG00000123119.12, ENSG00000102104.9, ENSG00000157654.19, ENSG00000176771.17, ENSG00000135604.10, ENSG00000161649.13, ENSG00000109205.16`

  (Genes appearing in both lists, e.g. `ENSG00000234281.5`, `ENSG00000224215.1`, `ENSG00000150625.16`, are consistently discriminative in both directions — this is a meaningful result under the fixed method, not an error.)
- Figures: `figures/TUMOR_FEATURE_IMP_NN_SELECTED.png`, `figures/NORMAL_FEATURE_IMP_NN_SELECTED.png` — SHAP beeswarm plots for the 20 markers per class; both confirmed visually to show proper distributions (not blank).
- **This selection-logic change is the one deliberate methods deviation from the original script in this entire project, and was made only after this specific fix was explicitly reviewed and approved.** If the preprint's Methods section describes SHAP marker selection, it should describe the signed-contribution method above — not the original magnitude+exclusivity method, which is not what produced the numbers in this project's final results.

### 5.4 Stage 4 — GO/KEGG pathway enrichment (`scripts/04_pathway_go_enrichment.R`)

- Input: the 20 `TUMOR_MARKER` + 20 `NORMAL_MARKER` genes above (analyzed as two separate 20-gene lists).
- Ensembl gene IDs mapped to HGNC gene symbols via `biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', ...)`, connecting via `useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")`. **20/20 tumor genes and 20/20 normal genes mapped successfully — no genes were dropped.**
- **Infrastructure note:** the direct ancestor script (`legacy/.../BIOMARKER_PATH_GO_SCRIPT.R`) used `useMart("ensembl", ...)`, which returned an HTTP 404 in this environment (Ensembl's BioMart discovery endpoint has moved since that script was written). Replaced with `useEnsembl()`, the current recommended connection function, which also automatically retried via an Asia mirror when the primary Ensembl site was unresponsive during this run. This is a connection-layer fix only — the underlying `getBM()` query (attributes, filters, values) is identical.
- Enrichment performed via `enrichR::enrichr()` against four Enrichr libraries: `KEGG_2021_Human`, `GO_Biological_Process_2023`, `GO_Molecular_Function_2023`, `GO_Cellular_Component_2023`.
- Results (non-empty for every database, both marker sets):

  | Marker set | KEGG terms | GO Biological Process terms | GO Molecular Function terms | GO Cellular Component terms |
  |---|---|---|---|---|
  | SHAP_TUMOR (20 genes) | 32 | 320 | 33 | 49 |
  | SHAP_NORMAL (20 genes) | 7 | 144 | 15 | 43 |

  Full term-level results (term name, P-value, adjusted P-value, overlap, combined score, genes) are in `publication_v1/results/pathway_go/SHAP_TUMOR_PATHWAY_GO_RES.xlsx` and `SHAP_NORMAL_PATHWAY_GO_RES.xlsx`, one sheet per database. These are reported here as generated, without biological interpretation — the term tables and plots are the primary source for deciding what's worth discussing.
- Plots: for each of the 8 (marker set × database) combinations, a dot plot and a bar plot of the top 20 terms by raw P-value — 16 PNGs total in `publication_v1/figures/pathway_go/`.
- **Cosmetic labeling issue, inherited unchanged from the legacy script:** the dot plots' color-legend is labeled "p.adj" but the plotted/colored value is actually Enrichr's raw `P.value` column — the plotting code never references a separately adjusted-p column. This is a plot-label-only issue; the underlying xlsx result tables contain Enrichr's own `P.value` and `Adjusted.P.value` columns correctly, both usable if the preprint needs adjusted values.
- RF-ranked genes were **not** run through this enrichment step (only the SHAP markers were) — this matches the legacy script, which had an RF-enrichment call present but commented out and never executed, in every version of this project.

### 5.5 QC diagnostics (`scripts/05_qc_diagnostics.R`) — supplementary, not part of the DEG-calling methodology

- A **separate, single combined DESeq2 model** across all 598 samples at once (`design = ~Group`, `ref = "Control"`) — distinct from the 539 independent per-sample comparisons in stage 1, which remain the basis for all actual DEG calls and gene selection. This combined model exists only to produce standard whole-cohort QC figures; it is not used anywhere else in the pipeline.
- Figures (`publication_v1/figures/`):
  - `PCA_plot.png` — `vst()`-transformed PCA of all 598 samples, colored by Group. Visually confirmed: tight, well-separated Control cluster vs. a much more dispersed Treated cluster.
  - `MA_plot.png` — DESeq2's `plotMA()`. **Its default blue/"significant" coloring uses DESeq2's own default `alpha = 0.1`, not this project's actual 0.05 padj threshold used everywhere else.** Cosmetic only; does not affect any filtering or the volcano plot.
  - `volcano_plot.png` — custom `ggplot2` volcano plot, explicitly using `padj <= 0.05` and `|log2FoldChange| >= 1` for the significance coloring (consistent with the actual thresholds used in stages 1–2). A warning was logged that 15,179 rows (genes) were dropped from the plot due to `NA` padj values — these are genes DESeq2's independent filtering excluded (typically very low mean count), a standard and expected DESeq2 behavior, not a script error.
  - `dispersion_plot.png` — DESeq2's `plotDispEsts()`, standard shape confirmed visually (gene-wise estimates shrunk to fitted curve, few outliers circled).

---

## 6. Known issues, limitations, and caveats (comprehensive)

### 6.1 Why the original run's raw p-value filtering was wrong (the reason `publication_v1` exists)

The original scripts (`scripts/01_diff_exp_analysis.R`, `scripts/02_data_process.ipynb`) filtered on DESeq2's **raw p-value** (`pvalue <= 0.05`) at two separate steps: per-sample DEG calling, and the 90%-of-comparisons feature-selection threshold. With **60,660 genes** tested in every one of the 539 comparisons, this is a multiple-testing problem: at raw p ≤ 0.05, roughly 3,000 genes are expected to look "significant" by chance alone in any single comparison, before any real biological signal is considered. This inflated the apparent significant-gene count (309 genes survived the original run's filter) and undermines any claim about "significant" biomarkers built on it. `publication_v1` filters on **padj** (BH/FDR-adjusted p-value) instead, at both steps, and 110 genes survive — everything else in the DESeq2 design and downstream pipeline is unchanged (§5.1–5.2).

### 6.2 Pre-existing path/format mismatch in the original scripts (fixed as a path correction, not an analysis change)

Before this session, `scripts/01_diff_exp_analysis.R` and `scripts/02_data_process.ipynb` referenced flat files named `Normal_LUAD_Counts.txt` / `Tumor_LUAD_Counts.txt` directly under a `COUNTS/` directory, while the actual data on disk was one level deeper, in `NORMAL/`/`TUMOR/` subdirectories, with a `.tsv` extension. This mismatch pre-dates this project's active session (the root `data/` output files are dated January–February 2025, roughly a year before this session's work) — it was corrected as part of a folder reorganization earlier in this project (updating file locations/extensions the scripts look for), not as part of any statistical/analysis change. `publication_v1`'s scripts read from the correct, current locations from the start.

### 6.3 SHAP marker selection bug — see §5.3 for the full account

Summarized here for completeness: the original (unmodified) SHAP marker-selection algorithm produced 0 markers on `publication_v1`'s 110-gene pool, root-caused to high correlation (r=0.87) between the two class-conditioned SHAP rankings, and fixed by switching to signed-SHAP-contribution ranking with no forced cross-class exclusivity. The original 309-gene run's SHAP list is *not* empty (it has 20/20 markers, using the same flawed selection algorithm) but must not be cited, because it sits downstream of the flawed raw-p-value gene selection (§6.1) — see §3.2.

### 6.4 No model performance metrics existed until this project's fixes

The original `scripts/03_biomarker_shap.py` imported `classification_report`/`accuracy_score` but never called them, and never called `model.evaluate()` on the Keras model — so the original run has **no accuracy, precision, recall, F1, or loss figures anywhere**, for either model. `publication_v1/scripts/03_biomarker_shap_padj.py` added both (§5.3); those are the only model-performance numbers that exist in this project.

### 6.5 Class imbalance

538/539 tumor samples vs. 59 normal samples at the full-cohort level (~9:1), and 104/16 in the specific test split used for stage 3 (§5.3). The Random Forest's 0.88 precision / 0.94 recall for the normal class (n=16 in test) should be read with this small denominator in mind — headline accuracy (driven mostly by the much larger tumor class) is not, on its own, a reliable summary of performance on the minority (normal) class.

### 6.6 NN reproducibility (see §5.3 for full detail)

`tf.random.set_seed(42)` reduces but does not eliminate run-to-run variance in the Keras model's results; a same-seed, same-code re-run produced test accuracy 0.9667 vs. the final reported 0.975 (test loss 0.0743 vs. 0.0957). Treat NN performance as "approximately 97% on held-out data," not a single bit-exact figure, unless the preprint explicitly wants to describe this as a range.

### 6.7 No clinical/sample metadata (see §4)

Only a binary Group (Control/Treated) label exists for any sample anywhere in this project. Age, sex, tumor stage, smoking status, batch, or any other clinical covariate is not available in any file in this project.

### 6.8 GO enrichment only covers SHAP markers, not the full RF ranking

Pathway enrichment (§5.4) was run only on the 20+20 SHAP marker genes, not on the full 110-gene RF-ranked list (or any subset of it) — this matches every prior version of the GO enrichment script in this project's history, none of which ever executed an RF-based enrichment call (present in the legacy script as commented-out code, never run).

### 6.9 Cosmetic/labeling issues that do not affect any data or statistics

- `MA_plot.png` colors by DESeq2's default `alpha=0.1`, not the project's actual 0.05 threshold (§5.5).
- The pathway_go dot plots' "p.adj" color-legend label actually shows raw `P.value`, not an adjusted value (§5.4). The underlying data tables have both columns correctly.

### 6.10 DESeq2 independent filtering and `NA` padj values

DESeq2's independent filtering assigns `NA` to `padj` (and sometimes `pvalue`) for genes it excludes from multiple-testing correction (typically very low mean normalized count, or extreme Cook's-distance outliers). Every `padj <= 0.05` comparison in this project (R's `subset()`, pandas' `.le()`) correctly treats `NA` as "not significant" (both silently exclude `NA` rows/values rather than erroring) — this is standard, expected DESeq2/analysis behavior, not a bug, and is why, e.g., the volcano plot (§5.5) reports 15,179 genes dropped from that specific figure.

---

## 7. Environment and reproducibility

**R:** version 4.6.1 (aarch64-apple-darwin23), running under macOS Tahoe 26.5.2, on an Apple M4 (4 Performance cores + 6 Efficiency cores, 16 GB RAM). Key package versions: DESeq2 1.52.0, tibble 3.3.1, dplyr 1.2.1, purrr 1.2.2, writexl 2.0.0, ggplot2 4.0.3, stringr 1.6.0, readxl 1.5.0, biomaRt 2.68.0, enrichR 3.4. Full `sessionInfo()` output: `publication_v1/environment/r_sessionInfo.txt`.

**Python:** version 3.10.15, run inside an isolated virtual environment (`BIOENV/`, created via `python3 -m venv` from a pre-existing Python 3.10 interpreter — the machine's default Python was 3.13, incompatible with the pinned `tensorflow==2.13.0`). Key package versions (exactly matching `requirements.txt`): pandas 2.0.3, numpy 1.24.3, scikit-learn 1.3.2, shap 0.44.1, tensorflow 2.13.0 (with `tensorflow-macos` 2.13.0 auto-installed as the arm64 CPU backend — no GPU/Metal acceleration), matplotlib 3.7.5, openpyxl 3.1.5. Full `pip freeze` output: `publication_v1/environment/python_pip_freeze.txt`.

**Random seeds set:** `train_test_split(random_state=42)`, `RandomForestClassifier(random_state=52)`, `tf.random.set_seed(42)` (the last one added during this project; see §5.3/§6.6 for its reproducibility limits).

**Internet-dependent steps:** stage 4 (§5.4) requires live connectivity to Ensembl's BioMart service and the Enrichr web API (`https://maayanlab.cloud/Enrichr/`) — both were reachable and used successfully during this run, but re-running stage 4 later requires the same external services to be available.

---

## 8. Suggested scope check before drafting

Everything needed for a Methods and Results section describing the DESeq2 → feature-selection → RF/NN classification → SHAP → GO-enrichment pipeline is in `publication_v1/` and documented above. Before drafting, note that this project's scope does **not** currently include: clinical/cohort characteristics (§6.7), a validation cohort or external test set (only a single 80/20 split of the same TCGA samples was used), or any biological interpretation of the GO/KEGG enrichment terms (§5.4) — the term tables and plots exist, but deciding which are worth discussing is an open task for whoever drafts the Discussion.
