# publication_v1 — Corrected LUAD DEG/ML Re-run (padj-based)

A parallel, corrected re-run of the original pipeline (`scripts/`, `data/`,
`results/` at the project root — untouched by this run). This folder is
self-contained: scripts, data, results, figures, and environment records
specific to this corrected run all live under `publication_v1/`.

## What changed vs. the original run

The original pipeline (`scripts/01_diff_exp_analysis.R` and
`scripts/02_data_process.ipynb`) filtered on DESeq2's **raw p-value**
(`pvalue <= 0.05`) at two separate steps: per-sample DEG calling, and the
90%-of-comparisons feature-selection threshold. With **60,660 genes** tested
per comparison, this is a multiple-testing problem — at a raw p ≤ 0.05
threshold, roughly 3,000 genes are expected to look "significant" by chance
alone in any single comparison, before any real biological signal is
considered. This inflates the apparent significant-gene count and undermines
any downstream claim about "significant" biomarkers.

This run instead filters on **`padj`** — DESeq2's Benjamini-Hochberg (BH)
FDR-adjusted p-value — at both of those steps. The DESeq2 design (`~Group`,
`ref = "Control"`), the per-tumor-sample-vs-pooled-control comparison
structure, the `log2FoldChange >= 1` / `<= -1` up/down thresholds, the
90%-of-comparisons feature-selection logic, and the RF/NN modeling
approach are all unchanged from the original. The **SHAP marker
selection logic was changed** partway through this run, after the
original (unmodified) logic was found to produce zero markers on the
smaller corrected gene pool — see "SHAP marker selection" section below
for the full investigation and the fix that was applied, with your
explicit sign-off.

Additionally, this run adds things the original never computed:
performance metrics for both models, NN training history, a fixed
TensorFlow seed, and four DESeq2 QC diagnostic figures.

## Final sample sizes

- **N normal (Control) = 59**
- **N tumor (Treated) = 539**
- **N total = 598**

(Identical to the original run — sample selection was not touched.)

## Final gene count after feature selection — compared to the original

| | Threshold basis | Genes selected (padj/pvalue ≤ 0.05 in ≥ 90% of 539 comparisons) |
|---|---|---|
| Original run (`data/FEATURE_DATA_MATRIX.csv`) | raw p-value | **309** |
| This run (`publication_v1/data/FEATURE_DATA_MATRIX.csv`) | BH-adjusted p-value (padj) | **110** |

The gene count dropped by roughly two-thirds once multiple-testing
correction was applied — consistent with the raw-p-value run having
included a substantial number of false positives.

## Model performance (now present — the original never computed this)

`classification_report`/`accuracy_score` were imported in the original
`03_biomarker_shap.py` but never called, and the Keras model was never
`.evaluate()`d. Both are now computed on the held-out 20% test split
(120 samples: 16 normal / 104 tumor, `random_state=42`):

**Random Forest** (`publication_v1/results/rf_metrics.txt`): accuracy
**0.975**. Per-class: normal (class 0) precision 0.88 / recall 0.94 / F1
0.91 (n=16); tumor (class 1) precision 0.99 / recall 0.98 / F1 0.99 (n=104).
`RandomForestClassifier` is deterministic on CPU, so this result is
identical across every run of this pipeline.

**Keras NN** (`publication_v1/results/nn_metrics.txt`, final values):
test loss **0.0957**, test accuracy **0.975**. Full 50-epoch
training/validation loss and accuracy history is in
`publication_v1/results/nn_training_history.csv`.

**Reproducibility caveat, observed directly in this run:** an earlier
run of this exact script (same code, same `tf.random.set_seed(42)`,
same data) produced test loss 0.0743 / test accuracy 0.9667 — different
from the final numbers above. `RandomForestClassifier` reproduced
identically both times; the Keras model did not.
`tf.random.set_seed()` seeds TensorFlow's global RNG, but does not by
itself force bit-exact determinism for all ops on CPU (thread-level
floating-point reduction order can still vary between runs). The two
runs' numbers are close (96.7–97.5% accuracy) but should be treated as
"reproducible to within about a percentage point," not bit-identical,
unless further determinism settings (e.g.
`tf.config.experimental.enable_op_determinism()`, single-threaded
execution) are added — not done here, since that wasn't requested and
changes runtime behavior beyond seeding.

These numbers should be read cautiously: the normal class has only 16
test samples (of 120), so its precision/recall estimates carry wide
uncertainty despite the high headline accuracy. Accuracy alone is not a
reliable summary given the ~9:1 tumor:normal class imbalance.

## SHAP marker selection — bug found and fixed (with your sign-off)

The original selection logic (kept identical to `scripts/03_biomarker_shap.py`
in the first version of this run) builds two marker lists by walking two
SHAP-score-sorted gene rankings — one ranked by mean **|SHAP|** (absolute
value) on tumor test samples, one on normal test samples — in lockstep;
whenever a gene appears near the top of **both** rankings, it gets added
to one list then immediately evicted when encountered in the other. On
the padj-corrected **110-gene pool**, the two rankings turned out to be
highly correlated (Pearson r = 0.87 across all 110 genes; 70–88% top-K
overlap at every K), which made sense in hindsight: both rankings come
from explaining the *same* binary classifier, just evaluated on
differently-labeled subsets of test points. This caused the add-then-evict
cycle to empty **both lists down to zero** by the time all 110 gene pairs
were exhausted (traced step-by-step; the lists actually peaked at 12/12
around the loop's midpoint before the tail-end churn erased them) —
`SHAP_BIOMARKERS` had 0 rows and both SHAP figures were blank in the
first version of this run.

**Fix applied** (`scripts/03_biomarker_shap_padj.py`, "BIOMARKER SELECTION
LOGIC" section): instead of ranking by SHAP magnitude and then forcing
mutual exclusivity between the two lists, each gene's **signed** mean
SHAP contribution is used directly to ask the more direct question: does
this gene push the model's prediction toward tumor or toward normal?

- `TUMOR_MARKER` = genes with positive mean SHAP value among
  tumor-labeled test samples (i.e. push predictions toward tumor),
  ranked by strongest positive push, top 20.
- `NORMAL_MARKER` = genes with negative mean SHAP value among
  normal-labeled test samples (push predictions toward normal), ranked
  by strongest negative push, top 20.
- No forced exclusivity: a gene can legitimately appear in both lists
  (e.g. `ENSG00000234281.5`, `ENSG00000224215.1` do here) if it's a
  consistently strong discriminator in both directions. That is a
  meaningful result, not an error.

This run of the fix produced **20/20 markers** as intended, and
`TUMOR_FEATURE_IMP_NN_SELECTED.png` / `NORMAL_FEATURE_IMP_NN_SELECTED.png`
now show proper SHAP beeswarm plots (tumor markers cluster positive,
normal markers cluster negative, as expected by construction).

**This is a methods change from the original script**, made only after
you explicitly reviewed the investigation and approved this specific
approach — flagged here so it's clearly documented as a deliberate,
authorized deviation, not a silent correction.

## GO/pathway enrichment (added after the SHAP fix above)

Run on the fixed SHAP marker list (20 tumor + 20 normal genes,
`publication_v1/results/ML_BIOMARKERS_RES.xlsx`, `SHAP_BIOMARKERS` sheet)
via `scripts/04_pathway_go_enrichment.R`, adapted from
`legacy/step2_2025-02/SCRIPTS/BIOMARKER_PATH_GO_SCRIPT.R` — same analysis
logic (biomaRt Ensembl→HGNC mapping, same four Enrichr databases
`KEGG_2021_Human`, `GO_Biological_Process_2023`,
`GO_Molecular_Function_2023`, `GO_Cellular_Component_2023`, same
top-20-by-p-value dot/bar plots), only input/output paths changed to
point at this run's marker list instead of the old superseded one.
(The legacy script's `enrichplot` library import was unused in its own
code and was dropped here; a defensive filter removing blank HGNC-symbol
mappings before querying Enrichr was added, since biomaRt can return an
empty string for genes with no current symbol — neither changes which
real genes are analyzed.)

**One infrastructure fix was required:** `useMart("ensembl", ...)`
(the legacy script's connection call) returned an HTTP 404 — Ensembl's
BioMart discovery endpoint has moved since the legacy script was written.
Replaced with `useEnsembl(biomart = "genes", dataset =
"hsapiens_gene_ensembl")`, the current recommended connection function,
which also auto-retries via an Asia mirror when the main Ensembl site is
unresponsive (observed during this run). This is purely a connection-layer
fix — the `getBM()` query (attributes, filters, values) is unchanged.

**Gene mapping:** 20/20 tumor markers and 20/20 normal markers mapped to
valid HGNC gene symbols — no genes were dropped.

**Outputs:**
- `results/pathway_go/SHAP_TUMOR_PATHWAY_GO_RES.xlsx`,
  `SHAP_NORMAL_PATHWAY_GO_RES.xlsx` — one sheet per database (KEGG,
  BIOLOGICAL_PRO, MOLECULAR_FUN, CELLULAR_COMP) per marker set, full
  Enrichr result tables (all non-empty: e.g. 320 Biological Process terms
  for the tumor marker set, 144 for normal).
- `results/pathway_go/SHAP_Biomarkers_Gene_Symbols.xlsx` — the
  Ensembl ID → HGNC symbol mapping table for both marker sets.
- `figures/pathway_go/` — 16 PNGs (dot + bar plot × 4 databases × 2
  marker sets), each showing the top 20 terms by raw P-value.

**Cosmetic labeling issue inherited from the legacy script (not
introduced here, not fixed without your sign-off):** the dot plots'
color legend is labeled "p.adj" but the plotted/colored value is
actually the raw Enrichr `P.value` column, not a separately-adjusted
p-value — the plotting code never references an adjusted-p column at
all. This is a label-only issue; it doesn't affect the underlying data
in the xlsx files, which contain Enrichr's own `P.value` and
`Adjusted.P.value` columns both, correctly.

These enrichment results are reported as generated, without biological
interpretation — deciding which terms are meaningful for the preprint is
your call to make from the full tables and plots.

## Where each file came from

| File | Produced by |
|---|---|
| `data/LUAD_MERGED_ALL_DEG_RESULT_padj.csv` | `scripts/01_diff_exp_analysis_padj.R` |
| `data/DEG_RES_ALL/*.xlsx` (539 files) | `scripts/01_diff_exp_analysis_padj.R`, one per tumor sample |
| `data/FEATURE_DATA_MATRIX.csv` | `scripts/02_data_process_padj.py` |
| `results/ML_BIOMARKERS_RES.xlsx` | `scripts/03_biomarker_shap_padj.py` |
| `results/rf_metrics.txt` | `scripts/03_biomarker_shap_padj.py` (RF `classification_report`/`accuracy_score`) |
| `results/nn_metrics.txt` | `scripts/03_biomarker_shap_padj.py` (Keras `model.evaluate()`) |
| `results/nn_training_history.csv` | `scripts/03_biomarker_shap_padj.py` (Keras `history.history`) |
| `figures/Top_MARKERS_RandomForest.png` | `scripts/03_biomarker_shap_padj.py` |
| `figures/TUMOR_FEATURE_IMP_NN_SELECTED.png`, `NORMAL_FEATURE_IMP_NN_SELECTED.png` | `scripts/03_biomarker_shap_padj.py` — SHAP beeswarm plots for the 20 signed-SHAP markers per class, see fix described above |
| `figures/PCA_plot.png`, `MA_plot.png`, `volcano_plot.png`, `dispersion_plot.png` | `scripts/05_qc_diagnostics.R` — a **single combined** DESeq2 model across all 598 samples (design `~Group`, ref="Control"), separate from the 539 pairwise per-sample comparisons used for the actual DEG calls. This is the standard way to produce whole-cohort QC figures and does not replace the per-sample DEG methodology. |
| `results/pathway_go/*.xlsx` (3 files), `figures/pathway_go/*.png` (16 files) | `scripts/04_pathway_go_enrichment.R` — GO/KEGG enrichment on the SHAP marker list, see section above |
| `environment/r_sessionInfo.txt` | `sessionInfo()`, captured at run time |
| `environment/python_pip_freeze.txt` | `pip freeze` inside `BIOENV` (the project's isolated Python 3.10 venv) |
| `environment/stage1_run.log`, `stage2_run.log`, `stage3_run.log`, `stage5_qc_run.log` | raw console output of each script, for provenance |

Note on `MA_plot.png`: `plotMA()`'s default blue/"significant" coloring
uses DESeq2's own default `alpha = 0.1`, not this study's actual 0.05
padj threshold. This is cosmetic only — it does not affect the volcano
plot (which explicitly uses padj ≤ 0.05) or any filtering/gene selection
elsewhere in the pipeline.

`scripts/02_data_process_padj.py` re-implements the original
`scripts/02_data_process.ipynb` as a plain script (not a notebook) so the
corrected run could be executed non-interactively and logged; the
filtering logic itself is unchanged apart from padj vs. raw p-value.

## Reproducibility additions

- `tf.random.set_seed(42)` added in `03_biomarker_shap_padj.py`, on top of
  the original's existing `train_test_split(random_state=42)` and
  `RandomForestClassifier(random_state=52)`. This does not fully
  guarantee bit-identical re-runs (some ops remain non-deterministic on
  CPU), but substantially reduces run-to-run variance in the NN/SHAP
  results, which was previously unseeded entirely.
- Stage 1 (539 independent per-tumor-sample DESeq2 calls) was
  parallelized with `parallel::mclapply(mc.cores = 4)`, matching this
  machine's Performance-core count (Apple M4: 4 P-cores + 6 slower
  E-cores). This is a purely computational change — mc.cores=4 was
  empirically the fastest setting tested (beating both mc.cores=6 and
  mc.cores=10, which suffer from P/E-core scheduling contention when
  oversubscribed) — and does not alter any DESeq2 inputs, parameters, or
  outputs; results are identical to what a serial loop would produce.
  Total stage 1 runtime: ~103 minutes for all 539 samples.

## What is NOT included in this run (deferred, not blocking)

- **Clinical/sample metadata** (age, sex, stage, batch, etc.) — not
  available anywhere in the project as a saved table; only `Group`
  (Control/Treated) exists. Deferred per your instruction; needs a
  separate data-sourcing step before it can be added to any manuscript
  table.
- GO/pathway enrichment is now done — see the new section below. It was
  run on the SHAP marker list only, not the RF ranking (RF enrichment
  was already commented out / never run in the legacy script either).

## Pipeline summary (plain English)

This pipeline starts from raw RNA-seq gene counts for 59 normal lung
tissue samples and 539 lung adenocarcinoma (LUAD) tumor samples from
TCGA. Step 1 uses DESeq2 to statistically compare each tumor sample
individually against the pooled group of 59 normal samples, testing all
60,660 genes each time and correcting for multiple testing (BH/padj) to
control the false-discovery rate — this produces a fold-change and
adjusted p-value for every gene, for every one of the 539 comparisons.
Step 2 keeps only the genes that were consistently significant (padj ≤
0.05) in at least 90% of those 539 comparisons — 110 genes survive this
filter — and builds a table of all 598 samples by these 110 genes. Step 3
uses that table to train two independent classifiers (a Random Forest and
a small neural network) to distinguish tumor from normal samples purely
from gene expression, both of which perform well on held-out test data
(97.5% accuracy for both models). It then uses SHAP, a model
explainability technique, to identify which of the 110 genes most
strongly push each model's predictions toward "tumor" versus toward
"normal" — producing a shortlist of 20 tumor-associated and 20
normal-associated genes (a gene can appear in both lists if it's a
strong discriminator in both directions). Step 4 takes those two
20-gene shortlists and asks, for each one, which broader biological
pathways and functions (KEGG pathways, GO Biological Process/Molecular
Function/Cellular Component terms) they're statistically enriched for —
producing ranked term lists and plots for each marker set. Four
additional diagnostic plots (PCA, MA, volcano, dispersion) provide
standard quality-control views of the overall tumor/normal separation
and the statistical behavior of the DESeq2 model.
