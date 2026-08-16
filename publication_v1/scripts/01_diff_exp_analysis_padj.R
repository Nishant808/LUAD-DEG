
# Corrected re-run of scripts/01_diff_exp_analysis.R
# CHANGE vs original: filters on DESeq2's BH-adjusted p-value (padj), not
# raw pvalue, at the per-sample significance step. Everything else
# (design, per-sample-vs-pooled-control structure, log2FC thresholds)
# is unchanged.
#
# Parallelized across the 539 independent per-tumor-sample DESeq2 calls
# using parallel::mclapply. mc.cores = 4, matching this machine's
# Performance-core count (Apple M4: 4 P-cores + 6 E-cores) -- empirically
# the fastest setting (154.6s/20 samples), beating mc.cores=6 (202.6s) and
# mc.cores=10 (174.7s), which suffer from P/E-core scheduling contention
# when oversubscribed.

library(DESeq2)
library(tibble)
library(dplyr)
library(purrr)
library(writexl)
library(parallel)

PROJECT_DIR = "/Users/nishantthalwal/Documents/Data ALL/TCGA_LUAD"
COUNTS_PATH = paste0(PROJECT_DIR, "/data/counts")
OUT_DATA_PATH = paste0(PROJECT_DIR, "/publication_v1/data")
N_CORES = 4

Normal_PATH = paste0(COUNTS_PATH, "/NORMAL/Normal_LUAD_Counts.tsv")
Tumor_PATH = paste0(COUNTS_PATH, "/TUMOR/Tumor_LUAD_Counts.tsv")
print(Normal_PATH)
print(Tumor_PATH)

Normal_Counts = read.table(Normal_PATH, sep = "\t")
Tumor_Counts = read.table(Tumor_PATH, sep = "\t")

Metadata_norm = data.frame(Sample_ID = colnames(Normal_Counts), Group = "Control")
Metadata_Tumor = data.frame(Sample_ID = colnames(Tumor_Counts), Group = "Treated")
Metadata = rbind(Metadata_norm, Metadata_Tumor)
rownames(Metadata) = Metadata$Sample_ID

Normal_Counts_rc = rownames_to_column(Normal_Counts, var = "Sample")

dir.create(OUT_DATA_PATH, recursive = TRUE, showWarnings = FALSE)
setwd(OUT_DATA_PATH)
dir.create("DEG_RES_ALL", showWarnings = FALSE)

process_sample = function(Sample){
  TUM_SAM = data.frame(Sample = row.names(Tumor_Counts), Count = Tumor_Counts[ , Sample])
  colnames(TUM_SAM) = NULL
  colnames(TUM_SAM) = c("Sample", Sample)
  COUNT_DF = merge(Normal_Counts_rc, TUM_SAM, by = "Sample")
  COUNT_DF  = column_to_rownames(COUNT_DF, var = "Sample")
  METADATA = Metadata[colnames(COUNT_DF), ]
  dds <- DESeqDataSetFromMatrix(countData=COUNT_DF,
                                colData=METADATA,
                                design=~Group)
  dds$Group <- relevel(dds$Group, ref = "Control")
  dds <- DESeq(dds)
  res = data.frame(results(dds))
  res = rownames_to_column(res, var = "ID")
  # CHANGE: padj (BH-adjusted) instead of raw pvalue
  sig = subset(res, res$padj <= 0.05)
  up = subset(sig, sig$log2FoldChange >= 1)
  down = subset(sig, sig$log2FoldChange <= -1)
  DGE_RES_LIST = list(ALL_DEG = res, SIGNIFICANT = sig, UPREGULATED = up, DOWNREGULATED = down)
  write_xlsx(DGE_RES_LIST, file.path("DEG_RES_ALL", paste0(Sample, "_vs_Control_DEG_RESULT.xlsx")))
  # CHANGE: carry padj (not pvalue) into the merged per-sample table
  RES = res[ , c("ID", "log2FoldChange", "padj")]
  colnames(RES) = NULL
  colnames(RES) = c("ID",  paste0("Control_vs_",Sample, "_LOG2FC"), paste0("Control_vs_",Sample, "_PADJ"))
  RES
}

t0 = Sys.time()
RES_LIST = mclapply(colnames(Tumor_Counts), process_sample, mc.cores = N_CORES)
t1 = Sys.time()
cat("DESeq2 loop across", length(RES_LIST), "samples took:", as.numeric(t1 - t0, units = "mins"), "min\n")

errors = Filter(function(x) inherits(x, "try-error"), RES_LIST)
if (length(errors) > 0) {
  stop(paste("mclapply worker failures:", length(errors), "- aborting before merge."))
}

Merged_DF = reduce(RES_LIST, inner_join, by = "ID")
write.csv(Merged_DF, "LUAD_MERGED_ALL_DEG_RESULT_padj.csv", row.names = FALSE)
cat("Merged output written:", nrow(Merged_DF), "genes x", ncol(Merged_DF), "columns\n")
