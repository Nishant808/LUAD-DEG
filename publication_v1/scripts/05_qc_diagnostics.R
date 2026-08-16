
# QC diagnostics for the LUAD DESeq2 analysis: PCA, MA plot, volcano plot
# (padj-based), and dispersion plot.
#
# These are computed from ONE combined DESeqDataSet across all 598 samples
# (design = ~Group, ref = "Control") -- the standard way to produce
# whole-cohort QC figures. This is separate from, and does not replace,
# the per-tumor-sample-vs-pooled-control DESeq2 structure used in
# 01_diff_exp_analysis_padj.R for the actual DEG calls.

library(DESeq2)
library(ggplot2)

PROJECT_DIR = "/Users/nishantthalwal/Documents/Data ALL/TCGA_LUAD"
COUNTS_PATH = paste0(PROJECT_DIR, "/data/counts")
FIG_PATH = paste0(PROJECT_DIR, "/publication_v1/figures")
dir.create(FIG_PATH, recursive = TRUE, showWarnings = FALSE)

Normal_Counts = read.table(paste0(COUNTS_PATH, "/NORMAL/Normal_LUAD_Counts.tsv"), sep = "\t")
Tumor_Counts  = read.table(paste0(COUNTS_PATH, "/TUMOR/Tumor_LUAD_Counts.tsv"), sep = "\t")

Count_Mat = cbind(Normal_Counts, Tumor_Counts)
Metadata = data.frame(
  Sample_ID = colnames(Count_Mat),
  Group = c(rep("Control", ncol(Normal_Counts)), rep("Treated", ncol(Tumor_Counts)))
)
rownames(Metadata) = Metadata$Sample_ID

dds <- DESeqDataSetFromMatrix(countData = Count_Mat, colData = Metadata, design = ~Group)
dds$Group <- relevel(dds$Group, ref = "Control")
dds <- DESeq(dds)
res <- results(dds)

cat("Combined model: N samples =", ncol(dds), " N genes =", nrow(dds), "\n")

# --- PCA plot ---
vsd <- vst(dds, blind = TRUE)
png(file.path(FIG_PATH, "PCA_plot.png"), width = 8, height = 6, units = "in", res = 300)
print(plotPCA(vsd, intgroup = "Group") + ggtitle("PCA - LUAD Tumor vs Normal (all samples)"))
dev.off()

# --- MA plot ---
png(file.path(FIG_PATH, "MA_plot.png"), width = 8, height = 6, units = "in", res = 300)
plotMA(res, main = "MA Plot - LUAD Tumor vs Normal", ylim = c(-8, 8))
dev.off()

# --- Volcano plot (padj-based) ---
res_df <- as.data.frame(res)
res_df$Significant <- with(res_df, !is.na(padj) & padj <= 0.05 & abs(log2FoldChange) >= 1)
volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot - LUAD Tumor vs Normal (padj)",
       x = "log2 Fold Change", y = "-log10(padj)")
ggsave(file.path(FIG_PATH, "volcano_plot.png"), volcano, width = 8, height = 6, dpi = 300)

# --- Dispersion plot ---
png(file.path(FIG_PATH, "dispersion_plot.png"), width = 8, height = 6, units = "in", res = 300)
plotDispEsts(dds, main = "DESeq2 Dispersion Estimates - LUAD Tumor vs Normal")
dev.off()

cat("QC diagnostics complete.\n")
