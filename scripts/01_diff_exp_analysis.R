
library(DESeq2)
library(tibble)
library(dplyr)
library(purrr)
library(writexl)


PROJECT_DIR = "/Users/nishantthalwal/Documents/Data ALL/TCGA_LUAD"
PATH = paste0(PROJECT_DIR, "/data/counts")
DATA_PATH = paste0(PROJECT_DIR, "/data")


DEG_CALCULATION = function(Normal_PATH, Tumor_PATH, cancer_type){
  Normal_Counts = read.table(Normal_PATH, sep = "\t")
  Tumor_Counts = read.table(Tumor_PATH, sep = "\t")
  Metadata_norm = data.frame(Sample_ID = colnames(Normal_Counts),
                             Group = "Control")
  
  Metadata_Tumor = data.frame(Sample_ID = colnames(Tumor_Counts),
                              Group = "Treated")
  Metadata = rbind(Metadata_norm, Metadata_Tumor)
  rownames(Metadata) = Metadata$Sample_ID
  
  #------------------------------------------------------------------------------------------
  
  # Count_Mat = cbind(Normal_Counts, Tumor_Counts)
  # 
  # dds <- DESeqDataSetFromMatrix(countData=Count_Mat,
  #                               colData=Metadata,
  #                               design=~Group)
  # 
  # dds$Group <- relevel(dds$Group, ref = "Control")
  # dds <- DESeq(dds)
  # res = data.frame(results(dds))
  # res = rownames_to_column(res, var = "ID")
  # sig = subset(res, res$pvalue <= 0.05)
  # up = subset(sig, sig$log2FoldChange >= 1)
  # down = subset(sig, sig$log2FoldChange <= -1)
  # DGE_RES_LIST = list(ALL_DEG = res, SIGNIFICANT = sig, UPREGULATED = up, DOWNREGULATED = down)
  # write_xlsx(DGE_RES_LIST, paste0(cancer_type, "_Case_vs_Control_DEG_RESULT_ALL.xlsx"))
  # write.csv(res, paste0(cancer_type, "_Case_vs_Control_DEG_RESULT_ALL.csv"))
  #------------------------------------------------------------------------------------------
  Normal_Counts  = rownames_to_column(Normal_Counts, var = "Sample")

  dir.create("DEG_RES_ALL")
  setwd("DEG_RES_ALL")
  
  RES_LIST = list()
  for (Sample in colnames(Tumor_Counts)){
    TUM_SAM = data.frame(Sample = row.names(Tumor_Counts), Count = Tumor_Counts[ , Sample])
    colnames(TUM_SAM) = NULL
    colnames(TUM_SAM) = c("Sample", Sample)
    COUNT_DF = merge(Normal_Counts, TUM_SAM, by = "Sample")
    COUNT_DF  = column_to_rownames(COUNT_DF, var = "Sample")
    METADATA = Metadata[colnames(COUNT_DF), ]
    dds <- DESeqDataSetFromMatrix(countData=COUNT_DF,
                                  colData=METADATA,
                                  design=~Group)
    dds$Group <- relevel(dds$Group, ref = "Control")
    dds <- DESeq(dds)
    res = data.frame(results(dds))
    res = rownames_to_column(res, var = "ID")
    sig = subset(res, res$pvalue <= 0.05)
    up = subset(sig, sig$log2FoldChange >= 1)
    down = subset(sig, sig$log2FoldChange <= -1)
    DGE_RES_LIST = list(ALL_DEG = res, SIGNIFICANT = sig, UPREGULATED = up, DOWNREGULATED = down)
    RES = res[ , c("ID", "log2FoldChange", "pvalue")]
    colnames(RES) = NULL
    colnames(RES) = c("ID",  paste0("Control_vs_",Sample, "_LOG2FC"), paste0("Control_vs_",Sample, "_P_VALUE"))
    RES_LIST = c(RES_LIST, list(RES))
    write_xlsx(DGE_RES_LIST, paste0(Sample, "_vs_Control_DEG_RESULT.xlsx"))
    # write.csv(res, paste0(Sample, "_vs_Control_DEG_RESULT.csv"), row.names = FALSE)
  }
  Merged_DF = reduce(RES_LIST,inner_join, by = "ID")
  setwd("..")
  write.csv(Merged_DF, paste0(cancer_type,"_MERGED_ALL_DEG_RESULT.csv"), row.names = FALSE)
  
}

#-------------------------------------------------------------------------------

Normal_PATH = paste0(PATH,"/NORMAL/Normal_LUAD_Counts.tsv")
Tumor_PATH = paste0(PATH, "/TUMOR/Tumor_LUAD_Counts.tsv")

print(Normal_PATH)
print(Tumor_PATH)

setwd(DATA_PATH)
DEG_CALCULATION(Normal_PATH, Tumor_PATH, "LUAD")
