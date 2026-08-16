
library(biomaRt)
library(readxl)
library(writexl)
library(enrichR)
library(enrichplot)
library(ggplot2)
library(stringr)

PATH = "BIOMARKERS_RES/ML_BIOMARKERS_RES.xlsx"

setwd(dirname(PATH))

# Connect to the Ensembl database

# RF_DF = data.frame(read_excel("ML_BIOMARKERS_RES.xlsx", sheet = "RANDOM_FOREST_BIOMARKERS"))
# RF_TOP = head(RF_DF[order(RF_DF$Importance, decreasing = TRUE), 1], 20)

SHAP_DF  = data.frame(read_excel("ML_BIOMARKERS_RES.xlsx", sheet = "SHAP_BIOMARKERS"))
SHAP_TUMOR = SHAP_DF$TUMOR_MARKER
SHAP_NORMAL = SHAP_DF$NORMAL_MARKER

dir.create("PATHWAY_GO")
dir.create("PATHWAY_GO/PLOTS")

PLOT_PATH_ENRICH = function(ANN_DF, DB, MODEL){
  TOP_ANN = head(ANN_DF[order(ANN_DF$P.value, decreasing = FALSE), ], 20)
  TOP_ANN$GENE_COUNT <- sapply(strsplit(TOP_ANN$Overlap, "/"), function(x) x[1])
  TOP_ANN$Term <- str_wrap(TOP_ANN$Term, width = 100)
  TOP_ANN$Term <- factor(TOP_ANN$Term, levels = TOP_ANN$Term[order(TOP_ANN$P.value, decreasing = FALSE)])
  
  DOTPLOT = ggplot(TOP_ANN, aes(x = P.value,  y = Term, size = GENE_COUNT, color = P.value)) +
    geom_point() +
    labs(x = "P-Value", y = "Term", size = "Gene Count", color = "p.adj") +
    theme_minimal() +
    theme(axis.text.y = element_text(face = 'bold', size = 12))+
    scale_color_gradient(low = "blue", high = "red")
  png(paste0("PATHWAY_GO/PLOTS/",MODEL, "_", DB,"_DOT.png"), height = 6, width = 14, units = "in", res = 300)
  print(DOTPLOT)
  dev.off()
  
  BARPLOT = ggplot(TOP_ANN, aes(x = Term,  y = P.value, fill = P.value)) +
    geom_bar(stat="identity", width = 0.7) +
    labs(x = "Terms", y = "P-Value")+
    coord_flip()+
    theme_minimal()+
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(face = 'bold'),
          axis.text = element_text(face = 'bold', size = 12),
          # axis.line = element_line(color = "black", size = 0.3),
          axis.line = element_blank())
  png(paste0("PATHWAY_GO/PLOTS/", MODEL, "_", DB,"_BAR.png"), height = 6, width = 18, units = "in", res = 300)
  print(BARPLOT)
  dev.off()
  
}

ENRICH_PATH = function(ENS_ID_TOP, MODEL){
  ENS_ID = sapply(strsplit(ENS_ID_TOP, "\\."), function(ENS_ID_TOP) ENS_ID_TOP[1])
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # Retrieve gene symbols
  GENES_DF = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                   filters = 'ensembl_gene_id', 
                   values = ENS_ID, 
                   mart = ensembl)

  GENES = GENES_DF$hgnc_symbol
  setEnrichrSite("Enrichr")
  
  # DATABASES <- listEnrichrDbs()
  # sort(DATABASES$libraryName)
  
  BASES <- c("KEGG_2021_Human", "GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023")
  enriched <- enrichr(GENES, BASES)
  KEGG_DF = enriched$KEGG_2021_Human
  BP_DF = enriched$GO_Biological_Process_2023
  MF_DF = enriched$GO_Molecular_Function_2023
  CC_DF = enriched$GO_Cellular_Component_2023
  
  ANN_LIST = list(KEGG = KEGG_DF, BIOLOGICAL_PRO = BP_DF, MOLECULAR_FUN = MF_DF, CELLULAR_COMP = CC_DF)
  write_xlsx(ANN_LIST, paste0("PATHWAY_GO/", MODEL, "_PATHWAY_GO_RES.xlsx"))
  PLOT_PATH_ENRICH(KEGG_DF, "KEGG", MODEL)
  PLOT_PATH_ENRICH(BP_DF, "BP", MODEL)
  PLOT_PATH_ENRICH(MF_DF, "MF", MODEL)
  PLOT_PATH_ENRICH(CC_DF, "CC", MODEL)
  
  return(GENES_DF)
}

# ENRICH_PATH(RF_TOP, "RANDOM_FOREST")
TUMAR_SYMBOL_DF = ENRICH_PATH(SHAP_TUMOR, "SHAP_TUMOR")
NORMAL_SYMBOL_DF = ENRICH_PATH(SHAP_NORMAL, "SHAP_NORMAL")



RES_LIST = list(TumorFeatures = TUMAR_SYMBOL_DF, Normal_Features = NORMAL_SYMBOL_DF)
write_xlsx(RES_LIST, "SHAP_Biomarkers_Gene_Symbols.xlsx")


