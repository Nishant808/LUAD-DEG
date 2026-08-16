
# GO/KEGG pathway enrichment on the signed-SHAP tumor/normal marker lists
# (publication_v1/results/ML_BIOMARKERS_RES.xlsx, SHAP_BIOMARKERS sheet).
#
# Adapted from legacy/step2_2025-02/SCRIPTS/BIOMARKER_PATH_GO_SCRIPT.R --
# same analysis logic (biomaRt Ensembl ID -> HGNC symbol mapping, same
# four Enrichr databases, same top-20-by-p-value plotting), only the
# input/output paths changed to point at this run's (fixed, 20/20) SHAP
# marker list instead of the old superseded one. (legacy's `enrichplot`
# import was unused in the original script and is dropped here.)
#
# Requires internet access: biomaRt queries Ensembl's public BioMart
# service, enrichR queries the Enrichr web API (https://maayanlab.cloud).

library(readxl)
library(writexl)
library(enrichR)
library(ggplot2)
library(stringr)
library(biomaRt)

PROJECT_DIR = "/Users/nishantthalwal/Documents/Data ALL/TCGA_LUAD"
MARKERS_XLSX = paste0(PROJECT_DIR, "/publication_v1/results/ML_BIOMARKERS_RES.xlsx")
RESULTS_DIR = paste0(PROJECT_DIR, "/publication_v1/results/pathway_go")
FIG_DIR = paste0(PROJECT_DIR, "/publication_v1/figures/pathway_go")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

SHAP_DF = data.frame(read_excel(MARKERS_XLSX, sheet = "SHAP_BIOMARKERS"))
SHAP_TUMOR = SHAP_DF$TUMOR_MARKER
SHAP_NORMAL = SHAP_DF$NORMAL_MARKER
cat("SHAP_TUMOR genes:", length(SHAP_TUMOR), "\n")
cat("SHAP_NORMAL genes:", length(SHAP_NORMAL), "\n")

PLOT_PATH_ENRICH = function(ANN_DF, DB, MODEL){
  if (is.null(ANN_DF) || nrow(ANN_DF) == 0) {
    cat("No enrichment results for", MODEL, DB, "- skipping plots.\n")
    return(invisible(NULL))
  }
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
  png(file.path(FIG_DIR, paste0(MODEL, "_", DB, "_DOT.png")), height = 6, width = 14, units = "in", res = 300)
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
          axis.line = element_blank())
  png(file.path(FIG_DIR, paste0(MODEL, "_", DB, "_BAR.png")), height = 6, width = 18, units = "in", res = 300)
  print(BARPLOT)
  dev.off()
}

ENRICH_PATH = function(ENS_ID_TOP, MODEL){
  ENS_ID = sapply(strsplit(ENS_ID_TOP, "\\."), function(x) x[1])
  # useMart("ensembl", ...) returned HTTP 404 (Ensembl's biomart discovery
  # endpoint has moved); useEnsembl() is the current recommended connection
  # function and auto-falls back to a mirror if the main site is unresponsive.
  # Same underlying data/query -- getBM() below is unchanged.
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  GENES_DF = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = ENS_ID,
                   mart = ensembl)

  GENES = GENES_DF$hgnc_symbol
  GENES = GENES[GENES != ""]
  cat(MODEL, "- Ensembl IDs:", length(ENS_ID), " mapped HGNC symbols:", length(GENES), "\n")

  setEnrichrSite("Enrichr")
  BASES <- c("KEGG_2021_Human", "GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023")
  enriched <- enrichr(GENES, BASES)
  KEGG_DF = enriched$KEGG_2021_Human
  BP_DF = enriched$GO_Biological_Process_2023
  MF_DF = enriched$GO_Molecular_Function_2023
  CC_DF = enriched$GO_Cellular_Component_2023

  ANN_LIST = list(KEGG = KEGG_DF, BIOLOGICAL_PRO = BP_DF, MOLECULAR_FUN = MF_DF, CELLULAR_COMP = CC_DF)
  write_xlsx(ANN_LIST, file.path(RESULTS_DIR, paste0(MODEL, "_PATHWAY_GO_RES.xlsx")))
  PLOT_PATH_ENRICH(KEGG_DF, "KEGG", MODEL)
  PLOT_PATH_ENRICH(BP_DF, "BP", MODEL)
  PLOT_PATH_ENRICH(MF_DF, "MF", MODEL)
  PLOT_PATH_ENRICH(CC_DF, "CC", MODEL)

  return(GENES_DF)
}

TUMOR_SYMBOL_DF = ENRICH_PATH(SHAP_TUMOR, "SHAP_TUMOR")
NORMAL_SYMBOL_DF = ENRICH_PATH(SHAP_NORMAL, "SHAP_NORMAL")

RES_LIST = list(TumorFeatures = TUMOR_SYMBOL_DF, Normal_Features = NORMAL_SYMBOL_DF)
write_xlsx(RES_LIST, file.path(RESULTS_DIR, "SHAP_Biomarkers_Gene_Symbols.xlsx"))

cat("GO/pathway enrichment complete.\n")
