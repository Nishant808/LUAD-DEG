# Gene Expression Profiling and Pathway Enrichment Analysis in Lung Adenocarcinoma  
*A Machine Learning Approach*  
 

---

## 📘 Overview

This project explores the gene expression landscape of Lung Adenocarcinoma (LUAD) using machine learning (ML) techniques and pathway enrichment analysis. It integrates neural network-based classification, SHAP-based feature importance scoring, and downstream functional annotation to identify key biomarkers and dysregulated biological pathways relevant to LUAD tumorigenesis.

---

## 🎯 Objectives

- Classify LUAD and normal samples based on gene expression profiles using deep learning.
- Identify potential tumor-specific and normal-specific biomarkers using SHAP values.
- Perform Gene Ontology (GO) and KEGG pathway enrichment analysis to interpret the functional roles of identified biomarkers.
- Validate the robustness of biomarkers across independent datasets.

---

## 🧪 Methods

### 🔬 Data Collection & Preprocessing
- **Source:** TCGA (LUAD dataset)
- **Preprocessing:**
  - Label encoding (tumor vs. normal)
  - Normalization using `StandardScaler`
  - Splitting into training and testing sets

### 🤖 Machine Learning
- **Neural Network** using **TensorFlow (Keras Sequential API)**:
  - 3 hidden layers (Dropout = 0.6)
  - Binary classification (LUAD vs. Normal)
  - Trained for 50 epochs with `Adam` optimizer

### 🧠 SHAP Feature Importance
- Used **SHAP DeepExplainer** for model interpretation
- Identified top predictive genes as candidate biomarkers

### 🧬 Pathway Enrichment
- **Tools:** `biomaRt`, `EnrichR`, `ShinyGO`
- **Analysis:** GO Biological Process and KEGG Pathways
- **Gene ID conversion:** Ensembl to HGNC symbols
- **Visualizations:** Dot plots & bar plots highlighting enriched pathways

---

## 📊 Results

### 🔑 Key Findings
- **Tumor Enriched Pathways:**  
  - Lipid metabolism (e.g., cholesterol biosynthesis, PPAR signaling)  
  - Immune pathways (e.g., Phagosome, Hematopoietic cell lineage)  
  - AMPK signaling and insulin resistance
![Tumor_GO](https://github.com/Nishant808/LUAD-DEG/blob/main/GO_BP_Tumor_BARPLOT.png)

- **Top Biomarkers:**  
  - Identified via SHAP value rankings  
  - Categorized into tumor and normal-specific sets  
  - Saved in `ML_BIOMARKERS_RES.xlsx`
![SHAP_Tumor](https://github.com/Nishant808/LUAD-DEG/blob/main/TUMOR_FEATURE_IMP_NN_SELECTED.png)

- **Statistical Significance:**  
  - High fold enrichment (up to 400×) in lipid pathways  
  - FDR-adjusted p-values confirming pathway relevance
![FDR_DOTPLOT](https://github.com/Nishant808/LUAD-DEG/blob/main/Tumor_KEGG_DOTPLOT.png)
---

## 🧰 Tools & Technologies

- **Languages:** Python, R
- **Libraries:** TensorFlow, SHAP, NumPy, Pandas, ggplot2
- **Bioinformatics Databases:** TCGA, Ensembl, EnrichR, KEGG, GEO
- **Gene Annotation:** `biomaRt` for Ensembl → HGNC mapping
- **Visualization:** SHAP summary plots, ggplot2 (bar/dot plots)

---
