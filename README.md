# 🧬 ClinicalGenomicStats: R Toolkit for Clinical & Molecular Data Analysis

> A robust, end-to-end R pipeline designed for molecular biologists, geneticists, and clinical researchers. This repository provides automated workflows for demographic extraction, gene expression analysis (qRT-PCR), statistical testing, and generating high-resolution, publication-ready (300 DPI) visualizations.

![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)
![Data Analysis](https://img.shields.io/badge/Data%20Analysis-Bioinformatics-008080?style=for-the-badge)

---

## 📑 Table of Contents
1. [Setup & Installation](#1-environment-setup-01requiredpackagesr)
2. [Normality Testing](#2-normality-testing-02normalitytestsr)
3. [Expression RainCloud Plots](#3-expression-raincloud-plots-03raincloudplotr)
4. [Biomarker ROC Analysis](#4-biomarker-roc-analysis-04roc-analysisr)
5. [Demographic Extraction](#5-demographic-extraction-05demographicdataextractionr)
6. [Clinical Features Correlation](#6-clinical-features-correlation-06analysisofclinicalfeaturesr)
7. [qRT-PCR Expression Analysis](#7-qrt-pcr-expression-analysis-07qpcranalyzer)
8. [Cell Cycle Flow Cytometry](#8-cell-cycle-flow-cytometry-08cellcycleanalyzer)

---

## 🚀 Quick Start

To begin utilizing this repository, clone it to your local machine or Linux server:

```bash
git clone [https://github.com/YourUsername/ClinicalGenomicStats.git](https://github.com/YourUsername/ClinicalGenomicStats.git)
cd ClinicalGenomicStats
```

---

## 📦 1. Environment Setup (01.RequiredPackages.R)

Before running any analytical scripts, you must configure your R environment. This script intelligently scans your system, identifies missing dependencies, installs them seamlessly from CRAN, and loads them for immediate use.

**Dependencies handled:** `ggplot2`, `dplyr`, `pROC`, `readxl`, `ggdist`, `writexl`, `nortest`, `gridExtra`, `fs`, `openxlsx`, `pacman`, `tidyr`.

### 💻 Usage

```r
source("01.RequiredPackages.R")
# Automatically runs init_packages() and prepares your workspace.
```

---

## 📊 2. Normality Testing (02.NormalityTests.R)

Statistical validity begins with assessing data distribution. This script performs Shapiro-Wilk, Anderson-Darling, and Lilliefors tests on a specified clinical variable. It outputs a comprehensive statistical summary (Excel) alongside high-quality Histograms and Q-Q plots.

### 💻 Usage

```r
source("02.NormalityTests.R")

# Run normality tests on a specific column
results <- normality(
  excel_path = "clinical_data.xlsx", 
  column = "Age", 
  preferred_name = "Patient_Age"
)
```

---

## 🌧️ 3. Expression RainCloud Plots (03.RainCloudPlot.R)

Compare relative gene expression between tissues (e.g., Tumor vs. Marginal) using aesthetically pleasing RainCloud plots. These plots combine a half-violin distribution, boxplot, and raw data jitter for maximum transparency. Automatically computes non-parametric tests (Mann-Whitney or Wilcoxon).

### 💻 Usage

```r
source("03.RainCloudPlot.R")

generate_raincloud(
  file_path = "expression_data.xlsx",
  gene_name = "P53",
  col_case  = "TumoralRelativeExp",
  col_ctrl  = "MarginalRelativeExp",
  label_case = "Tumor",
  label_control = "Marginal",
  test_method = "wilcoxon" # Options: "mann_whitney" or "wilcoxon"
)
```

---

## 🎯 4. Biomarker ROC Analysis (04.ROC-Analysis.R)

Evaluate the diagnostic potential of a gene. This script calculates the Area Under the Curve (AUC), 95% Confidence Intervals, Sensitivity, Specificity, PPV, NPV, and Likelihood Ratios based on Youden's Index. ROC curves are exported in `.tiff`, `.png`, `.jpeg`, and `.pdf` formats.

### 💻 Usage

```r
source("04.ROC-Analysis.R")
library(readxl)
library(writexl)

my_data <- read_excel("expression_data.xlsx")

roc_results <- analyze_roc_gene(
  data = my_data,
  case_col = "TumoralRelativeExp",
  control_col = "MarginalRelativeExp",
  gene_name = "P53",
  roc_color = "#1976D2"
)

write_xlsx(roc_results, "P53_ROC_Stats.xlsx")
```

---

## 📋 5. Demographic Extraction (05.DemographicDataExtraction.R)

Instantly generate Table 1 for your clinical paper. This function intelligently detects numeric vs. categorical columns, outputting Mean ± SD (and Ranges) for continuous variables, and Frequencies with Percentages for categorical variables.

### 💻 Usage

```r
source("05.DemographicDataExtraction.R")
library(readxl)

clinical_df <- read_excel("clinical_data.xlsx")

# Ensure categorical data (like Stage/Grade) are formatted as factors
clinical_df$Stage <- as.factor(clinical_df$Stage)

my_vars <- c("Age", "Histology", "Stage", "HPV")

extract_demographics(
  dataset = clinical_df, 
  selected_columns = my_vars,
  output_name = "Demographic_Table1.xlsx"
)
```

---

## 🧬 6. Clinical Features Correlation (06.AnalysisOfClinicalFeatures.R)

Investigate how gene expression correlates with clinical pathological features. This manager function handles both discrete data (Kruskal-Wallis/Mann-Whitney -> Boxplots) and continuous data (Spearman correlation -> Scatter plots with trendlines).

### 💻 Usage

```r
source("06.AnalysisOfClinicalFeatures.R")
library(readxl)

my_data <- read_excel("clinical_and_expression.xlsx")
my_data$Histology <- as.factor(my_data$Histology)

genes_config <- list(
  list(name = "P53", case = "TumoralRelativeExp", control = "MarginalRelativeExp")
)

clinical_vars <- c("Age", "Histology", "Tumor-size", "Stage")
my_cutoffs <- list("Age" = 50) # Binarizes age at 50

run_multi_gene_analysis(
  dataset = my_data,
  gene_config_list = genes_config,
  clinical_cols = clinical_vars,
  cutoff_list = my_cutoffs
)
```

---

## 🔬 7. qRT-PCR Expression Analysis (07.qpcrAnalyze.R)

A complete automated pipeline for qRT-PCR data. Input your Delta Ct values, and the script evaluates normality, checks variance equality, runs the appropriate T-test (Student's or Welch's), computes Fold Change ($2^{-\Delta\Delta Ct}$), and generates labeled boxplots with standard significance stars (***).

### 💻 Usage

```r
source("07.qpcrAnalyze.R")
library(readxl)

qpcr_data <- read_excel("qpcr_raw.xlsx")

analyze_qpcr(
  input_data = qpcr_data,
  case_col = "DeltaCtCase",
  control_col = "DeltaCtControl",
  case_label = "si-P53",
  control_label = "Scramble",
  y_label = "Expression ratio P53/B2M",
  palette = c("#FF6347", "#4682B4")
)
```

---

## 🔄 8. Cell Cycle Flow Cytometry (08.CellCycleAnalyze.R)

Analyze phase distributions (G1, S, G2/M) from flow cytometry experiments. This function calculates standard deviations, performs independent T-tests, outputs a detailed statistical report, and generates a stunning combined multi-group plot with automated significance brackets.

### 💻 Usage

```r
source("08.CellCycleAnalyze.R")
library(readxl)

flow_data <- read_excel("flow_cytometry.xlsx")

# Map phases to their respective Control and Treatment columns
cycle_phases <- list(
  "G1"   = c("G1-Control", "G1-Si"),
  "S"    = c("S-Control", "S-Si"),
  "G2/M" = c("G2-Control", "G2-Si")
)

analyze_cell_cycle(
  input_data    = flow_data,
  phase_list    = cycle_phases,
  control_label = "Control",
  case_label    = "Si-Treated",
  palette       = c("#512DA8", "#E91E63"),
  save_dir      = "CellCycle_Results"
)
```

---

### 💡 Notes on Scientific Visualization
All generated plots across these modules are automatically exported as **300 DPI or 600 DPI `.tiff` files** using LZW compression, ensuring they fully comply with the strict formatting guidelines required by high-impact scientific journals.
